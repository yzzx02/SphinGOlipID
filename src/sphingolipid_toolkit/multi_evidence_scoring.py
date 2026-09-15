"""Independent evidence gates and confidence scores over unchanged fragment masses."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
from math import ceil, isfinite
import re

import numpy as np
import pandas as pd

from .evidence import EVIDENCE_PRIORITY, POOL_PRIORITY, annotate_theoretical_fragments, scoring_pool_for_evidence
from .evidence_config import EvidenceScoringConfig, POOL_WEIGHTS, SUPPORT_SATURATION_COUNT, get_class_rule
from .scoring import match_fragments


class AnnotationLevel(str, Enum):
    CANDIDATE = "candidate"
    SPECIES = "species"
    MOLECULAR_SPECIES = "molecular_species"

    @property
    def display_name(self):
        return {self.CANDIDATE:"MS1 candidate", self.SPECIES:"Species level",
            self.MOLECULAR_SPECIES:"Molecular species level"}[self]


@dataclass(frozen=True)
class PrecursorEvidence:
    precursor_pass: bool
    precursor_error_ppm: float | None
    precursor_adduct: str | None
    source: str = "existing_ms1_match"
    observed_mz: float | None = None

    def __post_init__(self):
        if not isinstance(self.precursor_pass, (bool, np.bool_)):
            raise ValueError("precursor_pass must be an explicit boolean")
        if self.observed_mz is not None and (not isfinite(self.observed_mz) or self.observed_mz <= 0):
            raise ValueError("Observed precursor m/z must be finite and positive")
        if self.precursor_pass and (self.precursor_error_ppm is None or
            not isfinite(self.precursor_error_ppm) or not self.precursor_adduct):
            raise ValueError("Passing precursor evidence requires actual error ppm and adduct")

    @classmethod
    def from_existing_match(cls, passed, observed_mz, theoretical_mz, adduct):
        """Record the existing MS1 decision; introduce no new mass tolerance."""
        if not isfinite(observed_mz) or not isfinite(theoretical_mz) or theoretical_mz <= 0:
            raise ValueError("Finite precursor masses and positive theoretical mass required")
        return cls(passed, (observed_mz-theoretical_mz)/theoretical_mz*1e6, adduct, observed_mz=observed_mz)


@dataclass
class CandidateEvidenceResult:
    summary: dict
    fragments: pd.DataFrame


def _canonical_evidence(table, rule):
    """Collapse metadata-bearing duplicate masses before matching and denominators."""
    records = []
    for _, row in table.iterrows():
        alternatives = row.get("evidence_annotations")
        if isinstance(alternatives, str) and alternatives:
            records.extend({**row.to_dict(), **entry} for entry in json.loads(alternatives))
        else:
            records.append(row.to_dict())
    data = pd.DataFrame(records) if records else table.copy()
    # Stable migration for prior metadata tables; no chemistry inference added.
    nl = data.fragment_type.eq("NL")
    old_role = data.evidence_role.eq("structure_informative")
    data.loc[nl & old_role,"evidence_role"] = data.loc[nl & old_role,"diagnostic_strength"].map(
        lambda value: "diagnostic_nl" if value == "diagnostic" else "supporting_nl")
    data["scoring_pool"] = [scoring_pool_for_evidence(row.fragment_type,row.evidence_role,rule) for row in data.itertuples()]
    data.loc[data.fragment_type.isin(["common","Precursor"]) | (nl & data.evidence_role.ne("diagnostic_nl")),"gate_eligible"] = False
    if not data.empty and "anno" in data and data.anno.nunique(dropna=False) > 1:
        raise ValueError("Evaluate one candidate and one spectrum at a time")
    if not data.fragment_type.isin(EVIDENCE_PRIORITY).all():
        raise ValueError("Only Precursor/HG/LCB/NL/common evidence categories are supported")
    for column in ("scoring_eligible", "gate_eligible"):
        if not data[column].map(lambda value: isinstance(value,(bool,np.bool_))).all():
            raise ValueError(f"{column} must contain booleans, not strings")
    masses = pd.to_numeric(data.theoretical_mz, errors="raise")
    if not np.isfinite(masses).all() or (masses <= 0).any():
        raise ValueError("Theoretical fragment masses must be finite and positive")
    data["theoretical_mz"] = masses
    rows = []
    for _, group in data.groupby("theoretical_mz", sort=True):
        selected = group.assign(_priority=group.scoring_pool.map(POOL_PRIORITY).fillna(3),
            _chemical_priority=group.fragment_type.map(EVIDENCE_PRIORITY)).sort_values(
            ["_priority","scoring_eligible","_chemical_priority","fragment_name"],ascending=[True,False,True,True],kind="stable").iloc[0].drop(["_priority","_chemical_priority"]).to_dict()
        selected["fragment_name"] = " | ".join(sorted({label for names in group.fragment_name for label in names.split(" | ")}))
        rows.append(selected)
    return pd.DataFrame(rows, columns=data.columns)


def evaluate_annotation_gates(rule, group_counts, precursor, config):
    """Each required group must cover ceil(fraction * eligible group size)."""
    passes, requirements = {}, {}
    for group in ("HG","LCB","NL"):
        matched, theoretical = group_counts.get(group,(0,0))
        required = ceil(config.minimum_group_fraction * theoretical)
        passes[group] = bool(precursor.precursor_pass and theoretical > 0 and matched >= required)
        requirements[group] = {"matched":matched,"theoretical":theoretical,"required":required}
    species = bool(precursor.precursor_pass and rule.species_groups and all(passes[g] for g in rule.species_groups))
    molecular = bool(species and rule.molecular_groups and all(passes[g] for g in rule.molecular_groups))
    level = AnnotationLevel.MOLECULAR_SPECIES if molecular else AnnotationLevel.SPECIES if species else AnnotationLevel.CANDIDATE
    return level, passes, requirements


def _fa_mode(annotation, total_composition, rule, lcb_supported):
    if rule.chain_count == 1:
        return "not_applicable_single_chain"
    if rule.chain_count != 2:
        return "UNSPECIFIED_multichain"
    if not lcb_supported:
        return "not_assigned"
    candidate = re.search(r"\(([mdtq])(\d+):(\d+)/(\d+):(\d+)\)", annotation)
    total = re.search(r"\(([mdtq])(\d+):(\d+)\)", str(total_composition))
    if not candidate or not total:
        return "UNSPECIFIED_total_or_chain_composition"
    hydroxyl, carbon, unsaturation, fa_carbon, fa_unsaturation = candidate.groups()
    if hydroxyl != total[1] or int(carbon)+int(fa_carbon) != int(total[2]) or int(unsaturation)+int(fa_unsaturation) != int(total[3]):
        return "UNSPECIFIED_inconsistent_composition"
    return "inferred_from_total_composition_minus_lcb"


def structural_quality(intensity, half_saturation):
    """LipidGate saturation quality, with the same [0,1] clipping."""
    if not isfinite(intensity) or not isfinite(half_saturation) or half_saturation <= 0:
        raise ValueError("Finite intensity and positive k required")
    intensity = min(max(float(intensity),0.),1.)
    return min((1.+half_saturation)*intensity/(intensity+half_saturation),1.)


def support_count_score(total, matched):
    if total < 0 or matched < 0 or matched > total:
        raise ValueError("Support counts require 0 <= matched <= theory")
    target = min(total,SUPPORT_SATURATION_COUNT)
    return (min(matched/target,1.) if target else None), target


def _precursor_mz(precursor, theory, explicit):
    value = explicit if explicit is not None else precursor.observed_mz
    if value is None and "target" in theory:
        targets = theory.target.dropna().unique()
        if len(targets) > 1:
            raise ValueError("Expected a single precursor")
        value = float(targets[0]) if len(targets) else None
    if value is not None and (not isfinite(value) or value <= 0):
        raise ValueError("Precursor m/z must be finite and positive")
    return value


def evaluate_candidate_evidence(
    lipid_class, annotation, theoretical_fragments, observed_fragments,
    precursor: PrecursorEvidence, *, ppm_tolerance=20., min_fragment_intensity=0.,
    config=EvidenceScoringConfig(), registry=None, origins=None,
    eligibility_overrides=None, total_composition=None, legacy_match_score=None,
    legacy_total_score=None, precursor_mz=None,
) -> CandidateEvidenceResult:
    """One spectrum/candidate; fixed structural 60/20 + count-based support 20.

    Score uses max structural quality, while gates independently retain half
    coverage. Supply the complete spectrum and actual precursor m/z for the
    LipidGate non-precursor normalization. MS1 decisions are never inferred.
    """
    if not isfinite(ppm_tolerance) or ppm_tolerance < 0 or not isfinite(min_fragment_intensity) or min_fragment_intensity < 0:
        raise ValueError("Matching tolerances and intensity cutoffs must be finite and nonnegative")
    rule = get_class_rule(lipid_class, registry)
    specified = rule.scoring_policy != "UNSPECIFIED"
    summary = {"lipid_class":rule.lipid_class,"annotation":annotation,
        "scoring_policy":rule.scoring_policy,"precursor_pass":bool(precursor.precursor_pass),
        "precursor_error_ppm":precursor.precursor_error_ppm,"precursor_adduct":precursor.precursor_adduct,
        "precursor_source":precursor.source,"hg_pass":False,"lcb_pass":False,"nl_pass":False,
        "S_HG":None,"S_LCB":None,"S_NL":None,"S_Supp":None,
        "S_primary":None,"S_secondary":None,"S_supp":None,
        "multi_evidence_score":None,"new_score":None,"confidence_score":None,
        "primary_evidence_type":rule.primary_evidence,"secondary_evidence_type":rule.secondary_evidence,
        "support_status":"not_evaluated","N_support_total":None,"N_support_matched":None,"support_target":None,
        "legacy_match_score":legacy_match_score,"legacy_total_score":legacy_total_score,
        "annotation_level":"candidate","annotation_level_display":"MS1 candidate",
        "lcb_directly_supported":False,"fa_directly_fragment_confirmed":False,
        "fa_assignment_mode":_fa_mode(annotation,total_composition,rule,False),
        "evidence_summary":"not evaluated: precursor gate failed", "gate_summary":"{}",
        "matched_fragment_names":"", "score_status":"precursor_failed", "rank_eligible":False,
        "class_rule_status":rule.status,"policy_notes":rule.notes}
    for pool in ("primary","secondary","support"):
        summary[pool+"_weight"] = POOL_WEIGHTS[pool] if specified else None
    for pool in ("primary","secondary"):
        summary.update({pool+"_gate_pass":False if specified else None,
            pool+"_matched_count":None,pool+"_theoretical_count":None})
    if not precursor.precursor_pass:
        return CandidateEvidenceResult(summary, pd.DataFrame())
    theory = theoretical_fragments.copy()
    if "fragment_type" not in theory:
        theory = annotate_theoretical_fragments(theory, lipid_class, origins, eligibility_overrides, registry=registry)
    theory = _canonical_evidence(theory,rule)
    theory = theory[theory.fragment_type.ne("Precursor")].copy()
    selected_precursor = _precursor_mz(precursor,theory,precursor_mz)
    observed = observed_fragments[["fragment_mz","fragment_intensity"]].copy().astype(float)
    numeric = observed.to_numpy()
    if not np.isfinite(numeric).all() or (numeric[:,0] <= 0).any() or (numeric[:,1] < 0).any():
        raise ValueError("Observed masses must be positive and intensities finite/nonnegative")
    base_peak = float(observed.fragment_intensity.max()) if len(observed) else 0.
    outside_cluster = observed.fragment_mz.sub(selected_precursor).abs().gt(config.precursor_cluster_exclusion_da) if selected_precursor is not None else pd.Series(False,index=observed.index)
    non_precursor_base = float(observed.loc[outside_cluster,"fragment_intensity"].max()) if outside_cluster.any() else 0.
    normalization_status = "non_precursor_base_peak" if non_precursor_base > 0 else "spectrum_relative_fallback" if selected_precursor is not None else "missing_precursor_mz_spectrum_relative_fallback"
    matches = match_fragments(observed,theory,ppm_tolerance,min_fragment_intensity)
    intensities = dict(zip(matches.get("theoretical_mz",[]),matches.get("observed_intensity",[])))
    observed_masses = dict(zip(matches.get("theoretical_mz",[]),matches.get("observed_mz",[])))
    theory["matched"] = theory.theoretical_mz.isin(intensities)
    theory["observed_mz"] = theory.theoretical_mz.map(observed_masses)
    theory["observed_intensity"] = theory.theoretical_mz.map(intensities).fillna(0.)
    theory["normalized_intensity"] = theory.observed_intensity / base_peak if base_peak > 0 else 0.
    theory["normalized_structural_intensity"] = theory.normalized_intensity
    # Like LipidGate: only peaks outside the isotope window get the override;
    # near-window non-precursor matches retain their spectrum-relative value.
    if non_precursor_base > 0:
        outside = theory.observed_mz.sub(selected_precursor).abs().gt(config.precursor_cluster_exclusion_da)
        theory.loc[outside,"normalized_structural_intensity"] = (theory.loc[outside,"observed_intensity"] / non_precursor_base).clip(0.,1.)
    # A centroid exactly at the selected precursor is not independent MS/MS
    # evidence, even if an erroneous legacy label calls it a structural ion.
    is_precursor = theory.observed_mz.eq(selected_precursor) if selected_precursor is not None else pd.Series(False,index=theory.index)
    theory["evidence_matched"] = theory.matched & theory.observed_intensity.gt(0) & ~is_precursor
    candidate_lcb = re.search(r"\(([mdtq]\d+:\d+)(?:/|\))", annotation)
    theory["gate_qualifies"] = theory.gate_eligible & theory.scoring_eligible
    lcb_rows = theory.fragment_type.eq("LCB")
    if candidate_lcb:
        theory.loc[lcb_rows,"gate_qualifies"] &= theory.loc[lcb_rows,"lcb_identity"].eq(candidate_lcb[1])
    else:
        theory.loc[lcb_rows,"gate_qualifies"] = False
    # Preserve the existing gate's spectrum-relative intensity minimum. The
    # non-precursor override belongs only to quality scoring, not coverage.
    positive = theory.evidence_matched & theory.normalized_intensity.ge(config.minimum_normalized_intensity)
    counts, evidence_text = {}, []
    for category in ("HG","LCB","NL","common"):
        eligible = theory.fragment_type.eq(category) & theory.scoring_eligible
        evidence_text.append(f"{category}: {int((eligible & theory.evidence_matched).sum())}/{int(eligible.sum())}")
        gate_theory = eligible & theory.gate_eligible
        counts[category] = (int((gate_theory & theory.gate_qualifies & positive).sum()),int(gate_theory.sum()))
    level, passes, requirements = evaluate_annotation_gates(rule,counts,precursor,config)
    theory["q"] = pd.Series(float("nan"), index=theory.index)
    for pool, evidence, k in (("primary",rule.primary_evidence,config.primary_half_intensity),
        ("secondary",rule.secondary_evidence,config.secondary_half_intensity)):
        mask = theory.scoring_pool.eq(pool) & theory.scoring_eligible
        matched = mask & theory.evidence_matched
        theory.loc[mask,"q"] = 0.
        theory.loc[matched,"q"] = [structural_quality(value,k) for value in theory.loc[matched,"normalized_structural_intensity"]]
        size = int(mask.sum())
        value = float(theory.loc[mask,"q"].max()) if size else None
        summary["S_"+pool] = value
        summary[pool+"_matched_count"] = int(matched.sum()) if evidence else None
        summary[pool+"_theoretical_count"] = size if evidence else None
        category = "NL" if evidence == "diagnostic_nl" else evidence
        summary[pool+"_gate_pass"] = passes[category] if category else None
        if category:
            summary["S_"+category] = value
    support_theory = theory.scoring_pool.eq("support") & theory.scoring_eligible
    support_total = int(support_theory.sum())
    support_matched = int((support_theory & theory.evidence_matched).sum())
    support_quality,target = support_count_score(support_total,support_matched)
    summary.update({"S_supp":support_quality,"S_Supp":support_quality,"N_support_total":support_total,
        "N_support_matched":support_matched,"support_target":target,
        "support_status":"available" if support_total else "no_support_theory"})
    missing = [pool for pool in ("primary","secondary") if summary["S_"+pool] is None]
    weighted = (POOL_WEIGHTS["primary"]*summary["S_primary"] + POOL_WEIGHTS["secondary"]*summary["S_secondary"] +
        POOL_WEIGHTS["support"]*(support_quality if support_quality is not None else 0.)) if specified and not missing else None
    summary.update({"annotation_level":level.value,"annotation_level_display":level.display_name,
        "hg_pass":passes["HG"],"lcb_pass":passes["LCB"],"nl_pass":passes["NL"],
        "multi_evidence_score":weighted,"new_score":weighted,"confidence_score":weighted if level != AnnotationLevel.CANDIDATE else None,
        "score_status":"available" if weighted is not None else "UNSPECIFIED_policy" if not specified else "missing_structural_theory:"+",".join(missing),
        "evidence_summary":"; ".join(evidence_text),"gate_summary":json.dumps(requirements,sort_keys=True),
        "matched_fragment_names":"; ".join(theory.loc[theory.matched,"fragment_name"]),
        "lcb_directly_supported":passes["LCB"],
        "fa_assignment_mode":_fa_mode(annotation,total_composition,rule,passes["LCB"]),
        "rank_eligible":level != AnnotationLevel.CANDIDATE and weighted is not None,
        "base_peak_intensity":base_peak,"non_precursor_base_peak_intensity":non_precursor_base,
        "structural_normalization_status":normalization_status,"precursor_mz":selected_precursor,
        "precursor_cluster_exclusion_da":config.precursor_cluster_exclusion_da,
        "minimum_group_fraction":config.minimum_group_fraction,"minimum_normalized_intensity":config.minimum_normalized_intensity})
    return CandidateEvidenceResult(summary,theory)


def add_rt_validation(result: CandidateEvidenceResult, status, *, error=None):
    """Attach orthogonal RT validation; never modify gates, score or level."""
    return CandidateEvidenceResult({**result.summary,"rt_validation_status":status,"rt_error":error},result.fragments.copy())


def rank_evidence_candidates(summaries):
    """Rank only gated, scoreable candidates without altering legacy order/score."""
    data = pd.DataFrame(summaries).copy()
    if data.empty:
        return data
    data["multi_evidence_rank"] = pd.Series(pd.NA,index=data.index,dtype="Int64")
    selected = data[data.rank_eligible].sort_values(["confidence_score","annotation"],ascending=[False,True],kind="stable")
    for rank,index in enumerate(selected.index,1):
        data.loc[index,"multi_evidence_rank"] = rank
    return data
