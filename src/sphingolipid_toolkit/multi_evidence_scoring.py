"""Independent evidence gates and confidence scores over unchanged fragment masses."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
from math import ceil, isfinite
import re

import numpy as np
import pandas as pd

from .evidence import EVIDENCE_PRIORITY, annotate_theoretical_fragments
from .evidence_config import EvidenceScoringConfig, SCORING_TEMPLATES, get_class_rule
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

    def __post_init__(self):
        if not isinstance(self.precursor_pass, (bool, np.bool_)):
            raise ValueError("precursor_pass must be an explicit boolean")
        if self.precursor_pass and (self.precursor_error_ppm is None or
            not isfinite(self.precursor_error_ppm) or not self.precursor_adduct):
            raise ValueError("Passing precursor evidence requires actual error ppm and adduct")

    @classmethod
    def from_existing_match(cls, passed, observed_mz, theoretical_mz, adduct):
        """Record the existing MS1 decision; introduce no new mass tolerance."""
        if not isfinite(observed_mz) or not isfinite(theoretical_mz) or theoretical_mz <= 0:
            raise ValueError("Finite precursor masses and positive theoretical mass required")
        return cls(passed, (observed_mz-theoretical_mz)/theoretical_mz*1e6, adduct)


@dataclass
class CandidateEvidenceResult:
    summary: dict
    fragments: pd.DataFrame


def _canonical_evidence(table):
    """Collapse metadata-bearing duplicate masses before matching and denominators."""
    data = table.copy()
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
        selected = group.assign(_priority=group.fragment_type.map(EVIDENCE_PRIORITY)).sort_values(
            ["_priority","scoring_eligible","fragment_name"],ascending=[True,False,True],kind="stable").iloc[0].drop("_priority").to_dict()
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


def evaluate_candidate_evidence(
    lipid_class, annotation, theoretical_fragments, observed_fragments,
    precursor: PrecursorEvidence, *, ppm_tolerance=20., min_fragment_intensity=0.,
    config=EvidenceScoringConfig(), registry=None, templates=None, origins=None,
    eligibility_overrides=None, total_composition=None, legacy_match_score=None,
    legacy_total_score=None,
) -> CandidateEvidenceResult:
    """Evaluate a single candidate in one complete MS/MS spectrum.

    Pass the complete observed spectrum, including unmatched peaks, for base
    peak normalization. The optional raw cutoff is the caller's existing matcher
    cutoff; the evidence gate intensity is separately configured in [0,1].
    Failed precursor evidence skips MS/MS evaluation entirely.
    """
    if not isfinite(ppm_tolerance) or ppm_tolerance < 0 or not isfinite(min_fragment_intensity) or min_fragment_intensity < 0:
        raise ValueError("Matching tolerances and intensity cutoffs must be finite and nonnegative")
    rule = get_class_rule(lipid_class, registry)
    summary = {"lipid_class":rule.lipid_class,"annotation":annotation,
        "template":rule.scoring_template or "UNSPECIFIED", "precursor_pass":bool(precursor.precursor_pass),
        "precursor_error_ppm":precursor.precursor_error_ppm,"precursor_adduct":precursor.precursor_adduct,
        "precursor_source":precursor.source,"hg_pass":False,"lcb_pass":False,"nl_pass":False,
        "S_HG":0.,"S_LCB":0.,"S_NL":0.,"S_common":0.,"multi_evidence_score":None,"confidence_score":None,
        "legacy_match_score":legacy_match_score,"legacy_total_score":legacy_total_score,
        "annotation_level":"candidate","annotation_level_display":"MS1 candidate",
        "lcb_directly_supported":False,"fa_directly_fragment_confirmed":False,
        "fa_assignment_mode":_fa_mode(annotation,total_composition,rule,False),
        "evidence_summary":"not evaluated: precursor gate failed", "gate_summary":"{}",
        "matched_fragment_names":"", "score_status":"precursor_failed", "rank_eligible":False,
        "class_rule_status":rule.status,"policy_notes":rule.notes}
    if not precursor.precursor_pass:
        return CandidateEvidenceResult(summary, pd.DataFrame())
    theory = theoretical_fragments.copy()
    if "fragment_type" not in theory:
        theory = annotate_theoretical_fragments(theory, lipid_class, origins, eligibility_overrides, registry=registry)
    theory = _canonical_evidence(theory)
    # MS1 evidence is never an MS/MS node, denominator, matched intensity or score.
    theory = theory[theory.fragment_type.ne("Precursor")].copy()
    observed = observed_fragments[["fragment_mz","fragment_intensity"]].copy()
    numeric = observed.to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or (numeric[:,0] <= 0).any() or (numeric[:,1] < 0).any():
        raise ValueError("Observed masses must be positive and intensities finite/nonnegative")
    observed = observed.astype(float)
    base_peak = float(observed.fragment_intensity.max()) if len(observed) else 0.
    matches = match_fragments(observed,theory,ppm_tolerance,min_fragment_intensity)
    matched_intensity = dict(zip(matches.get("theoretical_mz",[]),matches.get("observed_intensity",[])))
    theory["matched"] = theory.theoretical_mz.isin(matched_intensity)
    theory["normalized_intensity"] = theory.theoretical_mz.map(matched_intensity).fillna(0.) / base_peak if base_peak > 0 else 0.
    theory["q"] = [float(row.normalized_intensity)/(float(row.normalized_intensity)+config.k_for(row.fragment_type,row.diagnostic_strength))
        if row.scoring_eligible else 0. for row in theory.itertuples()]
    # A chain-specific label must refer to the candidate's stated LCB; an
    # unrelated LCB annotation cannot satisfy this candidate's molecular gate.
    candidate_lcb = re.search(r"\(([mdtq]\d+:\d+)(?:/|\))", annotation)
    theory["gate_qualifies"] = theory.gate_eligible & theory.scoring_eligible
    lcb_rows = theory.fragment_type.eq("LCB")
    if candidate_lcb:
        theory.loc[lcb_rows,"gate_qualifies"] &= theory.loc[lcb_rows,"lcb_identity"].eq(candidate_lcb[1])
    else:
        theory.loc[lcb_rows,"gate_qualifies"] = False
    positive = theory.normalized_intensity.gt(0) & theory.normalized_intensity.ge(config.minimum_normalized_intensity)
    counts, evidence_text = {}, []
    for category in ("HG","LCB","NL","common"):
        eligible = theory[theory.fragment_type.eq(category) & theory.scoring_eligible]
        size = len(eligible)
        matched = int(eligible.matched.sum())
        summary["S_"+category] = float(eligible.q.sum()/size) if size else 0.
        summary["N_"+category] = size
        evidence_text.append(f"{category}: {matched}/{size}")
        gate_theory = theory.fragment_type.eq(category) & theory.gate_eligible & theory.scoring_eligible
        counts[category] = (int((gate_theory & theory.gate_qualifies & positive).sum()),int(gate_theory.sum()))
    level, passes, requirements = evaluate_annotation_gates(rule,counts,precursor,config)
    template = (templates if templates is not None else SCORING_TEMPLATES).get(rule.scoring_template)
    missing = [category for category,_ in template.weights if summary["N_"+category] == 0] if template else []
    # Do not silently score an inapplicable class template by filling absent
    # categories with zeros or redistributing its weights.
    weighted = sum(weight*summary["S_"+category] for category,weight in template.weights) if template and not missing else None
    summary.update({"annotation_level":level.value,"annotation_level_display":level.display_name,
        "hg_pass":passes["HG"],"lcb_pass":passes["LCB"],"nl_pass":passes["NL"],
        "multi_evidence_score":weighted,"confidence_score":weighted if level != AnnotationLevel.CANDIDATE else None,
        "score_status":"available" if weighted is not None else "UNSPECIFIED_template" if not template else "missing_theoretical_groups:"+",".join(missing),
        "evidence_summary":"; ".join(evidence_text),"gate_summary":json.dumps(requirements,sort_keys=True),
        "matched_fragment_names":"; ".join(theory.loc[theory.matched,"fragment_name"]),
        "lcb_directly_supported":passes["LCB"],
        "fa_assignment_mode":_fa_mode(annotation,total_composition,rule,passes["LCB"]),
        "rank_eligible":level != AnnotationLevel.CANDIDATE and weighted is not None,
        "base_peak_intensity":base_peak,"minimum_group_fraction":config.minimum_group_fraction,
        "minimum_normalized_intensity":config.minimum_normalized_intensity})
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
