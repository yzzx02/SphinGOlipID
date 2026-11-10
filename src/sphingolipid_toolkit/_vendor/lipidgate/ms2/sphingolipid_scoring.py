"""Sphingolipid evidence evaluation, returning the shared CandidateScore contract."""
from __future__ import annotations

import re
from typing import Sequence
from .config import DEFAULT_SEARCH_CONFIG
from .models import CandidateScore, ExperimentalSpectrum, FragmentRecord, FragmentMatch, LibraryRecord
from .negative_gm3 import is_negative_gm3, has_negative_gm3_chain_evidence
from .esterified_ceramide import is_esterified_ceramide
from .positive_cer_sodium import is_positive_cer_sodium
from .positive_cer1p import is_positive_cer1p
from .positive_complex_ceramide import is_positive_complex_ceramide
from .positive_hexosylceramide import (
    has_chain_evidence as has_positive_lcb_chain_evidence,
    is_positive_hexosylceramide,
)
from .rules import DEFAULT_RULES, RuleSet
from .sphingolipid_rules import SphingoRule
from .gate_policy import evaluate_sphingolipid_gate, fah_only_low_confidence_gate_passes as _fah_only_low_confidence_gate_passes
from .scoring_policy import FAH_ONLY_FALLBACK_REASON, HG_ONLY_FALLBACK_REASON
from .record_facts import with_record_facts
from .scoring import (
    _empty_pool_scores, _calculate_pool_scores, _non_precursor_quality_overrides,
    _total_score_from_pool_scores, _with_negative_cer_match_bonus, _matched_intensity_sums,

)

_PRODUCTION_QUALITY_OVERRIDES = _non_precursor_quality_overrides


@with_record_facts
def score_sphingolipid_candidate(
    spectrum: ExperimentalSpectrum, record: LibraryRecord, *,
    matches: Sequence[FragmentMatch], rule: SphingoRule, series: str,
    rules: RuleSet = DEFAULT_RULES,
    precursor_tolerance_da: float | None = None,
    precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
    hg_only_classes: set[str], hg_only_min_score: float,
    hg_only_min_relative_intensity: float,
) -> CandidateScore:
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    precursor_out_of_tolerance = (
        abs(spectrum.precursor_mz - record.precursor_mz) > precursor_tolerance_da
        if precursor_tolerance_da is not None
        else abs(ppm_error) > precursor_tolerance_ppm
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["sphingo_rule"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )

    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["sphingo_rule"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )

    gate = evaluate_sphingolipid_gate(record, matches, rule, series)
    library_types, matched_types = gate.library_types, gate.matched_types
    adduct, class_key = gate.adduct, gate.class_key
    uses_negative_hg_fah_tiered_gate = gate.uses_negative_hg_fah_tiered_gate

    matched_intensity_sum, matched_relative_intensity_sum = _matched_intensity_sums(spectrum, matches)
    pool_scores = _calculate_pool_scores(
        matches,
        record,
        rules.get(record.compound_class),
        quality_relative_intensity_overrides=(None if _non_precursor_quality_overrides is _PRODUCTION_QUALITY_OVERRIDES
                                             else _non_precursor_quality_overrides(spectrum, record, matches)),
        spectrum=spectrum,
    )
    total_score = _total_score_from_pool_scores(pool_scores)
    total_score = _with_negative_cer_match_bonus(record, matches, total_score)
    standard_gate_passed = gate.passed
    has_library_hg_or_structural = bool(library_types & {"Diagnostic_HG", "C类碎片"})
    matched_hg_or_structural = bool(matched_types & {"Diagnostic_HG", "C类碎片"})
    hg_only_low_confidence_gate = (
        not standard_gate_passed
        and rule.allow_hg_only_fallback
        and re.sub(r"[^A-Za-z0-9]+", "", str(record.compound_class or "").upper()) in hg_only_classes
        and not uses_negative_hg_fah_tiered_gate
        and "Diagnostic_HG" in library_types
        and "Diagnostic_HG" in matched_types
        and "LCB碎片" not in matched_types
        and len(matches) >= 2
        and pool_scores["hg"].pool_score >= hg_only_min_score
        and any(
            match.fragment.fragment_type == "Diagnostic_HG"
            and match.experimental_peak.relative_intensity >= hg_only_min_relative_intensity
            for match in matches
        )
    )
    fah_only_low_confidence_gate = (
        not standard_gate_passed
        and rule.allow_fah_only_fallback
        and not uses_negative_hg_fah_tiered_gate
        and has_library_hg_or_structural
        and not matched_hg_or_structural
        and _fah_only_low_confidence_gate_passes(
            record,
            pool_scores,
            missing_groups=["sphingo_rule"],
            has_hg_or_structural_pool=has_library_hg_or_structural,
        )
    )
    passed = standard_gate_passed or fah_only_low_confidence_gate or hg_only_low_confidence_gate
    has_negative_chain_evidence = gate.has_negative_chain_evidence
    if standard_gate_passed and is_positive_complex_ceramide(record.compound_class,record.adduct):
        confirmed=has_positive_lcb_chain_evidence(record,matches)
        resolution_level='chain_level' if confirmed else 'species_level'
        downgrade_reason='' if confirmed else 'missing_lcb_chain_evidence'
    elif standard_gate_passed and is_positive_cer1p(record):
        confirmed = (adduct == "[M+H]+" and has_positive_lcb_chain_evidence(record,matches))
        resolution_level = "chain_level" if confirmed else "species_level"
        downgrade_reason = ("" if confirmed else "missing_lcb_chain_evidence" if adduct == "[M+H]+"
                            else "sum_composition_identity")
    elif standard_gate_passed and is_negative_gm3(record.compound_class, adduct):
        if has_negative_gm3_chain_evidence(record, matches):
            resolution_level = "chain_level"
            downgrade_reason = ""
        else:
            resolution_level = "species_level"
            downgrade_reason = "missing_lcb_chain_evidence"
    elif (standard_gate_passed
          and (is_positive_hexosylceramide(record.compound_class, adduct)
               or (class_key == "SM" and adduct == "[M+H]+")
               or (class_key == "CER" and adduct == "[M+H]+" and series in {"d", "t"}
                   and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name)))
          and not has_positive_lcb_chain_evidence(record, matches)):
        resolution_level = "species_level"
        downgrade_reason = "missing_lcb_chain_evidence"
    elif standard_gate_passed and uses_negative_hg_fah_tiered_gate and not has_negative_chain_evidence:
        resolution_level = "species_level"
        downgrade_reason = "missing_fah_chain_evidence"
    elif standard_gate_passed and is_positive_cer_sodium(record):
        resolution_level = "species_level"
        downgrade_reason = "sum_composition_identity"
    elif standard_gate_passed and class_key == "SM" and adduct == "[M+Na]+":
        resolution_level = "species_level"
        downgrade_reason = "headgroup_losses_only"
    elif standard_gate_passed and class_key not in {"SPB", "LSM"} and re.fullmatch(
        rf"{re.escape(record.compound_class)}\([mdt]?\d+:\d+\)", record.lipid_chain_name
    ):
        resolution_level = "species_level"
        downgrade_reason = "sum_composition_identity"
    elif standard_gate_passed:
        resolution_level = "chain_level"
        downgrade_reason = ""
    elif fah_only_low_confidence_gate:
        resolution_level = "tentative_chain_level"
        downgrade_reason = FAH_ONLY_FALLBACK_REASON
    elif hg_only_low_confidence_gate:
        resolution_level = "tentative_species_level"
        downgrade_reason = HG_ONLY_FALLBACK_REASON
    else:
        resolution_level = "class_level"
        downgrade_reason = "sphingo_rule_failed"

    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=passed,
        missing_required_groups=[] if passed else ["sphingo_rule"],
        ppm_error=ppm_error,
        resolution_level=resolution_level,
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason=downgrade_reason,
    )

