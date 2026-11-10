"""Numeric evidence-pool scoring, with compatible exports for policy helpers."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
from . import native_bridge
from .positive_sphingolipid import is_positive_sphingolipid, OTHER_MIN_RELATIVE_INTENSITY
from typing import (
    Dict,
    List,
    Sequence,
    Tuple,
)
from .chain_utils import extract_chain_tokens as _extract_chain_tokens
from .esterified_ceramide import is_esterified_ceramide
from .positive_cer_sodium import is_positive_cer_sodium, POOL_WEIGHTS as CER_SODIUM_POOL_WEIGHTS
from .positive_cer1p import is_positive_cer1p, POOL_WEIGHTS as CER1P_POOL_WEIGHTS
from .evidence_policy import (
    CANDIDATE_HG_FRAGMENT_TYPE,
    CLASS_SPECIFIC_HG_MZ,
    FREE_SPHINGOID_BASE_CLASS_KEYS,
    LOSS_FRAGMENT_TYPES,
    NEGATIVE_CER_CHAIN_FRAGMENT_TYPES,
    NEGATIVE_FA_LOSS_SUPPORT_ONLY_CLASSES,
    POSITIVE_FA_FRAG_AS_LOSS_CLASSES,
    POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
    POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES,
    POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES,
    POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES,
    POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES,
    PRECURSOR_FRAGMENT_REQUIRED_CLASSES,
    SPB_UNASSIGNED_DIAGNOSTIC_NAMES,
    SPHINGOLIPID_LCB_FRAGMENT_TYPE,
    _chain_evidence_count_for_fragments,
    _chain_evidence_count_for_matches,
    _count_fah_expected_chains,
    _count_positive_nonzero_chains,
    _fa_frag_counts_as_effective_loss,
    _fragment_counts_as_chain_resolving_loss,
    _fragment_counts_as_effective_loss,
    _fragment_counts_as_fa_loss_gate,
    _fragment_counts_as_hg,
    _is_ahexcer_overlapping_lcb_loss,
    _is_cl_chain_info_fragment,
    _is_cl_double_negative_record,
    _is_class_specific_hg_fragment,
    _is_fa_record,
    _is_negative_ahexcer_record,
    _is_positive_adduct,
    _is_positive_glyceride_chain_gate_fragment,
    _is_positive_glyceride_rco_c3h6o2_fragment,
    _is_positive_glyceride_rco_fragment,
    _is_positive_glyceride_rco_gate_record,
    _is_positive_hg_dominant_phospholipid,
    _is_positive_mg_record,
    _is_positive_oxtg_full_fah_record,
    _is_rco_fragment_name,
    _matched_cl_chain_info_fragment_count,
    _matched_fah_tokens,
    _matched_hg_fragment_count,
    _matched_hg_max_relative_intensity,
    _matched_positive_glyceride_rco_fragment_count,
    _matched_positive_signature_fragment_count,
    _matched_precursor_ion_fragment_count,
    _matched_required_fah_count,
    _normalized_class_key,
    _normalized_compound_class,
    _pool_for_fragment,
    _pool_for_scoring_fragment,
    _positive_glyceride_chain_gate_group,
    _positive_required_loss_hits,
    _positive_single_chain_species_hg_only,
    _record_expected_fah_tokens,
    _record_fa_loss_fragment_count,
    _record_has_candidate_hg,
    _record_has_diagnostic_hg,
    _record_hg_fragment_count,
    _record_loss_fragment_count,
    _record_positive_glyceride_rco_fragment_count,
    _record_positive_signature_fragment_count,
    _record_precursor_ion_fragments,
    _requires_precursor_fragment,
    _scoring_pools_for_fragment,
    _uses_positive_support_pool_scoring,
)
from .gate_policy import (
    sodium_phospholipid_missing_hg,
    glyceride_substituent_evidence_group,
    glyceride_substituent_group_multiplicity,
    is_positive_three_substituent_glyceride,
    ETHER_GLYCERIDE_HALF_RCO_GATE_CLASSES,
    FAH_ONLY_FALLBACK_BLOCKED_CLASSES,
    POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY,
    POSITIVE_GLYCERIDE_RCO_MIN_HITS,
    STRICT_NEGATIVE_HG_CLASSES,
    cl_required_fah_hit_count as _cl_required_fah_hit_count,
    derive_required_groups as _derive_required_groups,
    fah_only_low_confidence_gate_passes as _fah_only_low_confidence_gate_passes,
    hg_only_low_confidence_gate_passes as _hg_only_low_confidence_gate_passes,
    missing_required_groups as _missing_required_groups,
    positive_glyceride_rco_gate_available as _positive_glyceride_rco_gate_available,
    positive_glyceride_rco_gate_passes as _positive_glyceride_rco_gate_passes,
    positive_glyceride_rco_required_hits as _positive_glyceride_rco_required_hits,
    positive_oxtg_full_fah_gate_passes as _positive_oxtg_full_fah_gate_passes,
    required_fah_gate_hits as _required_fah_gate_hits,
    required_fractional_hits as _required_fractional_hits,
    required_negative_hg_hits as _required_negative_hg_hits,
    required_positive_hg_hits as _required_positive_hg_hits,
)
from .matching import match_fragments
from .models import (
    CandidateScore,
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    PoolScore,
)
from .negative_gm3 import is_negative_gm3
from .negative_lpe import is_negative_lpe
from .negative_lpc import is_negative_lpc
from .positive_lpc_ether import is_positive_ether_lpc
from .negative_hexcer import is_negative_hexcer
from .oxidized_chain_evidence import oxidized_chain_requirements, oxidized_chain_missing_groups
from .positive_hexosylceramide import is_positive_ganglioside, is_positive_hexosylceramide
from .positive_catalog_glycosphingolipid import is_positive_catalog_glycosphingolipid
from .positive_catalog_phosphoceramide import is_positive_catalog_phosphoceramide
from .positive_catalog_legacy_sphingolipid import is_positive_catalog_legacy_sphingolipid,hg_minimum as catalog_legacy_hg_minimum
from .positive_ganglioside import is_positive_ganglioside as is_positive_extended_ganglioside
from .positive_pe_cer import is_positive_pe_cer
from .positive_complex_ceramide import is_positive_complex_ceramide
from .positive_spb import is_positive_spb
from .positive_pe_sodium import is_positive_pe_sodium
from .positive_pc_sodium import is_positive_pc_sodium
from .positive_single_chain_h import is_protonated_single_chain_nagly, SINGLE_CHAIN_POOL_WEIGHTS, single_chain_rule
from .record_facts import with_record_facts, cached_record_fact, cached_rule_fact, cached_spectrum_fact, get_native_query
from .resolution_policy import (
    CHAIN_LEVEL_INFO_MISSING_REASON,
    PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON,
    profile_chain_resolution,
    positive_phospholipid_chain_resolution,
)
from .rules import ClassRule
from .scoring_policy import FAH_ONLY_FALLBACK_REASON, HG_ONLY_FALLBACK_REASON


# Compatibility exports of helpers now owned by the policy modules.
from .gate_policy import (
    is_positive_ps as _is_positive_ps_record,
    is_positive_pep,
    positive_pep_chain_gate_counts,
    is_positive_ps_chain_ketene_loss as _is_positive_ps_headgroup_ketene_loss,
    is_positive_ps_headgroup_loss as _is_positive_ps_headgroup_loss,
    is_positive_tg_est_full_loss_record as _is_positive_tg_est_full_loss_gate_record,
    is_positive_tg_full_chain_record as _is_positive_tg_full_chain_gate_record,
    is_positive_tg_o_full_loss_record as _is_positive_tg_o_full_loss_gate_record,
    positive_tg_est_full_loss_gate_passes as _positive_tg_est_full_loss_gate_passes,
    positive_tg_full_chain_gate_passes as _positive_tg_full_chain_gate_passes,
    positive_tg_o_full_loss_gate_passes as _positive_tg_o_full_loss_gate_passes,
)

POOL_NAMES = ("fah", "hg", "other")
FRAGMENT_QUALITY_FULL_SCORE_RELATIVE_INTENSITY = 0.10
PRIMARY_POOL_SATURATION_HALF_INTENSITY = 0.10
SECONDARY_POOL_SATURATION_HALF_INTENSITY = 0.05
IDENTICAL_THREE_CHAIN_SECOND_EVIDENCE_MULTIPLIER = 0.50
# Exclude the precursor and M through M+4 when normalizing structural intensity.
PRECURSOR_CLUSTER_EXCLUSION_DA = 4.1
NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS = frozenset({"FA", "OXFA", "NAE"})
SCORE_MIN = 0.0
SCORE_MAX = 100.0
_NATIVE_EVIDENCE_PRECHECK = object()
POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES = {"ASM", "LPC", "LPCO", "SM", "LSM"}
POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS = {"fah": 20.0, "hg": 60.0, "other": 20.0}
AHEXCER_NEGATIVE_POOL_WEIGHTS = {"fah": 60.0, "hg": 20.0, "other": 20.0}
SM_SODIUM_POOL_WEIGHTS = {"fah": 0.0, "hg": 75.0, "other": 25.0}
LPE_SODIUM_POOL_WEIGHTS = {"fah": 0.0, "hg": 75.0, "other": 25.0}
ADGGA_POSITIVE_POOL_WEIGHTS = {"fah": 20.0, "hg": 60.0, "other": 20.0}
POSITIVE_HG_CHAIN_LEVEL_CLASSES = {
    "CL", "MLCL", "NAPE", "LNAPE",
}


def _empty_pool_scores() -> Dict[str, PoolScore]:
    return {
        pool_name: PoolScore(
            pool_name=pool_name,
            matched_count=0,
            total_count=0,
            count_ratio=0.0,
            intensity_ratio=0.0,
            weight_ratio=0.0,
            pool_score=0.0,
        )
        for pool_name in POOL_NAMES
    }


def _match_fragments(
    peaks: Sequence[ExperimentalPeak],
    fragments: Sequence[FragmentRecord],
    mz_tolerance: float | None,
    ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
) -> List[FragmentMatch]:
    return match_fragments(peaks, fragments, lambda mz: _fragment_window_da(mz, mz_tolerance, ppm_tolerance), experimental_mz)


def _fragment_window_da(fragment_mz: float, mz_tolerance: float | None, ppm_tolerance: float | None) -> float:
    if mz_tolerance is not None:
        return float(mz_tolerance)
    if ppm_tolerance is not None:
        return abs(float(fragment_mz)) * float(ppm_tolerance) * 1e-6
    return abs(float(fragment_mz)) * 10.0 * 1e-6


def _pool_weights_for_record(record: LibraryRecord, rule: ClassRule) -> Dict[str, float]:
    if (is_positive_catalog_phosphoceramide(record.compound_class, record.adduct)
            or is_positive_catalog_legacy_sphingolipid(record.compound_class,record.adduct)):
        return {"hg":60.0,"fah":0.0,"other":40.0}
    if is_positive_catalog_glycosphingolipid(record.compound_class, record.adduct):
        return {"hg": 50.0, "fah": 0.0, "b": 25.0, "other": 25.0}
    if is_positive_complex_ceramide(record.compound_class,record.adduct):
        return {'fah':0.,'hg':60.,'lcb':20.,'other':20.}
    if is_positive_cer1p(record) or (record.compound_class == "S1P" and record.adduct == "[M+H]+"):
        return CER1P_POOL_WEIGHTS
    if is_positive_cer_sodium(record):
        return CER_SODIUM_POOL_WEIGHTS
    if (record.compound_class == "Cer" and record.adduct == "[M+H]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name)):
        return {"fah": 0.0, "hg": 0.0, "lcb": 60.0, "other": 40.0}
    if is_positive_spb(record):
        return {"fah": 0.0, "hg": 0.0, "lcb": 75.0, "other": 25.0}
    if is_positive_ether_lpc(record):
        return {"fah": 20.0, "hg": 60.0, "other": 20.0}
    if is_negative_lpc(record):
        return {"fah": 0.0, "hg": 75.0, "other": 25.0}
    class_key = _normalized_class_key(record.compound_class)
    adduct = str(record.adduct or "").strip()
    if is_protonated_single_chain_nagly(record):
        return SINGLE_CHAIN_POOL_WEIGHTS
    if is_negative_hexcer(record.compound_class, adduct):
        return {"hg": 60.0, "fah": 20.0, "lcb": 20.0, "other": 0.0}
    if is_positive_ganglioside(record.compound_class, adduct):
        return {"hg": 50.0, "fah": 0.0, "b": 25.0, "other": 25.0}
    if record.compound_class == "Gb4" and adduct == "[M+H]+":
        return {"hg": 75.0, "fah": 0.0, "other": 25.0}
    if is_positive_hexosylceramide(record.compound_class, adduct):
        return {"hg": 60.0, "fah": 0.0, "lcb": 20.0, "other": 20.0}
    if is_negative_gm3(record.compound_class, adduct) or is_positive_pe_cer(record.compound_class, adduct):
        return {"hg": 60.0, "fah": 0.0, "lcb": 20.0, "other": 20.0}
    if is_esterified_ceramide(record.compound_class, record.lipid_chain_name):
        if _is_positive_adduct(adduct):
            return {"fah": 60.0, "hg": 0.0, "lcb": 20.0, "other": 20.0}
        return {"fah": 60.0, "hg": 20.0, "other": 20.0}
    if _is_negative_ahexcer_record(record):
        return AHEXCER_NEGATIVE_POOL_WEIGHTS
    if class_key == "SM" and adduct in {"[M+H]+", "[M+Na]+"}:
        return SM_SODIUM_POOL_WEIGHTS
    if class_key == "ADGGA" and _is_positive_adduct(record.adduct):
        return ADGGA_POSITIVE_POOL_WEIGHTS
    if class_key == "PEP" and adduct == "[M+Na]+":
        return POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS
    if class_key in {"LPE", "LPC"} and adduct == "[M+Na]+":
        return LPE_SODIUM_POOL_WEIGHTS
    if (
        adduct == "[M-H]-" and class_key in {"CER1P", "CERP"}
    ) or (
        adduct in {"[M+CH3COO]-", "[M+HCOO]-"} and class_key == "SM"
    ):
        return {"fah": 20.0, "hg": 60.0, "other": 20.0}
    if _is_positive_hg_dominant_phospholipid(record) or (
        _is_positive_adduct(record.adduct)
        and class_key in POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES
    ):
        return POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS
    return rule.score_profile.pool_weights


def _effective_pool_weights(
    pool_weights: Dict[str, float],
    active_pools: Sequence[str],
) -> Dict[str, float]:
    """Align the dominant structural weight with the available evidence.

    The score profiles describe a dominant structural pool (FAH or HG), a
    secondary structural pool and ordinary support.  Some lipid classes only
    have one structural pool in their library records.  In that case the sole
    structural pool is the dominant pool; leaving it at the secondary 20-point
    weight would incorrectly turn a structural/support score into 50/50.
    """

    effective = dict(pool_weights)
    active_structural_pools = [
        pool_name for pool_name in ("fah", "hg") if pool_name in active_pools
    ]
    if len(active_structural_pools) == 1:
        dominant_structural_weight = max(
            max(float(pool_weights.get("fah", 0.0)), 0.0),
            max(float(pool_weights.get("hg", 0.0)), 0.0),
        )
        effective[active_structural_pools[0]] = dominant_structural_weight
    return effective


def _positive_three_substituent_top_half_quality(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    quality_relative_intensity_overrides: Dict[int, float] | None,
    half_saturation: float,
) -> float:
    """Average the strongest half of three substituents using distinct peaks.

    TG, TG-O and TG-EST each use two evidence groups. Repeated substituents can
    share one physical diagnostic loss, while ether and FAHFA groups remain
    distinct even when their chain compositions equal an ordinary FA chain.
    """

    substituent_multiplicity = glyceride_substituent_group_multiplicity(record)
    expected_groups = frozenset(substituent_multiplicity)
    if not expected_groups:
        return 0.0
    strongest_quality_by_group: Dict[str, float] = {}
    for match in matches:
        group = glyceride_substituent_evidence_group(record, match.fragment)
        if group is None or group not in expected_groups:
            continue
        quality = _saturation_fragment_quality(
            _match_relative_intensity(
                match,
                quality_relative_intensity_overrides=quality_relative_intensity_overrides,
            ),
            half_saturation,
        )
        strongest_quality_by_group[group] = max(
            strongest_quality_by_group.get(group, 0.0),
            quality,
        )
    distinct_group_qualities = sorted(strongest_quality_by_group.values(), reverse=True)
    if len(expected_groups) == 1 and sum(substituent_multiplicity.values()) == 3:
        # One physical loss legitimately represents all three identical acyl
        # substituents, so the gate must not demand duplicate peaks.  It also
        # must not receive the same confidence as two independently observed
        # chains: use one full-quality observation and one half-quality virtual
        # replicate for the top-half (two-of-three) structural score.
        strongest_quality = distinct_group_qualities[0] if distinct_group_qualities else 0.0
        return (
            strongest_quality
            + strongest_quality * IDENTICAL_THREE_CHAIN_SECOND_EVIDENCE_MULTIPLIER
        ) / 2.0
    selected = distinct_group_qualities[:2]
    return sum(selected) / len(selected) if selected else 0.0


@cached_record_fact
def _library_fragments_by_pool(record):
    result = defaultdict(list)
    for fragment in record.fragments:
        for pool_name in _scoring_pools_for_fragment(record, fragment, matched=False):
            result[pool_name].append(fragment)
    return {name: tuple(fragments) for name, fragments in result.items()}


@cached_rule_fact
def _library_pool_configuration(record, rule):
    return _build_pool_configuration(record, rule, _library_fragments_by_pool(record))


def _build_pool_configuration(record, rule, total_fragments_by_pool):
    pool_weights = _pool_weights_for_record(record, rule)
    pool_names = POOL_NAMES + tuple(name for name in ("lcb", "b") if name in pool_weights)
    active_pools = [name for name in pool_names if total_fragments_by_pool.get(name)]
    effective_pool_weights = _effective_pool_weights(pool_weights, active_pools)
    active_sum = sum(max(float(effective_pool_weights.get(name, 0.0)), 0.0) for name in active_pools)
    hexosylceramide = is_positive_hexosylceramide(record.compound_class, record.adduct)
    if (is_negative_gm3(record.compound_class, record.adduct) or hexosylceramide or is_positive_spb(record)
            or (record.compound_class == "Cer" and record.adduct == "[M+H]+"
                and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name))):
        active_sum = sum(pool_weights.values())
    equal_weight = 100.0 / len(active_pools) if active_pools else 0.0
    structural = [name for name in ("fah", "hg", "lcb", "b") if total_fragments_by_pool.get(name)]
    primary = max(structural, key=lambda name: float(effective_pool_weights.get(name, 0.0))) if structural else None
    result = []
    for name in pool_names:
        total_count = len(total_fragments_by_pool.get(name, ()))
        weight = (max(float(effective_pool_weights.get(name, 0.0)), 0.0) / active_sum * 100.0
                  if total_count and active_sum > 0.0 else equal_weight if total_count else 0.0)
        half = PRIMARY_POOL_SATURATION_HALF_INTENSITY if name == primary else SECONDARY_POOL_SATURATION_HALF_INTENSITY
        kind = (4 if name == "other" and is_positive_sphingolipid(record.compound_class, record.adduct) else
                0 if name == "other" and not _is_negative_ahexcer_record(record) else
                2 if name == "hg" and hexosylceramide and record.compound_class not in {"HexCer", "Hex2Cer", "LacCer", "SHexCer", "SHexCer+O"}
                    and not is_positive_extended_ganglioside(record.compound_class,record.adduct)
                    and not is_positive_catalog_glycosphingolipid(record.compound_class,record.adduct)
                    and not is_positive_catalog_phosphoceramide(record.compound_class,record.adduct)
                    and not (is_positive_catalog_legacy_sphingolipid(record.compound_class,record.adduct)
                             and catalog_legacy_hg_minimum(record.compound_class)<3) else
                2 if name == "fah" and oxidized_chain_requirements(record) is not None else
                3 if name == "fah" and is_positive_three_substituent_glyceride(record) else 1)
        result.append((name, total_count, weight, half, kind))
    return tuple(result)


@cached_rule_fact
def _numeric_pool_plan(record, rule):
    configuration = _library_pool_configuration(record, rule)
    pool_bits = {name: 1 << i for i, (name, *_rest) in enumerate(configuration)}
    triple = any(item[4] == 3 for item in configuration)
    multiplicity = glyceride_substituent_group_multiplicity(record) if triple else {}
    groups = {name: i for i, name in enumerate(multiplicity)}
    routes = {}
    for fragment in record.fragments:
        mask = sum(pool_bits.get(name, 0) for name in _scoring_pools_for_fragment(record, fragment, matched=True))
        group = glyceride_substituent_evidence_group(record, fragment) if triple else None
        routes[id(fragment)] = (mask, groups.get(group, -1))
    return configuration, routes, (len(groups), len(groups) == 1 and sum(multiplicity.values()) == 3,
                                   IDENTICAL_THREE_CHAIN_SECOND_EVIDENCE_MULTIPLIER)


@cached_rule_fact
def _prepared_numeric_pool_plan(record, rule):
    kernel = native_bridge.load_bridge()
    if kernel is None or not hasattr(kernel, 'create_pool_plan'):
        return None
    configuration, routes, groups = _numeric_pool_plan(record, rule)
    capsule = kernel.create_pool_plan(routes, configuration, groups)
    # Include a conservative allocation allowance in the bounded fact LRU.
    return None if capsule is None else (capsule, configuration, bytes(1024 + len(routes)*64))


def _calculate_pool_scores(
    matches: Sequence[FragmentMatch],
    record: LibraryRecord,
    rule: ClassRule,
    fragments: Sequence[FragmentRecord] | None = None,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
    *, spectrum: ExperimentalSpectrum | None = None,
) -> Dict[str, PoolScore]:
    if _normalized_class_key(record.compound_class) == "NAE":
        # Presence is the complete evidence under the explicit NAE policy.
        # Shared-peak reranking must keep this binary definition as well.
        matched = any(
            match.fragment.fragment_type == "Precursor Ion"
            and abs(match.fragment.mz - record.precursor_mz) < 1e-8
            for match in matches
        )
        result = _empty_pool_scores()
        value = float(matched)
        result["fah"] = PoolScore("fah", int(matched), 1, value, value, 1.0, 100.0 * value)
        return result
    if (fragments is None and _saturation_fragment_quality is _PRODUCTION_SATURATION
            and _match_relative_intensity is _PRODUCTION_MATCH_INTENSITY
            and _positive_three_substituent_top_half_quality is _PRODUCTION_TRIPLE_QUALITY
            and _top_two_fragment_quality is _PRODUCTION_TOP_TWO):
        kernel = native_bridge.load_bridge()
        if kernel is not None and hasattr(kernel, "pool_metrics"):
            prepared = _prepared_numeric_pool_plan(record, rule)
            query = get_native_query(spectrum)
            normalization = (_native_pool_normalization(spectrum, record)
                             if spectrum is not None and quality_relative_intensity_overrides is None
                             and hasattr(kernel, 'pool_metrics_normalized') else None)
            if normalization is not None and prepared is not None:
                numeric = kernel.pool_metrics_normalized(matches, prepared[0], normalization)
                if numeric is not None:
                    result = {item[0]: PoolScore(item[0], values[0], item[1], values[1], values[2], 0.0, values[3])
                              for item, values in zip(prepared[1], numeric)}
                    return _lcb_other_support(result, matches, record, rule, fragments,
                                                  quality_relative_intensity_overrides, spectrum)
            if spectrum is not None and quality_relative_intensity_overrides is None:
                quality_relative_intensity_overrides = _non_precursor_quality_overrides(spectrum, record, matches)
            if prepared is not None:
                configuration = prepared[1]
                numeric = (kernel.pool_metrics_query(matches, prepared[0], quality_relative_intensity_overrides, query)
                           if query is not None and hasattr(kernel, 'pool_metrics_query') else
                           kernel.pool_metrics(matches, prepared[0], quality_relative_intensity_overrides))
            else:
                configuration, routes, groups = _numeric_pool_plan(record, rule)
                numeric = kernel.pool_metrics(matches, routes, configuration, groups, quality_relative_intensity_overrides)
            if numeric is not None:
                result = {item[0]: PoolScore(item[0], values[0], item[1], values[1], values[2], 0.0, values[3])
                          for item, values in zip(configuration, numeric)}
                return _lcb_other_support(result, matches, record, rule, fragments,
                                              quality_relative_intensity_overrides, spectrum)
    if spectrum is not None and quality_relative_intensity_overrides is None:
        quality_relative_intensity_overrides = _non_precursor_quality_overrides(spectrum, record, matches)
    matched_by_pool: Dict[str, List[FragmentMatch]] = defaultdict(list)
    total_fragments_by_pool: Dict[str, List[FragmentRecord]] = defaultdict(list)
    if fragments is None:
        # Pool membership is explicitly matched-independent. The alternate
        # fragment subset used for FA keeps the original occurrence path.
        total_fragments_by_pool = _library_fragments_by_pool(record)
    else:
        matched_fragment_ids = {id(match.fragment) for match in matches}
        for fragment in fragments:
            for pool_name in _scoring_pools_for_fragment(
                record, fragment, matched=id(fragment) in matched_fragment_ids,
            ):
                total_fragments_by_pool[pool_name].append(fragment)
    for match in matches:
        for pool_name in _scoring_pools_for_fragment(record, match.fragment, matched=True):
            matched_by_pool[pool_name].append(match)
    configuration = (_library_pool_configuration(record, rule) if fragments is None else
                     _build_pool_configuration(record, rule, total_fragments_by_pool))
    result: Dict[str, PoolScore] = {}
    for pool_name, total_count, dynamic_pool_weight, half_saturation, quality_kind in configuration:
        pool_matches = matched_by_pool.get(pool_name, [])
        if quality_kind == 4:
            # Admission and the audit retain every match. Only support credit
            # uses the strict 0.5% threshold against the original base peak.
            pool_matches = [match for match in pool_matches
                            if match.experimental_peak.relative_intensity > OTHER_MIN_RELATIVE_INTENSITY]
        matched_count = len(pool_matches)
        count_ratio = matched_count / total_count if total_count else 0.0
        if quality_kind in (0, 4):
            # Supporting fragments saturate after three observations.  When a
            # library only contains one or two supporting fragments, matching
            # all available evidence still earns the full pool credit.
            support_target = min(total_count, 3)
            pool_quality = (
                min(matched_count / support_target, 1.0)
                if support_target > 0
                else 0.0
            )
        else:
            if quality_kind == 2:
                pool_quality = _top_two_fragment_quality(
                    pool_matches, quality_relative_intensity_overrides, half_saturation,
                )
            elif quality_kind == 3:
                pool_quality = _positive_three_substituent_top_half_quality(
                    record,
                    pool_matches,
                    quality_relative_intensity_overrides,
                    half_saturation,
                )
            else:
                pool_quality = max(
                    (
                        _saturation_fragment_quality(
                            _match_relative_intensity(
                                match,
                                quality_relative_intensity_overrides=quality_relative_intensity_overrides,
                            ),
                            half_saturation,
                        )
                        for match in pool_matches
                    ),
                    default=0.0,
                )
        intensity_ratio = pool_quality
        pool_score = dynamic_pool_weight * pool_quality
        result[pool_name] = PoolScore(
            pool_name=pool_name,
            matched_count=matched_count,
            total_count=total_count,
            count_ratio=count_ratio,
            intensity_ratio=intensity_ratio,
            weight_ratio=0.0,
            pool_score=pool_score,
        )
    return _lcb_other_support(result, matches, record, rule, fragments,
                                  quality_relative_intensity_overrides, spectrum)


@cached_record_fact
def _uses_lcb_other_support(record):
    return (record.compound_class in {"Cer", "HexCer", "SHexCer", "SHexCer+O"} and record.adduct == "[M+H]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name))


def _lcb_other_support(result, matches, record, rule, fragments, overrides, spectrum):
    """Share auxiliary support across the enabled positive ceramide families."""
    if not _uses_lcb_other_support(record):
        return result
    library = record.fragments if fragments is None else fragments
    memberships = [_scoring_pools_for_fragment(record, f, matched=False) for f in library]
    pools = {id(f): membership for f, membership in zip(library, memberships)}
    lcb_total = sum("lcb" in value for value in memberships)
    if not lcb_total or "lcb" not in result:
        return result
    def route(match):
        membership = pools.get(id(match.fragment))
        return membership if membership is not None else _scoring_pools_for_fragment(record, match.fragment, matched=True)
    lcb_matches = [m for m in matches if "lcb" in route(m)]
    if overrides is None and spectrum is not None and len(lcb_matches) > 1:
        overrides = _non_precursor_quality_overrides(spectrum, record, matches)
    primary = max(lcb_matches, key=lambda m: (
        _match_relative_intensity(m, overrides), float(m.experimental_peak.relative_intensity)), default=None)
    # A centroid can provide either LCB intensity or Other support, once.
    seen = {float(primary.experimental_peak.mz)} if primary else set()
    support, auxiliary = [], []
    for m in sorted(matches, key=lambda m: "other" not in route(m)):
        membership = route(m)
        if not ("other" in membership or "lcb" in membership and m is not primary):
            continue
        peak = m.experimental_peak
        key = float(peak.mz)
        if peak.relative_intensity > OTHER_MIN_RELATIVE_INTENSITY and key not in seen:
            seen.add(key)
            support.append(m)
            if "lcb" in membership:
                auxiliary.append(m)
    positions = {id(f): i for i, f in enumerate(record.fragments)}
    def indexes(items):
        return tuple(positions[id(m.fragment)] for m in items if id(m.fragment) in positions)
    result["lcb"] = replace(result["lcb"], scoring_fragment_indexes=indexes([primary] if primary else []))
    # Theoretical support capacity is fixed independently of this spectrum.
    total = sum("other" in value for value in memberships) + max(lcb_total - 1, 0)
    target = min(total, 3)
    quality = min(len(support) / target, 1.0) if target else 0.0
    weight = _pool_weights_for_record(record, rule)["other"]
    result["other"] = replace(result["other"], matched_count=len(support), total_count=total,
                              count_ratio=len(support) / total if total else 0.0,
                              intensity_ratio=quality, pool_score=weight * quality,
                              scoring_fragment_indexes=indexes(support),
                              auxiliary_fragment_indexes=indexes(auxiliary))
    return result


def _match_relative_intensity(
    match: FragmentMatch,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> float:
    if quality_relative_intensity_overrides is not None and id(match) in quality_relative_intensity_overrides:
        relative_intensity = quality_relative_intensity_overrides[id(match)]
    else:
        relative_intensity = match.experimental_peak.relative_intensity
    return min(max(float(relative_intensity), 0.0), 1.0)


def _saturation_fragment_quality(relative_intensity: float, half_saturation: float) -> float:
    intensity = min(max(float(relative_intensity), 0.0), 1.0)
    if intensity <= 0.0:
        return 0.0
    return min(
        ((1.0 + float(half_saturation)) * intensity) / (intensity + float(half_saturation)),
        1.0,
    )


def _top_two_fragment_quality(matches, quality_relative_intensity_overrides, half_saturation):
    """Average up to two matched physical peaks; absent peaks are not zeros."""
    by_peak = {}
    for match in matches:
        mz = float(match.experimental_peak.mz)
        quality = _saturation_fragment_quality(
            _match_relative_intensity(match, quality_relative_intensity_overrides),
            half_saturation,
        )
        by_peak[mz] = max(by_peak.get(mz, 0.0), quality)
    selected = sorted(by_peak.values(), reverse=True)[:2]
    return sum(selected) / len(selected) if selected else 0.0


@cached_spectrum_fact
def _non_precursor_base_intensity(spectrum: ExperimentalSpectrum) -> float:
    bridge = native_bridge.load_bridge()
    query = get_native_query(spectrum, exclusion=PRECURSOR_CLUSTER_EXCLUSION_DA)
    if query is not None and bridge is not None and hasattr(bridge, 'query_base'):
        return bridge.query_base(query)[0]
    if bridge is not None and hasattr(bridge, 'non_precursor_base'):
        native = bridge.non_precursor_base(spectrum, PRECURSOR_CLUSTER_EXCLUSION_DA)
        if native is not None:
            return native[0]
    return max(
        (
            peak.intensity
            for peak in spectrum.peaks
            if abs(float(peak.mz) - float(spectrum.precursor_mz)) > PRECURSOR_CLUSTER_EXCLUSION_DA
        ),
        default=0.0,
    )


def _non_precursor_quality_overrides(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> Dict[int, float]:
    # Free fatty acids are intentionally precursor-driven.  Their precursor
    # evidence must retain the spectrum-wide normalization used by the
    # dedicated FA/OxFA scoring path.
    if _normalized_class_key(record.compound_class) in NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS:
        return {}

    bridge = native_bridge.load_bridge()
    query = get_native_query(spectrum, exclusion=PRECURSOR_CLUSTER_EXCLUSION_DA)
    if query is not None and bridge is not None and hasattr(bridge, 'quality_overrides_query'):
        native = bridge.quality_overrides_query(matches, query)
        if native is not None:
            return native

    non_precursor_base_intensity = _non_precursor_base_intensity(spectrum)
    if non_precursor_base_intensity <= 0.0:
        return {}

    bridge = native_bridge.load_bridge()
    if bridge is not None and hasattr(bridge, "quality_overrides"):
        native = bridge.quality_overrides(
            matches, spectrum.precursor_mz, non_precursor_base_intensity,
            PRECURSOR_CLUSTER_EXCLUSION_DA,
        )
        if native is not None:
            return native

    overrides: Dict[int, float] = {}
    for match in matches:
        if match.fragment.fragment_type == "Precursor Ion":
            continue
        if abs(float(match.experimental_peak.mz) - float(spectrum.precursor_mz)) <= PRECURSOR_CLUSTER_EXCLUSION_DA:
            continue
        adjusted_relative_intensity = min(
            float(match.experimental_peak.intensity) / non_precursor_base_intensity,
            1.0,
        )
        if adjusted_relative_intensity > float(match.experimental_peak.relative_intensity):
            overrides[id(match)] = adjusted_relative_intensity
    return overrides


def _native_pool_normalization(spectrum, record):
    """Supply only this query's scalars; unusual values keep reference math."""
    if (_non_precursor_quality_overrides is not _PRODUCTION_QUALITY_OVERRIDES
            or _normalized_class_key(record.compound_class) in NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS):
        return None
    query = get_native_query(spectrum, exclusion=PRECURSOR_CLUSTER_EXCLUSION_DA)
    if query is not None:
        return query
    base = _non_precursor_base_intensity(spectrum)
    if (type(base) in (float, int) and type(spectrum.precursor_mz) in (float, int)
            and base > 0.):
        return spectrum.precursor_mz, base, PRECURSOR_CLUSTER_EXCLUSION_DA
    return None


def _matched_intensity_sums(spectrum, matches):
    query = get_native_query(spectrum)
    bridge = native_bridge.load_bridge()
    if query is not None and bridge is not None and hasattr(bridge, 'query_intensity_sums'):
        sums = bridge.query_intensity_sums(matches, query)
        if sums is not None:
            return sums
    return (sum(match.experimental_peak.intensity for match in matches),
            sum(match.experimental_peak.relative_intensity for match in matches))


def _total_score_from_pool_scores(pool_scores: Dict[str, PoolScore]) -> float:
    return min(SCORE_MAX, max(SCORE_MIN, sum(pool.pool_score for pool in pool_scores.values())))


def _with_negative_cer_match_bonus(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    total_score: float,
) -> float:
    """Add ten points for >=10 physical product-ion matches in negative Cer."""
    if (
        _normalized_class_key(record.compound_class) == 'CER'
        and not _is_positive_adduct(record.adduct)
        and len(matches) >= 10
    ):
        return min(SCORE_MAX, float(total_score) + 10.0)
    return float(total_score)


def _required_chain_level_fah_hits(
    record: LibraryRecord,
    expected_fah_tokens: Sequence[str],
    chain_tokens: Sequence[str],
    rule: ClassRule,
) -> int:
    token_count = len(set(expected_fah_tokens)) if expected_fah_tokens else _count_fah_expected_chains(chain_tokens)
    if _is_cl_double_negative_record(record):
        return _cl_required_fah_hit_count(token_count)
    return min(max(token_count, 1), rule.chain_level_min_fah)


def _determine_resolution(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    is_positive_mode: bool,
    positive_loss_can_resolve_chain: bool,
) -> Tuple[str, str]:
    if is_negative_lpe(record) or is_negative_lpc(record) or is_positive_ether_lpc(record):
        return "species_level", ""
    chain_tokens = _extract_chain_tokens(record.lipid_chain_name)
    if oxidized_chain_requirements(record) is not None:
        if not oxidized_chain_missing_groups(record, matches):
            return "chain_level", ""
        return "species_level", "insufficient_oxidized_chain_evidence"
    if not chain_tokens:
        return "species_level", "missing_chain_annotation"
    configured_resolution = profile_chain_resolution(record, matches)
    if configured_resolution is not None and configured_resolution[0] == "chain_level":
        return configured_resolution
    if is_positive_mode:
        if (
            expected_fah_tokens
            and _matched_required_fah_count(record, matches, expected_fah_tokens)
                >= _required_chain_level_fah_hits(record, expected_fah_tokens, chain_tokens, rule)
        ):
            return "chain_level", ""
        if _positive_glyceride_rco_gate_passes(record, matches):
            return "chain_level", ""
        if positive_loss_can_resolve_chain:
            required_loss_hits = _positive_required_loss_hits(chain_tokens, _record_fa_loss_fragment_count(record))
            matched_chain_loss_count = _chain_evidence_count_for_matches(
                record, matches, _fragment_counts_as_chain_resolving_loss)
            if matched_chain_loss_count >= required_loss_hits:
                return "chain_level", ""
        if _matched_hg_fragment_count(record, matches) >= 1:
            if configured_resolution is not None:
                return configured_resolution
            has_library_loss = _record_loss_fragment_count(record) > 0
            if record.compound_class in POSITIVE_HG_CHAIN_LEVEL_CLASSES:
                return "chain_level", ""
            if not expected_fah_tokens and not has_library_loss:
                if not rule.allow_hg_only_if_no_fah:
                    return "chain_level", ""
                if rule.positive_hg_complete_can_resolve_chain:
                    return "chain_level", ""
                required_signature_hits = min(
                    max(int(rule.positive_signature_min_matches_if_no_fah), 0),
                    _record_positive_signature_fragment_count(record),
                )
                if required_signature_hits > 0:
                    matched_signature_count = _matched_positive_signature_fragment_count(matches)
                    if matched_signature_count >= required_signature_hits:
                        return "chain_level", ""
                return "class_level", "lyso_hg_only_fallback"
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        if _chain_evidence_count_for_matches(record, matches, _fragment_counts_as_effective_loss) >= 1:
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        return "class_level", "missing_positive_gate"
    if not expected_fah_tokens and rule.allow_loss_only_if_no_fah:
        matched_loss_count = sum(1 for match in matches if match.fragment.fragment_type == "Diagnostic_FA_Loss")
        if matched_loss_count >= 1:
            return "chain_level", ""
        return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
    required_fah_hits = _required_chain_level_fah_hits(record, expected_fah_tokens, chain_tokens, rule)
    if _matched_required_fah_count(record, matches, expected_fah_tokens) >= required_fah_hits:
        return "chain_level", ""
    return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON


@cached_rule_fact
def _prepared_fa_precursor_pool_plan(record, rule):
    """Compile the unchanged FA precursor-subset configuration, not a score."""
    kernel = native_bridge.load_bridge()
    if kernel is None or not hasattr(kernel, 'create_pool_plan'):
        return None
    fragments = _record_precursor_ion_fragments(record)
    totals = defaultdict(list)
    for fragment in fragments:
        for name in _scoring_pools_for_fragment(record, fragment, matched=False):
            totals[name].append(fragment)
    configuration = _build_pool_configuration(record, rule, totals)
    bits = {item[0]: 1 << i for i, item in enumerate(configuration)}
    routes = {id(fragment): (sum(bits.get(name, 0) for name in
               _scoring_pools_for_fragment(record, fragment, matched=True)), -1)
              for fragment in record.fragments}
    capsule = kernel.create_pool_plan(routes, configuration, (0, False))
    return None if capsule is None else (capsule, configuration,
        tuple(id(f) for f in fragments), bytes(1024 + len(routes)*64))


def _fa_precursor_pool_scores(spectrum, record, rule, matches, fragments):
    if (not _is_fa_record(record) or _calculate_pool_scores is not _FA_REFERENCE_POOL_SCORER
            or _build_pool_configuration is not _FA_REFERENCE_CONFIGURATION
            or _scoring_pools_for_fragment is not _FA_REFERENCE_ROUTING
            or PoolScore is not _FA_REFERENCE_POOL_TYPE
            or PoolScore.__init__ is not _FA_REFERENCE_POOL_INITIALIZER
            or PoolScore.__init__.__code__ is not _FA_REFERENCE_POOL_CODE
            or _saturation_fragment_quality is not _PRODUCTION_SATURATION
            or _match_relative_intensity is not _PRODUCTION_MATCH_INTENSITY
            or _top_two_fragment_quality is not _PRODUCTION_TOP_TWO
            or type(record) is not LibraryRecord
            or type(record.compound_class) is not str or type(record.adduct) is not str
            or any(type(f) is not FragmentRecord or type(f.name) is not str
                   or type(f.fragment_type) is not str for f in record.fragments)):
        return None
    query = get_native_query(spectrum)
    kernel = native_bridge.load_bridge()
    if query is None or kernel is None or not hasattr(kernel, 'pool_metrics_query'):
        return None
    prepared = _prepared_fa_precursor_pool_plan(record, rule)
    if prepared is None or tuple(id(f) for f in fragments) != prepared[2]:
        return None
    numeric = kernel.pool_metrics_query(matches, prepared[0], None, query)
    if numeric is None:
        return None
    return {item[0]: PoolScore(item[0], values[0], item[1], values[1], values[2], 0.0, values[3])
            for item, values in zip(prepared[1], numeric)}


def _score_fa_precursor_only_candidate(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    rule: ClassRule,
    precursor_fragments: Sequence[FragmentRecord],
    precursor_ppm_tolerance: float,
    precursor_mz_tolerance_da: float | None,
    fragment_mz_tolerance: float | None,
    fragment_ppm_tolerance: float | None,
    experimental_mz: Sequence[float] | None,
    precomputed_matches: Sequence[FragmentMatch] | None = None,
) -> CandidateScore:
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    mz_error_da = abs(spectrum.precursor_mz - record.precursor_mz)
    precursor_out_of_tolerance = (
        mz_error_da > precursor_mz_tolerance_da
        if precursor_mz_tolerance_da is not None
        else abs(ppm_error) > precursor_ppm_tolerance
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["precursor"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )

    matches = precomputed_matches if precomputed_matches is not None else _match_fragments(
        spectrum.peaks,
        precursor_fragments,
        fragment_mz_tolerance,
        ppm_tolerance=fragment_ppm_tolerance,
        experimental_mz=experimental_mz,
    )
    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["precursor"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )

    matched_intensity_sum, matched_relative_intensity_sum = _matched_intensity_sums(spectrum, matches)
    pool_scores = _fa_precursor_pool_scores(spectrum, record, rule, matches, precursor_fragments)
    if pool_scores is None:
        pool_scores = _calculate_pool_scores(matches, record, rule, fragments=precursor_fragments)
    total_score = _total_score_from_pool_scores(pool_scores)
    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=True,
        missing_required_groups=[],
        ppm_error=ppm_error,
        resolution_level="chain_level",
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason="",
    )


_prepared_scoring_module = None


@with_record_facts
def score_candidate(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    rule: ClassRule,
    precursor_ppm_tolerance: float = 10.0,
    precursor_mz_tolerance_da: float | None = None,
    fragment_mz_tolerance: float | None = 0.01,
    fragment_ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
    precomputed_matches: Sequence[FragmentMatch] | None = None,
    _search_only: bool = False,
    _native_evidence_checked=None,
) -> CandidateScore | None:
    rule = single_chain_rule(record, rule)
    if _normalized_class_key(record.compound_class) == "NAE":
        from .nae_precursor import nae_precursor_fragment
        # Re-match the single precursor, so legacy/precomputed HG or water
        # matches cannot satisfy this gate. Precursor metadata alone is not a hit.
        return _score_fa_precursor_only_candidate(
            spectrum, record, rule, [nae_precursor_fragment(record)],
            precursor_ppm_tolerance, precursor_mz_tolerance_da,
            fragment_mz_tolerance, fragment_ppm_tolerance, experimental_mz,
        )
    fa_precursor_fragments = _record_precursor_ion_fragments(record) if _is_fa_record(record) else []
    if fa_precursor_fragments:
        return _score_fa_precursor_only_candidate(
            spectrum=spectrum,
            record=record,
            rule=rule,
            precursor_fragments=fa_precursor_fragments,
            precursor_ppm_tolerance=precursor_ppm_tolerance,
            precursor_mz_tolerance_da=precursor_mz_tolerance_da,
            fragment_mz_tolerance=fragment_mz_tolerance,
            fragment_ppm_tolerance=fragment_ppm_tolerance,
            experimental_mz=experimental_mz,
            precomputed_matches=precomputed_matches,
        )

    if (_search_only and _native_evidence_checked is _NATIVE_EVIDENCE_PRECHECK
            and _missing_required_groups is _PRODUCTION_MISSING_GROUPS
            and _determine_resolution is _PRODUCTION_RESOLUTION):
        global _prepared_scoring_module
        if _prepared_scoring_module is None:
            from . import prepared_scoring as _prepared_scoring_module
        prepared = _prepared_scoring_module.score_prepared(
            spectrum, record, rule, precomputed_matches,
            precursor_mz_tolerance_da, precursor_ppm_tolerance)
        if prepared is not NotImplemented:
            return prepared

    # The public scorer retains complete rejected-candidate evidence. The
    # canonical searcher may omit work for candidates that cannot reach either
    # a strict result or an existing HG/FAH fallback. Matching is unchanged.
    search_matches = None
    search_support_only = False
    if _search_only and _native_evidence_checked is _NATIVE_EVIDENCE_PRECHECK:
        # Only the verified indexed plan may supply this token. Its program
        # already checked the search-only necessary conditions. The exact
        # mandatory gates, fallback intensities and resolution still run below.
        search_matches = precomputed_matches
        if not search_matches:
            return None
    elif _search_only:
        search_matches = precomputed_matches if precomputed_matches is not None else _match_fragments(
            spectrum.peaks, record.fragments, fragment_mz_tolerance,
            ppm_tolerance=fragment_ppm_tolerance, experimental_mz=experimental_mz,
        )
        if not search_matches:
            return None
        if not _is_positive_adduct(record.adduct) and oxidized_chain_requirements(record) is None:
            expected_fah = set(_record_expected_fah_tokens(record))
            if expected_fah and not expected_fah.issubset(_matched_fah_tokens(search_matches)):
                return None
        if (
            (is_positive_pep(record) or is_positive_pc_sodium(record) or is_positive_pe_sodium(record))
            and not any(_fragment_counts_as_hg(record, match.fragment) for match in search_matches)
        ):
            return None
        if sodium_phospholipid_missing_hg(record, search_matches):
            return None
        search_support_only = (
            _normalized_class_key(record.compound_class) != "OXFA"
            and all(
                match.fragment.fragment_type in {"Common", "Precursor Ion"}
                or (
                    match.fragment.fragment_type == "FA_Frag"
                    and not _is_positive_glyceride_chain_gate_fragment(record, match.fragment)
                    and not _fragment_counts_as_effective_loss(record, match.fragment)
                )
                for match in search_matches
            )
            and not _requires_precursor_fragment(record)
            and not any(_fragment_counts_as_hg(record, match.fragment) for match in search_matches)
        )
        if search_support_only and (
            _record_hg_fragment_count(record) > 0
            or _record_expected_fah_tokens(record)
            or (_is_positive_adduct(record.adduct) and _record_loss_fragment_count(record) > 0)
        ):
            return None

    (
        expected_fah_tokens,
        hg_fragment_count,
        loss_fragment_count,
        fa_loss_fragment_count,
        is_positive_mode,
        positive_loss_can_resolve_chain,
        allow_lyso_hg_only,
        allow_loss_only,
        require_loss_with_fah_only,
        supports_required_groups,
        required_group_names,
    ) = _derive_required_groups(record, rule)
    if (
        _search_only and required_group_names and search_support_only
    ):
        return None
    positive_required_hg_hits = _required_positive_hg_hits(
        rule=rule,
        expected_fah_tokens=expected_fah_tokens,
        hg_fragment_count=hg_fragment_count,
    )
    if not supports_required_groups and required_group_names:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=0.0,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="library_missing_required_fragments",
        )
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    mz_error_da = abs(spectrum.precursor_mz - record.precursor_mz)
    precursor_out_of_tolerance = (
        mz_error_da > precursor_mz_tolerance_da
        if precursor_mz_tolerance_da is not None
        else abs(ppm_error) > precursor_ppm_tolerance
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )
    matches = search_matches if search_matches is not None else precomputed_matches if precomputed_matches is not None else _match_fragments(
        spectrum.peaks,
        record.fragments,
        fragment_mz_tolerance,
        ppm_tolerance=fragment_ppm_tolerance,
        experimental_mz=experimental_mz,
    )
    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )
    chain_tokens = _extract_chain_tokens(record.lipid_chain_name)
    missing_groups = _missing_required_groups(
        record,
        rule,
        expected_fah_tokens,
        hg_fragment_count,
        positive_required_hg_hits,
        loss_fragment_count,
        chain_tokens,
        fa_loss_fragment_count,
        is_positive_mode,
        allow_loss_only,
        require_loss_with_fah_only,
        matches,
    )
    if _search_only and missing_groups:
        missing_set = set(missing_groups)
        # Negative fallback admits only missing HG. Positive fallbacks can
        # recover HG/FAH/loss, never failed precursor, TG all-chain, oxidized-
        # chain or other mandatory gates. Leave every recoverable case intact.
        recoverable = (
            missing_set.issubset({"hg", "sphingo_rule"})
            or missing_set.issubset({"fah", "loss"})
        ) if is_positive_mode else missing_set == {"hg"}
        if not recoverable:
            return None
    pool_scores = _calculate_pool_scores(
        matches,
        record,
        rule,
        spectrum=spectrum,
    )
    matched_intensity_sum, matched_relative_intensity_sum = _matched_intensity_sums(spectrum, matches)
    matched_hg_count = _matched_hg_fragment_count(record, matches)
    total_score = _total_score_from_pool_scores(pool_scores)
    positive_single_chain_species_hg_only = _positive_single_chain_species_hg_only(
        chain_tokens=chain_tokens,
        is_positive_mode=is_positive_mode,
        allow_lyso_hg_only=allow_lyso_hg_only,
        expected_fah_tokens=expected_fah_tokens,
    )
    resolution_level, downgrade_reason = _determine_resolution(
        record,
        matches,
        rule,
        expected_fah_tokens,
        is_positive_mode,
        positive_loss_can_resolve_chain,
    )
    fah_only_low_confidence_gate = _fah_only_low_confidence_gate_passes(
        record,
        pool_scores,
        missing_groups,
    )
    hg_only_low_confidence_gate = (
        resolution_level != "chain_level"
        and _hg_only_low_confidence_gate_passes(
            record,
            pool_scores,
            matches,
            missing_groups,
        )
    )
    weak_hg_only_candidate = (
        resolution_level != "chain_level"
        and not missing_groups
        and _is_positive_hg_dominant_phospholipid(record)
        and matched_hg_count > 0
        and not hg_only_low_confidence_gate
    )
    if fah_only_low_confidence_gate:
        missing_groups = []
        resolution_level = "tentative_chain_level"
        downgrade_reason = FAH_ONLY_FALLBACK_REASON
    elif hg_only_low_confidence_gate:
        missing_groups = []
        resolution_level = "tentative_species_level"
        downgrade_reason = HG_ONLY_FALLBACK_REASON
    elif weak_hg_only_candidate:
        missing_groups = ["hg"]
        resolution_level = "class_level"
        downgrade_reason = "missing_required_hg"
    elif missing_groups:
        resolution_level = "class_level"
        downgrade_reason = "missing_required_" + ",".join(missing_groups)
    elif allow_lyso_hg_only and not expected_fah_tokens:
        if resolution_level == "chain_level":
            downgrade_reason = ""
        elif positive_single_chain_species_hg_only:
            resolution_level = "species_level"
            downgrade_reason = "lyso_hg_only_fallback"
        elif (
            rule.positive_hg_complete_can_resolve_chain
            and matched_hg_count >= positive_required_hg_hits
        ):
            resolution_level = "chain_level"
            downgrade_reason = ""
        else:
            resolution_level = "class_level"
            downgrade_reason = "lyso_hg_only_fallback"
    chain_decision = positive_phospholipid_chain_resolution(record, matches)
    if not missing_groups and chain_decision is not None and chain_decision[1]:
        resolution_level, downgrade_reason = chain_decision
        if downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON:
            missing_groups = ["chain_confirmation"]
    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=not missing_groups,
        missing_required_groups=missing_groups,
        ppm_error=ppm_error,
        resolution_level=resolution_level,
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason=downgrade_reason,
    )


_PRODUCTION_SATURATION = _saturation_fragment_quality
_PRODUCTION_MATCH_INTENSITY = _match_relative_intensity
_PRODUCTION_TRIPLE_QUALITY = _positive_three_substituent_top_half_quality
_PRODUCTION_TOP_TWO = _top_two_fragment_quality
_PRODUCTION_MISSING_GROUPS = _missing_required_groups
_PRODUCTION_RESOLUTION = _determine_resolution
_PRODUCTION_QUALITY_OVERRIDES = _non_precursor_quality_overrides
_FA_REFERENCE_POOL_SCORER = _calculate_pool_scores
_FA_REFERENCE_CONFIGURATION = _build_pool_configuration
_FA_REFERENCE_ROUTING = _scoring_pools_for_fragment
_FA_REFERENCE_POOL_TYPE = PoolScore
_FA_REFERENCE_POOL_INITIALIZER = PoolScore.__init__
_FA_REFERENCE_POOL_CODE = PoolScore.__init__.__code__
