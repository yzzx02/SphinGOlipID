"""Class evidence gates shared by generic and sphingolipid scoring."""
from __future__ import annotations

import re
import math
from collections import Counter
from dataclasses import dataclass, replace
from typing import Sequence, Dict, List
from .chain_utils import extract_chain_tokens, extract_fragment_chain_token
from .evidence_policy import (
    is_positive_class,
    is_positive_pep,
    is_positive_ps,
    is_positive_ps_chain_ketene_loss,
    is_positive_ps_headgroup_loss,
    _chain_evidence_count_for_matches,
    _fragment_counts_as_effective_loss,
    _fragment_counts_as_fa_loss_gate,
    _fragment_counts_as_hg,
    _is_cl_double_negative_record,
    _is_positive_adduct,
    _is_positive_glyceride_rco_gate_record,
    _is_positive_hg_dominant_phospholipid,
    _is_positive_oxtg_full_fah_record,
    _matched_cl_chain_info_fragment_count,
    _matched_hg_fragment_count,
    _matched_hg_max_relative_intensity,
    _matched_positive_glyceride_rco_fragment_count,
    _matched_precursor_ion_fragment_count,
    _matched_required_fah_count,
    _normalized_class_key,
    _normalized_compound_class,
    _positive_glyceride_chain_gate_group,
    _positive_required_loss_hits,
    _record_expected_fah_tokens,
    _record_fa_loss_fragment_count,
    _record_hg_fragment_count,
    _record_loss_fragment_count,
    _record_positive_glyceride_rco_fragment_count,
    _requires_precursor_fragment,
    _uses_positive_support_pool_scoring,
    sphingolipid_fragment_types,
)
from .models import (
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    PoolScore,
    CandidateScore,
)
from .oxidized_chain_evidence import oxidized_chain_requirements, oxidized_chain_missing_groups
from .positive_hexosylceramide import is_positive_ganglioside
from .positive_ganglioside import is_positive_ganglioside as is_positive_extended_ganglioside, from_name as extended_ganglioside_from_name
from .negative_lpe import is_negative_lpe
from .negative_lpc import is_negative_lpc
from .positive_lpc_ether import is_positive_ether_lpc
from .positive_pc_sodium import is_positive_pc_sodium, HG_NAMES as PC_SODIUM_HG_NAMES
from .positive_pe_cer import is_positive_pe_cer, from_name as pe_cer_from_name
from .positive_pe_sodium import is_positive_pe_sodium, HG_NAMES as PE_SODIUM_HG_NAMES
from .resolution_policy import (
    fragment_is_chain_evidence,
    normalized_class_key,
    record_polarity,
    PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON,
)
from .rules import ClassRule
from .record_facts import cached_rule_fact, cached_record_fact, cached_fragment_fact, _current_facts
from .scoring_policy import (
    FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE,
    FAH_ONLY_FALLBACK_MIN_MATCHES,
    FAH_ONLY_FALLBACK_MIN_POOL_SCORE,
    HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY,
)
from .sphingolipid_rules import SphingoRule, validate_rule


POSITIVE_TG_FULL_CHAIN_CLASSES = frozenset({"TG"})
POSITIVE_TG_O_FULL_LOSS_CLASSES = frozenset({"TGO"})
POSITIVE_TG_EST_FULL_LOSS_CLASSES = frozenset({"TGEST"})
POSITIVE_THREE_SUBSTITUENT_GLYCERIDE_CLASSES = (
    POSITIVE_TG_FULL_CHAIN_CLASSES
    | POSITIVE_TG_O_FULL_LOSS_CLASSES
    | POSITIVE_TG_EST_FULL_LOSS_CLASSES
)


def candidate_passes_evidence_gate(result: CandidateScore) -> bool:
    """Admit the documented tentative chain result with its failed gate intact."""
    return result.passed_required_gates or (
        result.downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
        and result.missing_required_groups == ["chain_confirmation"]
    )


def sodium_phospholipid_missing_hg(
    record: LibraryRecord, matches: Sequence[FragmentMatch],
) -> bool:
    """Share the existing mandatory PE/PC sodium HG checks with search routing."""
    if is_positive_pe_sodium(record):
        names = {m.fragment.name for m in matches if m.fragment.fragment_type == "Diagnostic_HG"}
        return not PE_SODIUM_HG_NAMES.issubset(names)
    if is_positive_pc_sodium(record):
        peaks = {float(m.experimental_peak.mz) for m in matches
                 if m.fragment.fragment_type == "Diagnostic_HG" and m.fragment.name in PC_SODIUM_HG_NAMES}
        return len(peaks) < 2
    return False


def positive_pep_chain_gate_counts(
    record: LibraryRecord, matches: Sequence[FragmentMatch],
) -> tuple[int, int]:
    """Require half the distinct chain ions, never duplicate aliases or peaks."""
    expected_mz = {float(f.mz) for f in record.fragments if fragment_is_chain_evidence(record, f)}
    observed_mz = {float(m.experimental_peak.mz) for m in matches
                   if fragment_is_chain_evidence(record, m.fragment)}
    return (len(expected_mz) + 1) // 2, len(observed_mz)


def is_positive_tg_full_chain_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_FULL_CHAIN_CLASSES)


def is_positive_tg_o_full_loss_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_O_FULL_LOSS_CLASSES)


def is_positive_three_substituent_glyceride(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_THREE_SUBSTITUENT_GLYCERIDE_CLASSES)


@cached_fragment_fact
def glyceride_substituent_evidence_group(
    record: LibraryRecord,
    fragment: FragmentRecord,
) -> str | None:
    """Map one diagnostic loss to one physical glycerol substituent group."""

    if not is_positive_three_substituent_glyceride(record):
        return None
    if fragment.fragment_type != "Diagnostic_FA_Loss":
        return None
    class_key = normalized_class_key(record.compound_class)
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    token = extract_fragment_chain_token(fragment)

    if class_key == "TG":
        if "M-NH3-(ROOH)+NH4" in compact_name and token is not None:
            return f"acyl:{token}"
        return None

    if class_key == "TGO":
        if compact_name.startswith("[M-R1-OH+H]+(") and token is not None:
            return f"ether:{token}"
        if "M-NH3-(ROOH)+NH4" in compact_name and token is not None:
            return f"acyl:{token}"
        return None

    if class_key == "TGEST":
        if compact_name.startswith("[M+H-FAHFA]+("):
            return "fahfa"
        if compact_name.startswith("[M+H-FA]+(") and token is not None:
            return f"fa:{token}"
    return None


@cached_record_fact
def glyceride_substituent_group_multiplicity(record: LibraryRecord) -> Counter[str]:
    """Derive three structural positions while retaining repeated groups."""

    class_key = normalized_class_key(record.compound_class)
    if class_key == "TG":
        return Counter(
            f"acyl:{token}" for token in extract_chain_tokens(record.lipid_chain_name)
        )
    if class_key == "TGO":
        groups: list[str] = []
        for token in extract_chain_tokens(record.lipid_chain_name):
            groups.append(
                f"ether:{token}"
                if token.upper().startswith(("O-", "P-"))
                else f"acyl:{token}"
            )
        return Counter(groups)
    if class_key == "TGEST":
        matched = re.fullmatch(
            r"TG-EST(?:\s+|\()(?P<fa1>\d+:\d+)_(?P<fa2>\d+:\d+)_"
            r"(?P<hfa>\d+:\d+);O\(FA\s+(?P<attached>\d+:\d+)\)\)?",
            str(record.lipid_chain_name or "").strip(),
            flags=re.IGNORECASE,
        )
        if matched is None:
            return Counter()
        return Counter(
            [
                f"fa:{matched.group('fa1')}",
                f"fa:{matched.group('fa2')}",
                "fahfa",
            ]
        )
    return Counter()


def expected_glyceride_substituent_groups(record: LibraryRecord) -> frozenset[str]:
    return frozenset(glyceride_substituent_group_multiplicity(record))


def positive_three_substituent_glyceride_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require all physically distinct substituent-loss groups in the record."""

    if not is_positive_three_substituent_glyceride(record):
        return False
    expected_groups = expected_glyceride_substituent_groups(record)
    if not expected_groups:
        return False
    matched_groups = {
        group
        for match in matches
        for group in [glyceride_substituent_evidence_group(record, match.fragment)]
        if group is not None
    }
    return expected_groups.issubset(matched_groups)


def positive_tg_full_chain_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require evidence for all three TG chain positions.

    A repeated chain is represented by one physical fragment but contributes
    its structural multiplicity.  The other distinct chain must still be
    observed, so TG(18:2_18:2_8:0) cannot pass on 18:2 alone.
    """

    return (
        is_positive_tg_full_chain_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )


def is_positive_tg_est_full_loss_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_EST_FULL_LOSS_CLASSES)


def positive_tg_o_full_loss_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require the ether-chain loss and all distinct acyl-chain losses."""

    return (
        is_positive_tg_o_full_loss_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )


def positive_tg_est_full_loss_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require FA1, FA2 and FAHFA losses, without forcing the attached FA3."""

    return (
        is_positive_tg_est_full_loss_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )


POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY = 0.50
ETHER_GLYCERIDE_HALF_RCO_GATE_CLASSES = {"TGO"}
POSITIVE_GLYCERIDE_RCO_MIN_HITS = 2
FAH_ONLY_FALLBACK_BLOCKED_CLASSES = {"ADGGA", "ASM", "AMPS"}
STRICT_NEGATIVE_HG_CLASSES = {"NAGPS"}


def positive_oxtg_full_fah_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    if not _is_positive_oxtg_full_fah_record(record):
        return False
    required_fragments = [
        fragment
        for fragment in record.fragments
        if fragment.fragment_type == "Diagnostic_FA_Loss"
    ]
    matched_fragment_ids = {
        id(match.fragment)
        for match in matches
        if match.fragment.fragment_type == "Diagnostic_FA_Loss"
    }
    return bool(required_fragments) and all(
        id(fragment) in matched_fragment_ids for fragment in required_fragments
    )


def positive_glyceride_rco_required_hits(record: LibraryRecord) -> int:
    available_hits = _record_positive_glyceride_rco_fragment_count(record)
    if available_hits <= 0:
        return 0
    if _normalized_class_key(record.compound_class) in ETHER_GLYCERIDE_HALF_RCO_GATE_CLASSES:
        return max(1, math.ceil(available_hits / 2))
    return POSITIVE_GLYCERIDE_RCO_MIN_HITS


def positive_glyceride_rco_gate_available(record: LibraryRecord) -> bool:
    available_hits = _record_positive_glyceride_rco_fragment_count(record)
    required_hits = positive_glyceride_rco_required_hits(record)
    return required_hits > 0 and available_hits >= required_hits


def positive_glyceride_rco_gate_passes(record: LibraryRecord, matches: Sequence[FragmentMatch]) -> bool:
    required_hits = positive_glyceride_rco_required_hits(record)
    return (
        positive_glyceride_rco_gate_available(record)
        and _matched_positive_glyceride_rco_fragment_count(record, matches) >= required_hits
    )


def required_positive_hg_hits(
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    hg_fragment_count: int,
) -> int:
    if hg_fragment_count <= 0:
        return 0

    if not rule.allow_hg_only_if_no_fah or expected_fah_tokens:
        required_hg_hits = 1
    else:
        required_hg_hits = min(
            hg_fragment_count,
            max(int(rule.positive_hg_min_matches_if_no_fah), 1),
        )

    ratio_required_hg_hits = required_fractional_hits(
        total_count=hg_fragment_count,
        min_fraction=float(getattr(rule, "positive_hg_min_fraction", 0.0)),
        min_matches=int(getattr(rule, "positive_hg_min_matches", 0)),
    )
    if ratio_required_hg_hits > 0:
        required_hg_hits = max(required_hg_hits, ratio_required_hg_hits)
    return min(required_hg_hits, hg_fragment_count)


def required_negative_hg_hits(
    record: LibraryRecord,
    hg_fragment_count: int,
    rule: ClassRule,
) -> int:
    if hg_fragment_count <= 0:
        return 0
    if is_negative_lpe(record) or is_negative_lpc(record):
        return 1
    configured_minimum = max(int(getattr(rule, "negative_hg_min_matches", 0)), 0)
    if configured_minimum > 0:
        return min(hg_fragment_count, configured_minimum)
    if _normalized_compound_class(record.compound_class) in STRICT_NEGATIVE_HG_CLASSES:
        return hg_fragment_count
    return min(hg_fragment_count, max(1, math.ceil(hg_fragment_count / 2)))


def required_fractional_hits(total_count: int, min_fraction: float, min_matches: int) -> int:
    if total_count <= 0:
        return 0
    min_fraction = max(float(min_fraction), 0.0)
    min_matches = max(int(min_matches), 0)
    if min_fraction <= 0.0 and min_matches <= 0:
        return 0
    required_by_ratio = math.ceil(total_count * min_fraction) if min_fraction > 0.0 else 0
    required_hits = max(required_by_ratio, min_matches)
    required_hits = max(required_hits, 1)
    return min(required_hits, total_count)


@cached_rule_fact
def derive_required_groups(
    record: LibraryRecord,
    rule: ClassRule,
) -> tuple[List[str], int, int, int, bool, bool, bool, bool, bool, bool, List[str]]:
    fah_tokens = _record_expected_fah_tokens(record)
    hg_fragment_count = _record_hg_fragment_count(record)
    loss_fragment_count = _record_loss_fragment_count(record)
    fa_loss_fragment_count = _record_fa_loss_fragment_count(record)
    is_positive_mode = _is_positive_adduct(record.adduct)
    positive_rco_gate_available = positive_glyceride_rco_gate_available(record)
    positive_hg_or_loss_gate = is_positive_mode and (hg_fragment_count > 0 or loss_fragment_count > 0)
    allow_lyso_hg_only = rule.allow_hg_only_if_no_fah and not fah_tokens and hg_fragment_count > 0
    allow_loss_only = (
        rule.allow_loss_only_if_no_fah
        and not fah_tokens
        and hg_fragment_count == 0
        and loss_fragment_count > 0
    )
    require_loss_with_fah_only = (
        rule.require_loss_with_fah_only
        and bool(fah_tokens)
        and hg_fragment_count == 0
        and loss_fragment_count > 0
    )
    positive_loss_can_resolve_chain = is_positive_mode and fa_loss_fragment_count > 0
    if positive_hg_or_loss_gate:
        supports_required_groups = True
    else:
        supports_required_groups = bool(fah_tokens) or allow_lyso_hg_only or allow_loss_only or positive_rco_gate_available
    required_groups: List[str] = []
    if is_positive_ether_lpc(record):
        return ([], hg_fragment_count, 0, 0, True, False, False, False, False,
                True, ["hg"])
    if is_negative_lpc(record):
        return ([], hg_fragment_count, 0, 0, False, False, False, False, False,
                True, ["hg"])
    if is_negative_lpe(record):
        return ([], hg_fragment_count, 0, 0, False, False, False, False, False,
                True, ["precursor", "hg"])
    if positive_hg_or_loss_gate:
        if hg_fragment_count > 0:
            required_groups.append("hg")
        elif loss_fragment_count > 0:
            required_groups.append("loss")
    elif positive_rco_gate_available:
        required_groups.append(_positive_glyceride_chain_gate_group(record))
    elif is_positive_mode and not fah_tokens and _normalized_class_key(record.compound_class) == "DG":
        # Older imports can contain only RCO and precursor support.
        required_groups.append("loss")
    else:
        if fah_tokens:
            required_groups.append("fah")
        if hg_fragment_count > 0 and (fah_tokens or allow_lyso_hg_only):
            required_groups.append("hg")
        if require_loss_with_fah_only or allow_loss_only:
            required_groups.append("loss")
    return (
        fah_tokens,
        hg_fragment_count,
        loss_fragment_count,
        fa_loss_fragment_count,
        is_positive_mode,
        positive_loss_can_resolve_chain,
        allow_lyso_hg_only,
        allow_loss_only,
        require_loss_with_fah_only,
        supports_required_groups,
        required_groups,
    )


def fah_only_low_confidence_gate_passes(
    record: LibraryRecord,
    pool_scores: Dict[str, PoolScore],
    missing_groups: Sequence[str] | None = None,
    has_hg_or_structural_pool: bool | None = None,
) -> bool:
    if is_positive_pc_sodium(record) or is_positive_pe_sodium(record) or is_positive_pep(record):
        return False
    if not _is_positive_adduct(record.adduct):
        return False
    if _normalized_class_key(record.compound_class) in FAH_ONLY_FALLBACK_BLOCKED_CLASSES:
        return False
    if _is_positive_glyceride_rco_gate_record(record):
        return False
    if missing_groups is not None:
        missing_group_set = set(missing_groups)
        if not missing_group_set or "fah" in missing_group_set:
            return False
        if not missing_group_set.issubset({"hg", "sphingo_rule"}):
            return False

    fah_score = pool_scores.get("fah")
    if fah_score is None or fah_score.total_count <= 0:
        return False
    hg_score = pool_scores.get("hg")
    if has_hg_or_structural_pool is None:
        has_hg_or_structural_pool = hg_score is not None and hg_score.total_count > 0
    if not has_hg_or_structural_pool:
        return False
    if hg_score is not None and hg_score.matched_count > 0:
        return False
    if fah_score.matched_count < FAH_ONLY_FALLBACK_MIN_MATCHES:
        return False
    if (
        _uses_positive_support_pool_scoring(record)
        and fah_score.intensity_ratio < POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY
    ):
        return False
    min_pool_score = FAH_ONLY_FALLBACK_MIN_POOL_SCORE
    if _is_positive_hg_dominant_phospholipid(record):
        min_pool_score = FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE
    return fah_score.pool_score >= min_pool_score


def hg_only_low_confidence_gate_passes(
    record: LibraryRecord,
    pool_scores: Dict[str, PoolScore],
    matches: Sequence[FragmentMatch],
    missing_groups: Sequence[str] | None = None,
) -> bool:
    if is_positive_pc_sodium(record) or is_positive_pep(record):
        return False
    if not _is_positive_hg_dominant_phospholipid(record):
        return False
    missing_group_set = set(missing_groups or [])
    if "hg" in missing_group_set:
        return False
    if not missing_group_set.issubset({"fah", "loss"}):
        return False
    hg_score = pool_scores.get("hg")
    if hg_score is None or hg_score.total_count <= 0 or hg_score.matched_count <= 0:
        return False
    return _matched_hg_max_relative_intensity(record, matches) >= HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY


def missing_required_groups(
    record: LibraryRecord,
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    hg_fragment_count: int,
    positive_required_hg_hits: int,
    loss_fragment_count: int,
    chain_tokens: Sequence[str],
    fa_loss_fragment_count: int,
    is_positive_mode: bool,
    allow_loss_only: bool,
    require_loss_with_fah_only: bool,
    matches: Sequence[FragmentMatch],
) -> List[str]:
    missing = []
    if is_negative_lpc(record) or is_positive_ether_lpc(record):
        return [] if _matched_hg_fragment_count(record, matches) >= 1 else ["hg"]
    if is_negative_lpe(record):
        if not any(m.fragment.fragment_type == "Precursor Ion"
                   and abs(m.fragment.mz-record.precursor_mz) <= 1e-6 for m in matches):
            missing.append("precursor")
        if _matched_hg_fragment_count(record, matches) < 1:
            missing.append("hg")
        return missing
    if is_positive_pe_sodium(record):
        if sodium_phospholipid_missing_hg(record, matches):
            missing.append("hg")
        return missing
    if is_positive_pep(record):
        if _matched_hg_fragment_count(record, matches) < max(positive_required_hg_hits, 1):
            missing.append("hg")
        required, observed = positive_pep_chain_gate_counts(record, matches)
        if required <= 0 or observed < required:
            missing.append("fah")
        return missing
    if oxidized_chain_requirements(record) is not None:
        required_hg = positive_required_hg_hits if is_positive_mode else required_negative_hg_hits(record, hg_fragment_count, rule)
        if _matched_hg_fragment_count(record, matches) < required_hg:
            missing.append("hg")
        return missing + oxidized_chain_missing_groups(record, matches)
    if is_positive_pc_sodium(record):
        chain_peaks = {float(m.experimental_peak.mz) for m in matches
                       if fragment_is_chain_evidence(record, m.fragment)}
        if sodium_phospholipid_missing_hg(record, matches):
            missing.append("hg")
        if len(chain_peaks) < 2:
            missing.append("fah")
        return missing
    if is_positive_mode:
        matched_hg_count = _matched_hg_fragment_count(record, matches)

        if _requires_precursor_fragment(record) and _matched_precursor_ion_fragment_count(matches) < 1:
            missing.append("precursor")

        if is_positive_tg_est_full_loss_record(record):
            if not positive_tg_est_full_loss_gate_passes(record, matches):
                missing.append("tg_est_all_losses")
            return missing

        if is_positive_tg_o_full_loss_record(record):
            if not positive_tg_o_full_loss_gate_passes(record, matches):
                missing.append("tg_o_all_losses")
            return missing

        if is_positive_tg_full_chain_record(record):
            if not positive_tg_full_chain_gate_passes(record, matches):
                missing.append("tg_all_chains")
            return missing

        if _is_positive_oxtg_full_fah_record(record):
            if not positive_oxtg_full_fah_gate_passes(record, matches):
                missing.append("oxtg_all_fah")
            return missing

        if positive_glyceride_rco_gate_passes(record, matches):
            return missing

        if (
            positive_glyceride_rco_gate_available(record)
            and hg_fragment_count == 0
            and loss_fragment_count == 0
            and not expected_fah_tokens
        ):
            missing.append(_positive_glyceride_chain_gate_group(record))
            return missing

        if hg_fragment_count > 0 and matched_hg_count < positive_required_hg_hits:
            missing.append("hg")

        if hg_fragment_count == 0 and loss_fragment_count > 0:
            # Loss-only positive records require chain-count-consistent coverage.
            required_loss_hits = _positive_required_loss_hits(chain_tokens, fa_loss_fragment_count)
            matched_loss_count = _chain_evidence_count_for_matches(record, matches, _fragment_counts_as_effective_loss)
            if matched_loss_count < required_loss_hits:
                missing.append("loss")

        ratio_required_loss_hits = 0
        if hg_fragment_count > 0 and loss_fragment_count > 0:
            ratio_required_loss_hits = required_fractional_hits(
                total_count=fa_loss_fragment_count if fa_loss_fragment_count > 0 else loss_fragment_count,
                min_fraction=float(getattr(rule, "positive_loss_min_fraction", 0.0)),
                min_matches=int(getattr(rule, "positive_loss_min_matches", 0)),
            )
            if ratio_required_loss_hits > 0:
                predicate = _fragment_counts_as_fa_loss_gate if fa_loss_fragment_count > 0 else _fragment_counts_as_effective_loss
                matched_loss_for_ratio = _chain_evidence_count_for_matches(record, matches, predicate)
                if matched_loss_for_ratio < ratio_required_loss_hits:
                    missing.append("loss")

        if hg_fragment_count == 0 and loss_fragment_count == 0 and expected_fah_tokens:
            required_fah_hits = required_fah_gate_hits(record, expected_fah_tokens)
            if _matched_required_fah_count(record, matches, expected_fah_tokens) < required_fah_hits:
                missing.append("fah")
        return missing

    if expected_fah_tokens:
        required_fah_hits = (
            len(set(expected_fah_tokens))
            if bool(getattr(rule, "negative_require_all_fah", False))
            else required_fah_gate_hits(record, expected_fah_tokens)
        )
        if _matched_required_fah_count(record, matches, expected_fah_tokens) < required_fah_hits:
            missing.append("fah")
    if hg_fragment_count > 0:
        matched_hg_count = sum(1 for match in matches if _fragment_counts_as_hg(record, match.fragment))
        if matched_hg_count < required_negative_hg_hits(record, hg_fragment_count, rule):
            missing.append("hg")
    if _requires_precursor_fragment(record) and _matched_precursor_ion_fragment_count(matches) < 1:
        missing.append("precursor")
    if _is_cl_double_negative_record(record) and _matched_cl_chain_info_fragment_count(record, matches) < 1:
        missing.append("chain_info")
    if allow_loss_only or require_loss_with_fah_only:
        matched_loss_count = sum(1 for match in matches if _fragment_counts_as_effective_loss(record, match.fragment))
        if matched_loss_count < 1:
            missing.append("loss")
    if record.compound_class == "LNAPE" and not is_positive_mode:
        from .library_fragment_policy import is_lnape_n_acyl_headgroup

        if not expected_fah_tokens and "fah" not in missing:
            missing.append("fah")
        if not any(
            match.fragment.fragment_type == "Diagnostic_HG"
            and is_lnape_n_acyl_headgroup(match.fragment.name)
            for match in matches
        ) and "hg" not in missing:
            missing.append("hg")
    return missing


def cl_required_fah_hit_count(token_count: int) -> int:
    """Require every distinct CL acyl-chain FAH for [M-2H]2- assignments."""
    return max(token_count, 0)


def required_fah_gate_hits(record: LibraryRecord, expected_fah_tokens: Sequence[str]) -> int:
    token_count = len(set(expected_fah_tokens))
    if _is_cl_double_negative_record(record):
        return cl_required_fah_hit_count(token_count)
    return token_count


@dataclass(frozen=True)
class SphingolipidGateEvidence:
    passed: bool
    library_types: frozenset[str]
    matched_types: frozenset[str]
    adduct: str
    class_key: str
    uses_negative_hg_fah_tiered_gate: bool
    has_negative_chain_evidence: bool


@cached_record_fact
def _sphingolipid_library_evidence(record):
    types = {id(f): frozenset(sphingolipid_fragment_types(record, f)) for f in record.fragments}
    library_types = frozenset().union(*types.values())
    adduct = str(record.adduct or "").strip()
    class_key = normalized_class_key(record.compound_class)
    tiered = ((adduct == "[M-H]-" and class_key in {"CER1P", "CERP"})
              or (adduct in {"[M+CH3COO]-", "[M+HCOO]-"} and class_key == "SM"))
    identity_valid = True
    if is_positive_pe_cer(record.compound_class, record.adduct):
        identity_valid = pe_cer_from_name(record.lipid_chain_name) is not None
    if is_positive_extended_ganglioside(record.compound_class,record.adduct):
        identity_valid=identity_valid and extended_ganglioside_from_name(record.lipid_chain_name,record.compound_class) is not None
    elif is_positive_ganglioside(record.compound_class, record.adduct):
        identity_valid = identity_valid and re.fullmatch(
            rf"{record.compound_class}\(d\d+:\d+/\d+:\d+\)", record.lipid_chain_name) is not None
    hg_count = sum(f.fragment_type == "Diagnostic_HG" for f in record.fragments)
    return types, library_types, adduct, class_key, tiered, identity_valid, hg_count


def evaluate_sphingolipid_gate(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    rule: SphingoRule,
    series: str,
) -> SphingolipidGateEvidence:
    """Evaluate name/type gates without awarding scores or applying fallbacks."""
    if validate_rule is _REFERENCE_VALIDATE_RULE:
        compiled = _evaluate_compiled_sphingolipid_gate(record, matches, rule, series)
        if compiled is not None:
            return compiled
    matched_names = {match.fragment.name for match in matches}
    name_gate_passed = validate_rule(rule, matched_names, series)
    (types_by_id, library_types, adduct, class_key, uses_negative_hg_fah_tiered_gate,
     identity_valid, hg_total_count) = _sphingolipid_library_evidence(record)
    name_gate_passed = name_gate_passed and identity_valid

    def evidence_types(fragment: FragmentRecord) -> set[str]:
        found = types_by_id.get(id(fragment))
        return found if found is not None else sphingolipid_fragment_types(record, fragment)

    matched_types = set().union(*(evidence_types(m.fragment) for m in matches))
    def matched_type_count(type_group: set[str]) -> int:
        return sum(1 for match in matches if evidence_types(match.fragment) & type_group)

    type_gate_passed = True
    has_explicit_type_gate = bool(
        rule.required_type_any_groups
        or rule.required_type_count_groups
        or rule.required_type_count_any_groups
        or rule.required_type_fraction_groups
        or rule.required_name_pattern_fraction_groups
    )
    if rule.required_type_any_groups:
        for type_group in rule.required_type_any_groups:
            if not (type_group & matched_types):
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_count_groups:
        for type_group, minimum_count in rule.required_type_count_groups:
            if matched_type_count(type_group) < minimum_count:
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_count_any_groups:
        for alternatives in rule.required_type_count_any_groups:
            if not any(matched_type_count(type_group) >= minimum_count for type_group, minimum_count in alternatives):
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_fraction_groups:
        for type_group, minimum_fraction in rule.required_type_fraction_groups:
            library_count = sum(
                1 for fragment in record.fragments if evidence_types(fragment) & type_group
            )
            required_count = max(1, math.ceil(library_count * float(minimum_fraction)))
            if library_count <= 0 or matched_type_count(type_group) < required_count:
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_name_pattern_fraction_groups:
        for name_pattern, minimum_fraction in rule.required_name_pattern_fraction_groups:
            library_count = sum(
                1
                for fragment in record.fragments
                if re.search(name_pattern, str(fragment.name or ""), flags=re.IGNORECASE)
            )
            matched_count = sum(
                1
                for match in matches
                if re.search(name_pattern, str(match.fragment.name or ""), flags=re.IGNORECASE)
            )
            required_count = max(1, math.ceil(library_count * float(minimum_fraction)))
            if library_count <= 0 or matched_count < required_count:
                type_gate_passed = False
                break
    if uses_negative_hg_fah_tiered_gate:
        hg_matched_count = matched_type_count({"Diagnostic_HG"})
        if hg_total_count <= 0 or hg_matched_count * 2 < hg_total_count:
            type_gate_passed = False
    if not has_explicit_type_gate:
        if "Diagnostic_HG" in library_types and "Diagnostic_HG" not in matched_types:
            type_gate_passed = False
        if "C类碎片" in library_types and "C类碎片" not in matched_types:
            type_gate_passed = False
        if "LCB碎片" in library_types and "LCB碎片" not in matched_types:
            type_gate_passed = False

    return SphingolipidGateEvidence(
        passed=name_gate_passed and type_gate_passed,
        library_types=frozenset(library_types),
        matched_types=frozenset(matched_types),
        adduct=adduct,
        class_key=class_key,
        uses_negative_hg_fah_tiered_gate=uses_negative_hg_fah_tiered_gate,
        has_negative_chain_evidence=(uses_negative_hg_fah_tiered_gate
                                     and matched_type_count({"Diagnostic_FA", "Diagnostic_FA_Loss"}) >= 1),
    )


_REFERENCE_VALIDATE_RULE = validate_rule


@cached_record_fact
def _supports_compiled_sphingolipid_gate(record):
    return (len(record.fragments)<=64
            and len({id(f) for f in record.fragments})==len(record.fragments))


def _evaluate_compiled_sphingolipid_gate(record, matches, rule, series):
    current = _current_facts.get()
    if (current is None or current[0] is not record or not _supports_compiled_sphingolipid_gate(record)
            # A str subclass can compare equal to an earlier name while having
            # different conversion behavior. Keep this admission check fresh.
            or any(type(f.name) is not str for f in record.fragments)):
        return None
    from .equivalent_candidates import _sphingo_revision, _REFERENCE_SPHINGO_REVISION
    from .sphingolipid_rules import _series_key
    from .native_bridge import load_bridge
    last_key=(_evaluate_compiled_sphingolipid_gate,series)
    last=current[1].get(last_key)
    revision=None
    prepared=None
    if last is not None and last[0]==id(rule):
        bridge=load_bridge()
        unchanged=(bridge.sphingo_rule_unchanged(rule,last[1])
                   if bridge is not None and hasattr(bridge,'sphingo_rule_unchanged')
                   and _sphingo_revision is _REFERENCE_SPHINGO_REVISION else None)
        if unchanged is True:
            prepared=last[2]
        elif unchanged is None:
            revision=_sphingo_revision(rule)
            if revision==last[1]:
                prepared=last[2]
    if prepared is None:
        if revision is None:
            revision=_sphingo_revision(rule)
        key=(_evaluate_compiled_sphingolipid_gate,revision,series)
        prepared=current[1].get(key)
    if prepared is None:
        from .search_plan import _sphingo_program
        (types, library_types, adduct, class_key, tiered, identity_valid, hg_total) = _sphingolipid_library_evidence(record)
        groups = list(rule.required_any_groups)
        groups.extend(names for name,names in rule.required_any_by_series.items() if name==_series_key(series))
        exact_rule = replace(rule,required_any_groups=groups,required_any_by_series={},
                             allow_hg_only_fallback=False,allow_fah_only_fallback=False)
        program = _sphingo_program(record,exact_rule)
        routes = {id(f):(1 << i,types[id(f)]) for i,f in enumerate(record.fragments)}
        chain_mask = sum(1 << i for i,f in enumerate(record.fragments)
                         if types[id(f)] & {'Diagnostic_FA','Diagnostic_FA_Loss'})
        prepared = (program,routes,library_types,adduct,class_key,tiered,identity_valid,chain_mask)
        current[1][key] = prepared
    if revision is not None:
        current[1][last_key]=(id(rule),revision,prepared)
    program,routes,library_types,adduct,class_key,tiered,identity_valid,chain_mask = prepared
    mask, matched_types = 0,set()
    for match in matches:
        found = routes.get(id(match.fragment))
        if found is None or mask & found[0]:
            return None
        mask |= found[0]
        matched_types.update(found[1])
    passed = identity_valid and any(all((mask & bits).bit_count()>=minimum for bits,minimum in arm) for arm in program)
    return SphingolipidGateEvidence(bool(passed),frozenset(library_types),frozenset(matched_types),
                                   adduct,class_key,tiered,bool(tiered and mask & chain_mask))
