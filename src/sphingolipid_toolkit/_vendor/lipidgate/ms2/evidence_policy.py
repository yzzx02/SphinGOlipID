"""Fragment routing and read-only evidence facts; no scores or eligibility decisions."""
from __future__ import annotations

import re
from functools import lru_cache
from typing import Callable, List, Sequence
from .chain_utils import chain_token_counts as _chain_token_counts, extract_fragment_chain_token as _extract_fragment_chain_token, count_chain_evidence
from .esterified_ceramide import is_esterified_ceramide
from .positive_cer_sodium import is_positive_cer_sodium, GATE_NAME as CER_SODIUM_GATE
from .positive_cer1p import is_positive_cer1p, GATE_NAMES as CER1P_GATE_NAMES
from .positive_complex_ceramide import is_positive_complex_ceramide
from .models import FragmentMatch, FragmentRecord, LibraryRecord
from .negative_cer_hsu2016 import is_simple_negative_cer, logical_negative_cer_fragment_types
from .negative_gm3 import is_negative_gm3
from .negative_lpe import HG_NAMES as NEGATIVE_LPE_HG_NAMES, is_negative_lpe
from .negative_lpc import HG_NAME as NEGATIVE_LPC_HG_NAME, is_negative_lpc
from .positive_lpc_ether import HG_NAME as POSITIVE_LPC_HG_NAME, is_positive_ether_lpc
from .negative_hexcer import is_negative_hexcer, logical_fragment_types
from .oxidized_fragment_policy import parent_class_key, is_oxidized_precursor_water_loss
from .positive_hexosylceramide import (
    HG_TYPES as HEXOSYL_HG_TYPES, LCB_TYPE as HEXOSYL_LCB_TYPE,
    fragment_type as positive_hexosyl_fragment_type, is_positive_hexosylceramide,
)
from .positive_gm3 import B_FRAGMENT_TYPE
from .positive_pe_cer import is_positive_pe_cer
from .positive_spb import is_positive_spb, fragment_type as positive_spb_fragment_type, LCB_TYPE as SPB_LCB_TYPE
from .record_facts import cached_fragment_fact, cached_record_fact, cached_match_fact
from .resolution_policy import (
    normalized_class_key,
    record_polarity,
    fragment_is_chain_evidence,
    get_chain_evidence_profile,
)


LOSS_FRAGMENT_TYPES = {"Neutral_Loss", "Diagnostic_FA_Loss"}
SPHINGOLIPID_LCB_FRAGMENT_TYPE = "LCB碎片"
NEGATIVE_CER_CHAIN_FRAGMENT_TYPES = {"LCB碎片", "FA类碎片", "NAE碎片"}
NEGATIVE_FA_LOSS_SUPPORT_ONLY_CLASSES = {
    "PA",
    "PE",
    "OXPE",
    "PG",
    "PI",
    "PS",
    "PAO",
    "PEO",
    "PGO",
    "PIO",
    "PSO",
    "PEP",
    "DMPE",
    "MMPE",
    "PMEOH",
    "PETOH",
    "PGCN",
    "OXPGCN",
    "LPE",
    "LPETOH",
    "PLASMENYLPNE",
    "PNEP",
}
FREE_SPHINGOID_BASE_CLASS_KEYS = {"SPB", "SPH", "DHSPH", "PHYTOSPH"}
SPB_UNASSIGNED_DIAGNOSTIC_NAMES = {
    "SPB-Diagnostic-1",
    "SPB-Diagnostic-2",
}
POSITIVE_FA_FRAG_AS_LOSS_CLASSES = {"PA", "PE", "PG", "PI"}
POSITIVE_GLYCERIDE_RCO_GATE_CLASSES = {"TG", "DG", "TGO", "OXTG"}
CANDIDATE_HG_FRAGMENT_TYPE = "Candidate_HG"
POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES = {
    "PA",
    "PC",
    "PCO",
    "PCP",
    "PE",
    "PG",
    "PI",
    "PS",
}
POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES = {
    "PA",
    "PC",
    "PCO",
    "PCP",
    "PE",
    "PEO",
    "PEP",
    "PG",
    "PGO",
    "PI",
    "PIO",
    "PS",
    "PSO",
}
POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES = {
    "BMP",
    "CL",
    "DLCL",
    "DMPE",
    "ETHERLPG",
    "ETHERLPI",
    "ETHERPG",
    "HBMP",
    "LNAPE",
    "MLCL",
    "NAPE",
    "PETOH",
    "PIP",
    "PIP2",
    "PIP3",
    "PMEOH",
}
POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES = {
    "S1P",
    "ASM",
    "CERP",
    "CER1P",
    "DHSPH",
    "LSM",
    "PHYTOSPH",
    "SL",
    "SM",
    "SPB",
    "SPH",
}
CLASS_SPECIFIC_HG_MZ = {
    "PG": (152.9933, 171.0064, 209.0221, 227.0326),
    "PETOH": (181.0271, 181.0280),
    "PMEOH": (167.0109,),
    "DMPE": (168.0411, 168.0431),
}
PRECURSOR_FRAGMENT_REQUIRED_CLASSES = {"NAASP"}


@lru_cache(maxsize=256, typed=True)
def _is_positive_adduct(adduct: str) -> bool:
    return str(adduct or "").strip().endswith("+")


@lru_cache(maxsize=256, typed=True)
def _normalized_compound_class(compound_class: str) -> str:
    return str(compound_class or "").strip().upper()


@lru_cache(maxsize=256, typed=True)
def _normalized_class_key(compound_class: str) -> str:
    return _normalized_compound_class(compound_class).replace("-", "")


def _is_fa_record(record: LibraryRecord) -> bool:
    return _normalized_compound_class(record.compound_class) == "FA"


@cached_record_fact
def _is_positive_hg_dominant_phospholipid(record: LibraryRecord) -> bool:
    return (
        _is_positive_adduct(record.adduct)
        and parent_class_key(record.compound_class) in POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES
    )


@cached_record_fact
def _uses_positive_support_pool_scoring(record: LibraryRecord) -> bool:
    if not _is_positive_adduct(record.adduct):
        return False
    class_key = _normalized_class_key(record.compound_class)
    parent_key = parent_class_key(record.compound_class)
    phospholipid_key = parent_key[1:] if parent_key.startswith("L") else parent_key
    return (
        phospholipid_key in POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES
        or class_key in POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES
        or class_key in POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES
        or "CER" in class_key
        or class_key.startswith(("GM", "GD", "GT", "GQ"))
    )


@cached_record_fact
def _is_positive_glyceride_rco_gate_record(record: LibraryRecord) -> bool:
    cls = _normalized_class_key(record.compound_class)
    return _is_positive_adduct(record.adduct) and cls in POSITIVE_GLYCERIDE_RCO_GATE_CLASSES


def _is_positive_oxtg_full_fah_record(record: LibraryRecord) -> bool:
    return (
        _is_positive_adduct(record.adduct)
        and _normalized_class_key(record.compound_class) == "OXTG"
    )


def _is_rco_fragment_name(fragment_name: str) -> bool:
    name = str(fragment_name or "").strip().upper()
    return name.startswith("(R=O)+(") or name.startswith("RCO(") or name.startswith("[RCO]+")


def _is_positive_glyceride_rco_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    return (
        _is_positive_glyceride_rco_gate_record(record)
        and fragment.fragment_type == "FA_Frag"
        and _is_rco_fragment_name(fragment.name)
    )


def _is_positive_glyceride_rco_c3h6o2_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if not _is_positive_glyceride_rco_gate_record(record):
        return False
    if fragment.fragment_type != "Diagnostic_FA_Loss":
        return False
    name = str(fragment.name or "").strip().upper()
    return "C3H6O2" in name and "ROOH" not in name and "M-NH3" not in name


def _is_positive_glyceride_chain_gate_fragment(
    record: LibraryRecord,
    fragment: FragmentRecord,
) -> bool:
    return _is_positive_glyceride_rco_fragment(record, fragment) and _normalized_class_key(record.compound_class) != "DG"


def _positive_glyceride_chain_gate_group(record: LibraryRecord) -> str:
    return "rco"


@cached_record_fact
def _is_positive_mg_record(record: LibraryRecord) -> bool:
    return _is_positive_adduct(record.adduct) and _normalized_compound_class(record.compound_class) == "MG"


def _is_class_specific_hg_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    cls = _normalized_class_key(record.compound_class)
    if cls == "MG" and _is_positive_adduct(record.adduct):
        return str(fragment.name or "").strip() == "[M-H2O+H]+"
    if fragment.fragment_type not in {"Common", CANDIDATE_HG_FRAGMENT_TYPE}:
        return False
    if _is_positive_adduct(record.adduct):
        return False
    return any(abs(float(fragment.mz) - target_mz) <= 0.02 for target_mz in CLASS_SPECIFIC_HG_MZ.get(cls, ()))


def _record_precursor_ion_fragments(record: LibraryRecord) -> List[FragmentRecord]:
    return [fragment for fragment in record.fragments if fragment.fragment_type == "Precursor Ion"]


def _requires_precursor_fragment(record: LibraryRecord) -> bool:
    class_key = _normalized_class_key(record.compound_class)
    return (
        is_negative_lpe(record)
        or _normalized_compound_class(record.compound_class) in PRECURSOR_FRAGMENT_REQUIRED_CLASSES
        or (
            class_key == "ADGGA"
            and not _is_positive_adduct(record.adduct)
        )
    )


def _matched_precursor_ion_fragment_count(matches: Sequence[FragmentMatch]) -> int:
    return sum(1 for match in matches if match.fragment.fragment_type == "Precursor Ion")


@cached_record_fact
def _fa_frag_counts_as_effective_loss(record: LibraryRecord) -> bool:
    return (
        _is_positive_adduct(record.adduct)
        and parent_class_key(record.compound_class) in POSITIVE_FA_FRAG_AS_LOSS_CLASSES
    )


@cached_fragment_fact
def _fragment_counts_as_effective_loss(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if _is_positive_ps_record(record):
        return (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and _is_positive_ps_headgroup_ketene_loss(fragment)
        )
    if _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment):
        return False
    return (
        fragment.fragment_type in LOSS_FRAGMENT_TYPES
        or (_fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag")
    )


@cached_fragment_fact
def _fragment_counts_as_fa_loss_gate(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if _is_positive_ps_record(record):
        return (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and _is_positive_ps_headgroup_ketene_loss(fragment)
        )
    return (
        (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and not _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment)
        )
        or (_fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag")
    )


def _fragment_counts_as_chain_resolving_loss(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if get_chain_evidence_profile(record) is not None:
        return fragment_is_chain_evidence(record, fragment)
    return _fragment_counts_as_effective_loss(record, fragment)


@cached_record_fact
def _record_has_diagnostic_hg(record: LibraryRecord) -> bool:
    return any(
        fragment.fragment_type == "Diagnostic_HG"
        or _is_class_specific_hg_fragment(record, fragment)
        for fragment in record.fragments
    )


@cached_record_fact
def _record_has_candidate_hg(record: LibraryRecord) -> bool:
    return any(fragment.fragment_type == CANDIDATE_HG_FRAGMENT_TYPE for fragment in record.fragments)


@cached_fragment_fact
def _fragment_counts_as_hg(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if is_positive_spb(record):
        return False
    if is_positive_hexosylceramide(record.compound_class, record.adduct):
        return positive_hexosyl_fragment_type(record.compound_class, record.lipid_chain_name, fragment) in HEXOSYL_HG_TYPES
    if is_positive_ether_lpc(record):
        return fragment.fragment_type == "Diagnostic_HG" and fragment.name == POSITIVE_LPC_HG_NAME
    if is_negative_lpc(record):
        return fragment.fragment_type == "Diagnostic_HG" and fragment.name == NEGATIVE_LPC_HG_NAME
    if is_negative_lpe(record):
        return fragment.fragment_type == "Diagnostic_HG" and fragment.name in NEGATIVE_LPE_HG_NAMES
    if _is_positive_ps_record(record):
        return _is_positive_ps_headgroup_loss(fragment)
    if _is_positive_mg_record(record):
        return _is_class_specific_hg_fragment(record, fragment)
    if fragment.fragment_type == "Diagnostic_HG":
        return True
    if _is_class_specific_hg_fragment(record, fragment):
        return True
    if _record_has_diagnostic_hg(record):
        return False
    if fragment.fragment_type == CANDIDATE_HG_FRAGMENT_TYPE:
        return True
    return fragment.fragment_type == "Precursor Ion" and _record_has_candidate_hg(record)


def _pool_for_fragment(record: LibraryRecord, fragment: FragmentRecord) -> str:
    return _pool_for_scoring_fragment(record, fragment, matched=False)


@cached_record_fact
def _is_negative_ahexcer_record(record: LibraryRecord) -> bool:
    return (
        not _is_positive_adduct(record.adduct)
        and _normalized_class_key(record.compound_class) == "AHEXCER"
    )


def _is_ahexcer_overlapping_lcb_loss(fragment: FragmentRecord) -> bool:
    return (
        fragment.fragment_type == "Diagnostic_FA_Loss"
        and "LCB-C2H7NO" in str(fragment.name or "").upper()
    )


@cached_fragment_fact(matched_independent=True)
def _scoring_pools_for_fragment(
    record: LibraryRecord,
    fragment: FragmentRecord,
    *,
    matched: bool,
) -> tuple[str, ...]:
    primary_pool = _pool_for_scoring_fragment(record, fragment, matched=matched)
    if (
        is_negative_hexcer(record.compound_class, record.adduct)
        and fragment.fragment_type == "Diagnostic_FA"
        and SPHINGOLIPID_LCB_FRAGMENT_TYPE in logical_fragment_types(fragment)
    ):
        return ("fah", "lcb")
    if _is_negative_ahexcer_record(record) and _is_ahexcer_overlapping_lcb_loss(fragment):
        # The two AHexCer losses are both outer-acyl losses and the requested
        # LCB-loss evidence. One physical peak may satisfy both logical pools
        # without being duplicated by the peak matcher.
        return ("fah", "other")
    return (primary_pool,)


@cached_fragment_fact(matched_independent=True)
def _pool_for_scoring_fragment(record: LibraryRecord, fragment: FragmentRecord, *, matched: bool) -> str:
    if is_positive_complex_ceramide(record.compound_class,record.adduct):
        return 'hg' if fragment.fragment_type=='Diagnostic_HG' else 'lcb' if fragment.fragment_type==SPHINGOLIPID_LCB_FRAGMENT_TYPE else 'other'
    if is_positive_cer1p(record):
        return "hg" if fragment.name == CER1P_GATE_NAMES[record.adduct] else "other"
    if is_positive_cer_sodium(record):
        # Reuse the HG numeric slot for the precursor water-loss pool. This is
        # structural evidence, not a chemical headgroup or LCB assignment.
        return "hg" if fragment.name == CER_SODIUM_GATE else "other"
    if (_normalized_class_key(record.compound_class) == "CER" and str(record.adduct).strip() == "[M+H]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name)):
        role = positive_hexosyl_fragment_type("Cer", record.lipid_chain_name, fragment)
        return "lcb" if role == SPHINGOLIPID_LCB_FRAGMENT_TYPE else "other"
    if _normalized_class_key(record.compound_class) == "SM" and str(record.adduct).strip() == "[M+H]+":
        return "hg" if fragment.fragment_type == "Diagnostic_HG" else "other"
    if is_positive_spb(record):
        return "lcb" if positive_spb_fragment_type(fragment) == SPB_LCB_TYPE else "other"
    if is_positive_hexosylceramide(record.compound_class, record.adduct):
        role = positive_hexosyl_fragment_type(record.compound_class, record.lipid_chain_name, fragment)
        if role == HEXOSYL_LCB_TYPE:
            return "lcb"
        if role == B_FRAGMENT_TYPE:
            return "b"
        return "hg" if role in HEXOSYL_HG_TYPES else "other"
    if is_positive_ether_lpc(record):
        if fragment.fragment_type == "Precursor Ion" and abs(fragment.mz-record.precursor_mz) <= 1e-6:
            return "fah"  # Secondary precursor pool; no chain evidence or gate.
        return "hg" if _fragment_counts_as_hg(record, fragment) else "other"
    if is_negative_lpc(record):
        return "hg" if _fragment_counts_as_hg(record, fragment) else "other"
    if is_negative_lpe(record):
        if fragment.fragment_type == "Precursor Ion" and abs(fragment.mz-record.precursor_mz) <= 1e-6:
            return "fah"  # Existing primary-pool field; the evidence is precursor, not a chain ion.
        return "hg" if _fragment_counts_as_hg(record, fragment) else "other"
    if _normalized_class_key(record.compound_class).startswith("OX") and is_oxidized_precursor_water_loss(fragment):
        return "other"
    if (
        fragment.fragment_type == SPHINGOLIPID_LCB_FRAGMENT_TYPE
        and (
            (is_esterified_ceramide(record.compound_class, record.lipid_chain_name) and _is_positive_adduct(record.adduct))
            or is_negative_hexcer(record.compound_class, record.adduct)
            or is_positive_pe_cer(record.compound_class, record.adduct)
            or is_negative_gm3(record.compound_class, record.adduct)
        )
    ):
        return "lcb"
    if _is_positive_ps_record(record):
        if _is_positive_ps_headgroup_loss(fragment):
            return "hg"
        if (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and _is_positive_ps_headgroup_ketene_loss(fragment)
        ):
            return "fah"
        return "other"
    if (
        _normalized_class_key(record.compound_class) in {"FA", "OXFA"}
        and fragment.fragment_type == "Precursor Ion"
        and "RCOO" in str(fragment.name or "").upper()
    ):
        return "fah"
    if _is_positive_glyceride_rco_fragment(record, fragment):
        return "other"
    if _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment):
        return "other"
    if fragment.fragment_type == "Diagnostic_FA_Loss" and (
        not _is_positive_adduct(record.adduct)
        and parent_class_key(record.compound_class) in NEGATIVE_FA_LOSS_SUPPORT_ONLY_CLASSES
    ):
        return "other"
    if fragment.fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}:
        return "fah"
    if _fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag":
        return "fah"
    if _normalized_class_key(record.compound_class) in FREE_SPHINGOID_BASE_CLASS_KEYS:
        if fragment.fragment_type == "C类碎片":
            return "fah"
        if str(fragment.name or "").strip() in SPB_UNASSIGNED_DIAGNOSTIC_NAMES:
            return "other"
    if (
        _normalized_class_key(record.compound_class) == "CER"
        and not _is_positive_adduct(record.adduct)
        and fragment.fragment_type in NEGATIVE_CER_CHAIN_FRAGMENT_TYPES
    ):
        return "fah"
    if fragment.fragment_type == SPHINGOLIPID_LCB_FRAGMENT_TYPE:
        return "fah"
    if _fragment_counts_as_hg(record, fragment):
        return "hg"
    return "other"


def _chain_evidence_count_for_fragments(
    record: LibraryRecord,
    fragments: Sequence[FragmentRecord],
    predicate: Callable[[LibraryRecord, FragmentRecord], bool],
) -> int:
    return count_chain_evidence(record, fragments, predicate)


@cached_match_fact
def _chain_evidence_count_for_matches(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    predicate: Callable[[LibraryRecord, FragmentRecord], bool],
) -> int:
    return count_chain_evidence(record, matches, predicate)


def _record_expected_fah_tokens(record: LibraryRecord) -> List[str]:
    if is_negative_lpe(record) or is_negative_lpc(record) or is_positive_ether_lpc(record):
        return []
    tokens: List[str] = []
    seen = set()
    for fragment in record.fragments:
        if fragment.fragment_type != "Diagnostic_FA":
            continue
        token = _extract_fragment_chain_token(fragment.name)
        if token is None or token in seen:
            continue
        seen.add(token)
        tokens.append(token)
    return tokens


@cached_record_fact
def _record_hg_fragment_count(record: LibraryRecord) -> int:
    return sum(1 for fragment in record.fragments if _fragment_counts_as_hg(record, fragment))


@cached_match_fact
def _matched_hg_fragment_count(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> int:
    matched_hg_fragments = [
        match.fragment
        for match in matches
        if _fragment_counts_as_hg(record, match.fragment)
    ]
    if not matched_hg_fragments:
        return 0

    if (
        _is_positive_adduct(record.adduct)
        and _normalized_class_key(record.compound_class) == "BMP"
    ):
        # A repeated-chain BMP has two MAG library entries at the same m/z.
        # The peak matcher intentionally consumes an experimental peak once,
        # but that one MAG ion still identifies both identical chains.
        matched_mz = {float(fragment.mz) for fragment in matched_hg_fragments}
        return sum(
            1
            for fragment in record.fragments
            if _fragment_counts_as_hg(record, fragment)
            and any(abs(float(fragment.mz) - mz) <= 1e-6 for mz in matched_mz)
        )

    return len(matched_hg_fragments)


@cached_record_fact
def _record_loss_fragment_count(record: LibraryRecord) -> int:
    return sum(1 for fragment in record.fragments if _fragment_counts_as_effective_loss(record, fragment))


@cached_record_fact
def _record_fa_loss_fragment_count(record: LibraryRecord) -> int:
    if _is_positive_adduct(record.adduct):
        return _chain_evidence_count_for_fragments(record, record.fragments, _fragment_counts_as_fa_loss_gate)
    return sum(1 for fragment in record.fragments if _fragment_counts_as_fa_loss_gate(record, fragment))


@cached_record_fact
def _record_positive_glyceride_rco_fragment_count(record: LibraryRecord) -> int:
    if not _is_positive_glyceride_rco_gate_record(record):
        return 0
    return _chain_evidence_count_for_fragments(
        record,
        record.fragments,
        _is_positive_glyceride_chain_gate_fragment,
    )


@cached_match_fact
def _matched_positive_glyceride_rco_fragment_count(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> int:
    if not _is_positive_glyceride_rco_gate_record(record):
        return 0
    return _chain_evidence_count_for_matches(
        record,
        matches,
        _is_positive_glyceride_chain_gate_fragment,
    )


@cached_record_fact
def _record_positive_signature_fragment_count(record: LibraryRecord) -> int:
    return sum(
        1
        for fragment in record.fragments
        if not _fragment_counts_as_hg(record, fragment) and fragment.fragment_type != "Precursor Ion"
    )


def _count_positive_nonzero_chains(chain_tokens: Sequence[str]) -> int:
    return len([token for token in chain_tokens if token and token != "0:0"])


def _positive_required_loss_hits(chain_tokens: Sequence[str], fa_loss_fragment_count: int) -> int:
    chain_count = _count_positive_nonzero_chains(chain_tokens)
    expected = max(chain_count - 1, 1)
    if fa_loss_fragment_count > 0:
        expected = min(expected, fa_loss_fragment_count)
    return expected


@cached_match_fact
def _matched_hg_max_relative_intensity(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> float:
    return max(
        (
            float(match.experimental_peak.relative_intensity)
            for match in matches
            if _fragment_counts_as_hg(record, match.fragment)
        ),
        default=0.0,
    )


def _count_fah_expected_chains(chain_tokens: Sequence[str]) -> int:
    expected_tokens = set()
    for token in chain_tokens:
        if token == "0:0":
            continue
        if token.startswith("O-") or token.startswith("P-"):
            continue
        expected_tokens.add(token)
    return len(expected_tokens)


def _is_cl_double_negative_record(record: LibraryRecord) -> bool:
    return record.compound_class == "CL" and record.adduct == "[M-2H]2-"


def _is_cl_chain_info_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if not _is_cl_double_negative_record(record):
        return False
    if fragment.fragment_type in {"Diagnostic_FA", "Precursor Ion"}:
        return False
    name = str(fragment.name or "")
    if _extract_fragment_chain_token(name) is not None:
        return True
    return any(marker in name for marker in ("C3H5O4P+(R", "RCOO", "R=O", "ROOH", "RCH2"))


def _matched_cl_chain_info_fragment_count(record: LibraryRecord, matches: Sequence[FragmentMatch]) -> int:
    return sum(1 for match in matches if _is_cl_chain_info_fragment(record, match.fragment))


@cached_match_fact
def _matched_required_fah_count(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    expected_fah_tokens: Sequence[str],
) -> int:
    matched_fah = set(_matched_fah_tokens(matches))
    expected_fah = set(expected_fah_tokens)
    multiplicity = _chain_token_counts(record)
    if expected_fah:
        matched_expected = matched_fah.intersection(expected_fah)
        if not _is_positive_adduct(record.adduct):
            return len(matched_expected)
        if multiplicity:
            return sum(multiplicity.get(token, 1) for token in matched_expected)
        return len(matched_expected)
    if multiplicity:
        return sum(multiplicity.get(token, 1) for token in matched_fah)
    return len(matched_fah)


def _matched_fah_tokens(matches: Sequence[FragmentMatch]) -> List[str]:
    tokens = []
    for match in matches:
        if match.fragment.fragment_type != "Diagnostic_FA":
            continue
        token = _extract_fragment_chain_token(match.fragment.name)
        if token is not None:
            tokens.append(token)
    return tokens


def _matched_positive_signature_fragment_count(matches: Sequence[FragmentMatch]) -> int:
    return sum(
        1
        for match in matches
        if match.fragment.fragment_type not in {"Diagnostic_HG", CANDIDATE_HG_FRAGMENT_TYPE, "Precursor Ion"}
    )


def _positive_single_chain_species_hg_only(
    chain_tokens: Sequence[str],
    is_positive_mode: bool,
    allow_lyso_hg_only: bool,
    expected_fah_tokens: Sequence[str],
) -> bool:
    return (
        is_positive_mode
        and allow_lyso_hg_only
        and not expected_fah_tokens
        and _count_positive_nonzero_chains(chain_tokens) == 1
    )


def is_positive_class(record: LibraryRecord, class_keys: frozenset[str]) -> bool:
    return record_polarity(record) == "+" and normalized_class_key(
        record.compound_class
    ) in class_keys


@cached_record_fact
def is_positive_ps(record: LibraryRecord) -> bool:
    return is_positive_class(record, frozenset({"PS"}))


@cached_record_fact
def is_positive_pep(record: LibraryRecord) -> bool:
    return is_positive_class(record, frozenset({"PEP"}))


def is_positive_ps_headgroup_loss(fragment: FragmentRecord) -> bool:
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    return compact_name in {"[M-C3H8O6NP+H]+", "[M-C3H8NO6P+H]+"}


def is_positive_ps_chain_ketene_loss(fragment: FragmentRecord) -> bool:
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    return (
        fragment.fragment_type == "Diagnostic_FA_Loss"
        and (
            compact_name.startswith("[M-R=O-C3H8O6NP+H]+(")
            or compact_name.startswith("[M-R=O-C3H8NO6P+H]+(")
        )
    )


_is_positive_ps_record = is_positive_ps
_is_positive_ps_headgroup_loss = is_positive_ps_headgroup_loss
_is_positive_ps_headgroup_ketene_loss = is_positive_ps_chain_ketene_loss


def sphingolipid_fragment_types(record: LibraryRecord, fragment: FragmentRecord) -> set[str]:
    """Retain multiple logical roles without duplicating one physical peak."""
    if (record.compound_class == "Cer" and record.adduct == "[M+H]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name)):
        return {positive_hexosyl_fragment_type("Cer", record.lipid_chain_name, fragment)}
    if is_positive_spb(record):
        return {positive_spb_fragment_type(fragment)}
    if is_positive_hexosylceramide(record.compound_class, record.adduct):
        return {positive_hexosyl_fragment_type(record.compound_class, record.lipid_chain_name, fragment)}
    if is_negative_hexcer(record.compound_class, record.adduct):
        return logical_fragment_types(fragment)
    if record.compound_class == "Cer" and is_simple_negative_cer(record.lipid_chain_name, record.adduct):
        return logical_negative_cer_fragment_types(fragment)
    return {fragment.fragment_type}
