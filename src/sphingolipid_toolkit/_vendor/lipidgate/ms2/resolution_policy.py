from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

from .chain_utils import (
    FRAGMENT_CHAIN_TOKEN_RE,
    canonical_chain_token,
    chain_token_multiplicity,
    count_chain_evidence,
    extract_chain_tokens,
    extract_fragment_chain_token,
)
from .models import FragmentMatch, FragmentRecord, LibraryRecord
from .record_facts import cached_record_fact, cached_fragment_fact, cached_match_fact


CHAIN_LEVEL_INFO_MISSING_REASON = "missing_chain_level_information"
PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON = "insufficient_independent_chain_evidence"
SINGLE_CHAIN_COVERAGE_REASON = "multiple_fragments_single_chain_coverage"
POSITIVE_DIACYL_CHAIN_CONFIRMATION_CLASSES = frozenset({
    "PC", "PE", "PG", "PS", "PA", "PI", "PDPT", "MMPE", "DMPE", "PSC",
    "PMEOH", "PETOH", "PHEG", "CDPDG", "CEPE", "AMPS", "PIP", "PIP2", "PIP3",
    "OXPC", "OXPE",
})


@dataclass(frozen=True)
class FragmentEvidenceRule:
    fragment_types: frozenset[str]
    name_pattern: str = ""

    def matches(self, fragment: FragmentRecord) -> bool:
        if fragment.fragment_type not in self.fragment_types:
            return False
        if not self.name_pattern:
            return True
        return re.search(self.name_pattern, _compact_name(fragment.name), flags=re.IGNORECASE) is not None


@dataclass(frozen=True)
class ChainEvidenceProfile:
    name: str
    class_keys: frozenset[str]
    polarity: str
    evidence_rules: tuple[FragmentEvidenceRule, ...]
    required_coverage: str
    species_fallback_without_chain: bool = True
    applies_without_library_evidence: bool = False


_RCO_RULE = FragmentEvidenceRule(
    frozenset({"FA_Frag"}),
    r"^(?:\(R=O\)\+|RCO\(|\[RCO\]\+)",
)
_POSITIVE_HG_KETENE_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-R=O-[^]]+\+H\]\+\(",
)
_POSITIVE_PS_HG_KETENE_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-R=O-C3H8(?:O6NP|NO6P)\+H\]\+\(",
)
_POSITIVE_PC_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-\(?R(?:OOH|=O)\)?\+H\]\+\(",
)
_POSITIVE_PC_SODIUM_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}), r"^\[M\+NA-(?:59-)?(?:NA)?FA\]\+\(",
)
_POSITIVE_PE_SODIUM_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}), r"^\[M\+NA-(?:43-)?R(?:[12])?COOH\]\+\(",
)
_POSITIVE_PE_SODIUM_PROFILE = ChainEvidenceProfile(
    name="positive_pe_sodium", class_keys=frozenset({"PE"}), polarity="+",
    evidence_rules=(_RCO_RULE, _POSITIVE_PE_SODIUM_LOSS_RULE), required_coverage="infer_remaining_chain",
)
_POSITIVE_PC_P_VINYL_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}), r"^\[M-\(RCH2=CH-OH\)\+H\]\+\(P-",
)
_POSITIVE_LCB_RULE = FragmentEvidenceRule(frozenset({"LCB碎片"}))
_POSITIVE_LCB_LEGACY_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA"}),
    r"LCB",
)
_NEGATIVE_RCOO_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA"}),
    r"RCOO",
)
_NEGATIVE_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss", "Neutral_Loss"}),
    r"(?:R=O|ROOH)",
)
_NEGATIVE_CH3_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss", "Neutral_Loss"}),
    r"CH3.*(?:R=O|ROOH)|(?:R=O|ROOH).*CH3",
)


CHAIN_EVIDENCE_PROFILES: tuple[ChainEvidenceProfile, ...] = (
    ChainEvidenceProfile(
        name="positive_glycerophospholipid",
        class_keys=frozenset({"PA", "PE", "PG", "PI"}),
        polarity="+",
        evidence_rules=(_RCO_RULE, _POSITIVE_HG_KETENE_RULE),
        required_coverage="infer_remaining_chain",
    ),
    ChainEvidenceProfile(
        name="positive_ps",
        class_keys=frozenset({"PS"}),
        polarity="+",
        evidence_rules=(_POSITIVE_PS_HG_KETENE_RULE,),
        required_coverage="one",
        applies_without_library_evidence=True,
    ),
    ChainEvidenceProfile(
        name="positive_pc",
        class_keys=frozenset({"PC", "PCO"}),
        polarity="+",
        evidence_rules=(_POSITIVE_PC_CHAIN_LOSS_RULE, _POSITIVE_PC_SODIUM_LOSS_RULE),
        required_coverage="infer_remaining_chain",
    ),
    ChainEvidenceProfile(
        name="positive_pc_p", class_keys=frozenset({"PCP"}), polarity="+",
        evidence_rules=(_POSITIVE_PC_CHAIN_LOSS_RULE, _POSITIVE_PC_P_VINYL_LOSS_RULE),
        required_coverage="one",
    ),
    ChainEvidenceProfile(
        name="positive_sm",
        class_keys=frozenset({"SM", "LSM"}),
        polarity="+",
        evidence_rules=(_POSITIVE_LCB_RULE, _POSITIVE_LCB_LEGACY_RULE),
        required_coverage="one",
    ),
    ChainEvidenceProfile(
        name="negative_glycerophospholipid",
        class_keys=frozenset({"PA", "PE", "PG", "PI", "PS"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
    ChainEvidenceProfile(
        name="negative_pc",
        class_keys=frozenset({"PC", "PCO", "PCP"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CH3_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
    ChainEvidenceProfile(
        name="negative_sm",
        class_keys=frozenset({"SM", "LSM"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CH3_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
)


def normalized_class_key(compound_class: object) -> str:
    return _normalized_class_text(str(compound_class or ""))


@lru_cache(maxsize=256)
def _normalized_class_text(text: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", text.upper())


@cached_record_fact
def record_polarity(record: LibraryRecord) -> str:
    adduct = str(record.adduct or "").strip()
    if adduct.endswith("+"):
        return "+"
    if adduct.endswith("-"):
        return "-"
    return "+" if str(record.polarity or "").strip().startswith("+") else "-"


@cached_record_fact
def get_chain_evidence_profile(record: LibraryRecord) -> ChainEvidenceProfile | None:
    class_key = normalized_class_key(record.compound_class)
    if class_key == "PE" and str(record.adduct).strip() == "[M+Na]+":
        return _POSITIVE_PE_SODIUM_PROFILE
    polarity = record_polarity(record)
    return next(
        (
            profile
            for profile in CHAIN_EVIDENCE_PROFILES
            if profile.polarity == polarity and class_key in profile.class_keys
        ),
        None,
    )


@cached_fragment_fact
def fragment_is_chain_evidence(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    profile = get_chain_evidence_profile(record)
    if profile is None:
        return fragment.fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss", "FA_Frag"}
    return any(rule.matches(fragment) for rule in profile.evidence_rules)


def _compact_name(name: object) -> str:
    return _compact_name_text(str(name or ""))


@lru_cache(maxsize=4096)
def _compact_name_text(text: str) -> str:
    return re.sub(r"\s+", "", text).upper()


@cached_match_fact
def positive_phospholipid_chain_resolution(
    record: LibraryRecord, matches: Sequence[FragmentMatch],
) -> tuple[str, str] | None:
    """Separate the two-physical-peak gate from coverage confidence.

    Two losses of one chain retain chain-level reporting at low confidence.
    HG-only candidates and explicit ether/lyso exceptions keep their policies.
    """
    if (normalized_class_key(record.compound_class) not in POSITIVE_DIACYL_CHAIN_CONFIRMATION_CLASSES
            or record_polarity(record) != "+"):
        return None
    tokens = extract_chain_tokens(record.lipid_chain_name)
    if len(tokens) != 2 or any(token.upper().startswith(("O-", "P-")) for token in tokens):
        return None
    chain_matches = [match for match in matches if fragment_is_chain_evidence(record, match.fragment)]
    if not chain_matches:
        return None
    expected = {canonical_chain_token(match) for match in FRAGMENT_CHAIN_TOKEN_RE.finditer(record.lipid_chain_name)}
    observed = {
        extract_fragment_chain_token(match.fragment)
        for match in chain_matches
    }
    peaks = {
        float(match.experimental_peak.mz)
        for match in chain_matches
    }
    # Symmetric PE sodium carries two identical chain positions. One matching
    # chain ion covers both, without duplicating its measured peak/intensity.
    if (normalized_class_key(record.compound_class) == "PE"
            and str(record.adduct).strip() == "[M+Na]+"
            and len(expected) == 1 and expected.issubset(observed)
            and chain_evidence_count(record, matches) >= 2):
        return "chain_level", ""
    if len(peaks) < 2:
        return "tentative_species_level", PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
    if len(expected.intersection(observed)) < 2 or not expected.issubset(observed):
        return "chain_level", SINGLE_CHAIN_COVERAGE_REASON
    return "chain_level", ""


def positive_phospholipid_chain_confirmation_missing(record, matches) -> bool:
    decision = positive_phospholipid_chain_resolution(record, matches)
    return decision is not None and decision[1] == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON


def chain_evidence_count(
    record: LibraryRecord,
    fragments_or_matches: Sequence[FragmentRecord | FragmentMatch],
) -> int:
    return count_chain_evidence(record, fragments_or_matches, fragment_is_chain_evidence)


def required_chain_evidence_hits(record: LibraryRecord, profile: ChainEvidenceProfile) -> int:
    chain_count = len(extract_chain_tokens(record.lipid_chain_name))
    available_count = chain_evidence_count(record, record.fragments)
    if available_count <= 0:
        return 0
    if profile.required_coverage == "one":
        return 1
    if profile.required_coverage == "infer_remaining_chain":
        return min(max(chain_count - 1, 1), available_count)
    if profile.required_coverage == "all_observable_chains":
        return min(max(chain_count, 1), available_count)
    raise ValueError(f"Unsupported chain evidence coverage: {profile.required_coverage}")


@cached_record_fact
def _library_profile_required_hits(record):
    profile = get_chain_evidence_profile(record)
    return 0 if profile is None else required_chain_evidence_hits(record, profile)


def profile_chain_resolution(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> tuple[str, str] | None:
    """Return the configured chain/species decision, or None if not applicable."""

    profile = get_chain_evidence_profile(record)
    if profile is None:
        return None
    decision = positive_phospholipid_chain_resolution(record, matches)
    if decision is not None:
        return decision
    required_hits = _library_profile_required_hits(record)
    if required_hits <= 0:
        if profile.applies_without_library_evidence and profile.species_fallback_without_chain:
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        return None
    if chain_evidence_count(record, matches) >= required_hits:
        return "chain_level", ""
    if profile.species_fallback_without_chain:
        return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
    return None
