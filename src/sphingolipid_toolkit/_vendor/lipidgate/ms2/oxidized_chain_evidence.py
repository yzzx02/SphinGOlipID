"""Chain-specific count and localization gates for two-chain oxidized lipids."""

import re
from itertools import combinations

from .chain_utils import (
    FRAGMENT_CHAIN_TOKEN_RE,
    canonical_chain_token,
    extract_fragment_chain_token,
)


def oxygen_count(token):
    match = re.search(r";(?:O(\d+)|(\d+)OH)$", token)
    return int(match.group(1) or match.group(2)) if match else 0


def oxidized_chain_requirements(record):
    if not str(record.compound_class).upper().startswith("OX"):
        return None
    tokens = [
        canonical_chain_token(m)
        for m in FRAGMENT_CHAIN_TOKEN_RE.finditer(record.lipid_chain_name)
    ]
    if len(tokens) != 2 or not any(oxygen_count(t) for t in tokens):
        return None
    return tuple((token, 2 if oxygen_count(token) >= 2 else 1) for token in tokens)


def chain_evidence_groups(record, matches):
    requirements = oxidized_chain_requirements(record)
    if requirements is None:
        return {}
    groups = {token: {} for token, _ in requirements}
    for match in matches:
        fragment = match.fragment
        # Negative OxPC stores chain-labelled acid/ketene losses as Neutral_Loss.
        # They remain in other for scoring but can locate a chain.
        chain_loss = (
            fragment.fragment_type == "Neutral_Loss"
            and re.search(r"^\[M.*(?:R\d?(?:COOH|=O)|ROOH)", fragment.name)
        )
        if fragment.fragment_type not in {
            "Diagnostic_FA",
            "Diagnostic_FA_Loss",
            "FA_Frag",
        } and not chain_loss:
            continue
        token = extract_fragment_chain_token(fragment)
        if token not in groups:
            continue
        # A dehydrated oxidized FA can be exactly isobaric with an ordinary FA.
        # It supports a second observation but cannot alone localize oxygen.
        dehydrated = (
            re.search(r"-(?:\d*)H2O", fragment.name.replace(" ", ""), re.IGNORECASE)
            is not None
        )
        anchor = not dehydrated
        peak = float(match.experimental_peak.mz)
        current = groups[token].get(peak)
        if current is None or (anchor and not current[1]):
            groups[token][peak] = (match, anchor)
    return groups


def oxidized_chain_missing_groups(record, matches):
    requirements = oxidized_chain_requirements(record)
    if requirements is None:
        return []
    groups = chain_evidence_groups(record, matches)
    missing = []
    options = []
    for token, required in requirements:
        peaks = groups[token]
        oxidized = oxygen_count(token) > 0
        if len(peaks) < required:
            missing.append("oxidized_chain_fragments" if oxidized else "partner_chain")
        if oxidized and not any(anchor for _, anchor in peaks.values()):
            missing.append("oxidation_localization")
        options.append(
            [
                set(selected)
                for selected in combinations(peaks, required)
                if not oxidized or any(peaks[p][1] for p in selected)
            ]
        )
    if not missing and not any(a.isdisjoint(b) for a in options[0] for b in options[1]):
        missing.append("independent_chain_peaks")
    return list(dict.fromkeys(missing))
