"""Add chain identities only where positional labels or masses identify them."""

import re
from dataclasses import replace

from .chain_utils import (
    FRAGMENT_CHAIN_TOKEN_RE,
    extract_chain_tokens,
    extract_fragment_chain_token,
)
from .esterified_ceramide import PROTON, WATER, C, H, O


def annotate_chain_label(
    lipid_name: str, label: str, precursor_mz: float, mz: float, adduct: str
) -> str:
    if extract_fragment_chain_token(label) is not None:
        return label
    if not re.search(r"R[1-4]|ROOH|R=O|RCO", label, re.IGNORECASE):
        return label
    tokens = extract_chain_tokens(lipid_name)
    # Do not reinterpret sphingoid bases or N-acyl/FAHFA hierarchies as R1/R2.
    if not tokens or any(not re.match(r"^(?:O-|P-)?\d+:\d+", t) for t in tokens):
        return label
    if re.search(r"-N-|\(FA\s|NAPE|NAPS|LNAPE|HBMP", lipid_name, re.IGNORECASE):
        return label
    if re.search(r"(?<!\d)0:0", lipid_name):
        return label
    detailed = [
        match.group(0) for match in FRAGMENT_CHAIN_TOKEN_RE.finditer(lipid_name)
    ]
    if len(detailed) != len(tokens):
        return label
    tokens = [max(pair, key=len) for pair in zip(tokens, detailed)]
    positions = set(
        re.findall(r"R([1-4])(?=COOH|COO|CO|=O|OH|[-+\)])", label, re.IGNORECASE)
    )
    if (
        len(positions) == 1
        and len(tokens) == 2
        and not re.search(r"(?<!\d)0:0", lipid_name)
    ):
        index = int(next(iter(positions))) - 1
        if index < len(tokens):
            return f"{label}({tokens[index]})"
    if positions:
        return label
    # Generic R losses can be assigned only for these exact simple pathways.
    compact = re.sub(r"\s+", "", label)
    candidates = set()
    for token in tokens:
        if not re.fullmatch(r"\d+:\d+", token):
            continue
        carbon, db = map(int, token.split(":"))
        fa = carbon * C + (2 * carbon - 2 * db) * H + 2 * O
        expected = None
        if compact in {"[M-(ROOH)+H]+", "[M-ROOH+H]+"} and adduct == "[M+H]+":
            expected = precursor_mz - fa
        elif compact in {"[M-(R=O)+H]+", "[M-R=O+H]+"} and adduct == "[M+H]+":
            expected = precursor_mz - fa + WATER
        elif compact in {"[RCOO]-", "RCOO-"}:
            expected = fa - PROTON
        elif compact in {"(R=O)+", "[RCO]+", "RCO+"}:
            expected = fa - WATER + PROTON
        if expected is not None and abs(expected - mz) <= 0.002:
            candidates.add(token)
    return f"{label}({next(iter(candidates))})" if len(candidates) == 1 else label


def annotate_record_chain_labels(record):
    return replace(
        record,
        fragments=[
            replace(
                f,
                name=annotate_chain_label(
                    record.lipid_chain_name,
                    f.name,
                    record.precursor_mz,
                    f.mz,
                    record.adduct,
                ),
            )
            for f in record.fragments
        ],
    )
