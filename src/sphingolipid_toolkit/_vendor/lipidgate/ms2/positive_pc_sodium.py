"""Sodiated diacyl PC library: neutral-loss masses calculated from formulas."""

import re

from .esterified_ceramide import ELECTRON, C, H, N, O
from .models import FragmentRecord

P = 30.97376199842
NA = 22.9897692820
SODIUM_ION = NA - ELECTRON
TMA = 3 * C + 9 * H + N
PHOSPHOCHOLINE = 5 * C + 14 * H + N + 4 * O + P
SODIUM_PHOSPHOCHOLINE = PHOSPHOCHOLINE - H + NA
HG_NAMES = frozenset({"[M+Na-59]+", "[M+Na-183]+", "[M+Na-205]+"})


def is_positive_pc_sodium(record):
    return record.compound_class == "PC" and record.adduct == "[M+Na]+"


def sodium_pc_fragments(chain_name):
    match = re.fullmatch(r"PC(?:\s+|\()(\d+):(\d+)[_/](\d+):(\d+)\)?", chain_name)
    if not match:
        raise ValueError(f"Not a diacyl PC: {chain_name}")
    c1, db1, c2, db2 = map(int, match.groups())
    carbon = c1 + c2 + 8
    hydrogen = 2 * (c1 + c2) - 2 * (db1 + db2) + 16
    formula = f"C{carbon}H{hydrogen}NO8P"
    precursor = carbon * C + hydrogen * H + N + 8 * O + P + SODIUM_ION
    fragments = [
        FragmentRecord(precursor, "[M+Na]+", "Precursor Ion"),
        FragmentRecord(precursor - TMA, "[M+Na-59]+", "Diagnostic_HG"),
        FragmentRecord(precursor - PHOSPHOCHOLINE, "[M+Na-183]+", "Diagnostic_HG"),
        FragmentRecord(
            precursor - SODIUM_PHOSPHOCHOLINE, "[M+Na-205]+", "Diagnostic_HG"
        ),
        FragmentRecord(
            5 * C + 15 * H + N + 4 * O + P - ELECTRON, "[C5H15NO4P]+", "Common"
        ),
        FragmentRecord(
            2 * C + 5 * H + 4 * O + P + SODIUM_ION, "[C2H5O4PNa]+", "Common"
        ),
        FragmentRecord(5 * C + 12 * H + N - ELECTRON, "[C5H12N]+", "Common"),
    ]
    for c, db in sorted({(c1, db1), (c2, db2)}):
        fa = c * C + (2 * c - 2 * db) * H + 2 * O
        for loss, name in [
            (fa, "[M+Na-FA]+"),
            (fa - H + NA, "[M+Na-NaFA]+"),
            (TMA + fa, "[M+Na-59-FA]+"),
        ]:
            fragments.append(
                FragmentRecord(
                    precursor - loss, f"{name}({c}:{db})", "Diagnostic_FA_Loss"
                )
            )
    return precursor, formula, sorted(fragments, key=lambda f: f.mz)
