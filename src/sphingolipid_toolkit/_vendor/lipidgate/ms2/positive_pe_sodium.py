"""Formula-derived diacyl PE [M+Na]+ fragments.

The four acyl-loss channels follow Simoes et al., 2008, DOI 10.1002/rcm.3727.
Names encode chain composition, not experimentally established sn positions.
"""
from __future__ import annotations

import re
from .esterified_ceramide import C, H, N, O, ELECTRON
from .models import FragmentRecord
from .positive_pc_sodium import P, NA, SODIUM_ION

PE_HEADGROUP = 2 * C + 8 * H + N + 4 * O + P
VINYLAMINE = 2 * C + 5 * H + N
SODIUM_PE_HEADGROUP = PE_HEADGROUP - H + NA
SODIUM_PE_HEADGROUP_ION = PE_HEADGROUP + SODIUM_ION
HG_NAMES = frozenset({'[M+Na-43]+', '[M+Na-141]+'})


def is_positive_pe_sodium(record) -> bool:
    return str(record.compound_class).strip().upper() == 'PE' and str(record.adduct).strip() == '[M+Na]+'


def sodium_pe_fragments(chain_name: str) -> tuple[float, str, list[FragmentRecord]]:
    match = re.fullmatch(r"PE(?:\s+|\()(\d+):(\d+)[_/](\d+):(\d+)\)?", chain_name)
    if match is None:
        raise ValueError(f"Not a diacyl PE: {chain_name}")
    c1, db1, c2, db2 = map(int, match.groups())
    carbon = c1 + c2 + 5
    hydrogen = 2 * (c1 + c2) - 2 * (db1 + db2) + 10
    formula = f"C{carbon}H{hydrogen}NO8P"
    precursor = carbon * C + hydrogen * H + N + 8 * O + P + SODIUM_ION
    fragments = [
        FragmentRecord(precursor - PE_HEADGROUP, "[M+Na-141]+", "Diagnostic_HG"),
        FragmentRecord(precursor - VINYLAMINE, "[M+Na-43]+", "Diagnostic_HG"),
        FragmentRecord(precursor - SODIUM_PE_HEADGROUP, "[M+Na-163]+", "Common"),
        FragmentRecord(SODIUM_PE_HEADGROUP_ION, "[C2H8NO4PNa]+", "Common"),
        FragmentRecord(precursor, "[M+Na]+", "Precursor Ion"),
    ]
    # Equal chains give two physical loss ions, not four duplicate peaks.
    for carbon_count, db_count in sorted({(c1, db1), (c2, db2)}):
        acid = carbon_count * C + (2 * carbon_count - 2 * db_count) * H + 2 * O
        acylium = carbon_count * C + (2 * carbon_count - 2 * db_count - 1) * H + O - ELECTRON
        token = f"{carbon_count}:{db_count}"
        fragments.extend([
            FragmentRecord(acylium, f"(R=O)+({token})", "FA_Frag"),
            FragmentRecord(precursor - acid, f"[M+Na-RCOOH]+({token})", "Diagnostic_FA_Loss"),
            FragmentRecord(precursor - VINYLAMINE - acid,
                           f"[M+Na-43-RCOOH]+({token})", "Diagnostic_FA_Loss"),
        ])
    if any(fragment.mz <= 0 for fragment in fragments):
        raise ValueError(f"Invalid PE composition: {chain_name}")
    return precursor, formula, sorted(fragments, key=lambda fragment: fragment.mz)
