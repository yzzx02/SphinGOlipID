"""Negative Glc-GP and 6-O-acetyl-Glc-GP, with production MSP labels.

The acetylated class uses the user-supplied 18:0/20:0 measurements; its
deglycosylated HG is theoretical support, not a second mandatory HG hit.
"""
from __future__ import annotations

import re

from .esterified_ceramide import C, H, O, PROTON, WATER
from .models import FragmentRecord

P = 30.97376199842
GLYCOSYL_RESIDUE = 6 * C + 10 * H + 5 * O
ACETYL_INCREMENT = 2 * C + 2 * H + O
ACETYLGLYCOSYL_RESIDUE = GLYCOSYL_RESIDUE + ACETYL_INCREMENT
GLYCEROPHOSPHATE_ANION = 3 * C + 7 * H + 5 * O + P - PROTON
GLC_GP_CLASSES = frozenset({"PtdGlc", "Glc-GP", "6-Ac-Glc-GP"})
REFERENCE_TITLE = "Quantitative analysis of a novel glucosylated phospholipid by liquid chromatography–mass spectrometry"
REFERENCE_DOI = "10.1016/j.ab.2008.02.007"
REFERENCE_PMID = "18342611"


def canonicalize_glc_gp_identity(
    compound_class: str, lipid_name: str, lipid_chain_name: str,
) -> tuple[str, str, str]:
    """Read historical PtdGlc names but expose the requested Glc-GP identity."""
    if str(compound_class).strip() != "PtdGlc":
        return compound_class, lipid_name, lipid_chain_name
    def rename(value: str) -> str:
        return re.sub(r"^PtdGlc(?=\(|\s|$)", "Glc-GP", str(value).strip())
    return "Glc-GP", rename(lipid_name), rename(lipid_chain_name)


def chain_pair(name: str) -> tuple[tuple[int, int], tuple[int, int]] | None:
    match = re.fullmatch(
        r"(?:PtdGlc|Glc-GP|6-Ac-Glc-GP)(?:\s+(\d+):(\d+)[_/](\d+):(\d+)|\((\d+):(\d+)[_/](\d+):(\d+)\))",
        str(name).strip(),
    )
    if match is None:
        return None
    c1, db1, c2, db2 = map(int, (v for v in match.groups() if v is not None))
    if any(c < 2 or db < 0 or db >= c for c, db in ((c1, db1), (c2, db2))):
        return None
    return (c1, db1), (c2, db2)


def negative_glc_gp_fragments(
    name: str, *, precursor_mz: float | None = None,
) -> tuple[float, str, list[FragmentRecord]]:
    chains = chain_pair(name)
    if chains is None:
        raise ValueError(f"Not a diacyl Glc-GP composition: {name}")
    acetylated = str(name).strip().startswith("6-Ac-Glc-GP")
    total_c = sum(c for c, _ in chains)
    total_db = sum(db for _, db in chains)
    nc, nh = total_c + 9, 2 * total_c + 15 - 2 * total_db
    no = 13
    if acetylated:
        nc, nh, no = nc + 2, nh + 2, no + 1
    formula = f"C{nc}H{nh}O{no}P"
    exact_precursor = nc * C + nh * H + no * O + P - PROTON
    precursor = exact_precursor if precursor_mz is None else float(precursor_mz)
    residue = "C8H12O6" if acetylated else "C6H10O5"
    hg = precursor - (ACETYLGLYCOSYL_RESIDUE if acetylated else GLYCOSYL_RESIDUE)
    fragments = [
        FragmentRecord(hg, f"[M-{residue}-H]-", "Diagnostic_HG"),
        FragmentRecord(precursor, "[M-H]-", "Precursor Ion"),
    ]
    fragments.append(FragmentRecord(GLYCEROPHOSPHATE_ANION, "[C3H6O5P]-", "Diagnostic_HG"))
    # A symmetric ion has two chain positions but only one physical RCOO peak.
    for c, db in sorted(set(chains)):
        token = f"{c}:{db}"
        fatty_acid = c * C + (2 * c - 2 * db) * H + 2 * O
        fragments.extend([
            FragmentRecord(fatty_acid - PROTON, f"[RCOO]-({token})", "Diagnostic_FA"),
            FragmentRecord(hg - fatty_acid, f"[M-(ROOH)-{residue}-H]-({token})", "Common"),
        ])
        if not acetylated:
            fragments.append(FragmentRecord(precursor - (fatty_acid - WATER),
                                             f"[M-(R=O)-H]-({token})", "Common"))
    if any(fragment.mz <= 0 for fragment in fragments):
        raise ValueError(f"Nonpositive Glc-GP product mass: {name}")
    return precursor, formula, sorted(fragments, key=lambda fragment: fragment.mz)


# Compatibility for earlier local curation callers; exported names are Glc-GP.
negative_ptdglc_fragments = negative_glc_gp_fragments
