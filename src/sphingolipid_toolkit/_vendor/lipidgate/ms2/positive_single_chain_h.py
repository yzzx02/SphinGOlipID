"""Explicit protonated MG and single-acyl NAGly definitions.

Only the named [M+H]+ records are affected. Ammoniated MG and legacy
two-chain NAGly retain their existing fragments, gates and score profiles.
"""
from __future__ import annotations

from dataclasses import replace
import re

from .esterified_ceramide import C, H, N, O, ELECTRON, PROTON, WATER
from .models import FragmentRecord, LibraryRecord

GLYCEROL = 3 * C + 8 * H + 3 * O
GLYCINE = 2 * C + 5 * H + N + 2 * O
PROTONATED_GLYCINE = 2 * C + 6 * H + N + 2 * O - ELECTRON
SINGLE_CHAIN_POOL_WEIGHTS = {"fah": 0.0, "hg": 75.0, "other": 25.0}


def single_chain_composition(compound_class: str, name: str) -> tuple[int, int] | None:
    if compound_class not in {"MG", "NAGly"}:
        return None
    m = re.fullmatch(rf"{re.escape(compound_class)}(?:\s+|\()(\d+):(\d+)\)?", str(name).strip())
    if m is None:
        return None
    carbon, double_bonds = map(int, m.groups())
    if carbon < 2 or double_bonds < 0 or double_bonds >= carbon:
        return None
    return carbon, double_bonds


def is_protonated_single_chain_record(record: LibraryRecord) -> bool:
    return str(record.adduct).strip() == '[M+H]+' and single_chain_composition(
        record.compound_class, record.lipid_chain_name,
    ) is not None


def is_protonated_single_chain_nagly(record: LibraryRecord) -> bool:
    return record.compound_class == 'NAGly' and is_protonated_single_chain_record(record)


def protonated_single_chain_fragments(
    compound_class: str, name: str, *, precursor_mz: float | None = None,
) -> tuple[float, str, list[FragmentRecord]]:
    composition = single_chain_composition(compound_class, name)
    if composition is None:
        raise ValueError(f'Not a single-acyl {compound_class}: {name}')
    carbon, double_bonds = composition
    if compound_class == 'MG':
        nc, nh = carbon + 3, 2 * carbon - 2 * double_bonds + 6
        formula = f'C{nc}H{nh}O4'
        mass = nc * C + nh * H + 4 * O + PROTON
    else:
        nc, nh = carbon + 2, 2 * carbon - 2 * double_bonds + 3
        formula = f'C{nc}H{nh}NO3'
        mass = nc * C + nh * H + N + 3 * O + PROTON
    # Existing import metadata is preserved, including its precursor precision.
    precursor = mass if precursor_mz is None else float(precursor_mz)
    if compound_class == 'MG':
        fragments = [
            FragmentRecord(precursor-WATER, '[M-H2O+H]+', 'Diagnostic_HG'),
            FragmentRecord(precursor-GLYCEROL, '[M+H-92]+', 'Common'),
            FragmentRecord(precursor, '[M+H]+', 'Precursor Ion'),
        ]
    else:
        token = f'{carbon}:{double_bonds}'
        fragments = [
            FragmentRecord(PROTONATED_GLYCINE, '[C2H6NO2]+', 'Diagnostic_HG'),
            FragmentRecord(precursor-GLYCINE, f'[M+H-C2H5NO2]+({token})', 'Diagnostic_HG'),
            FragmentRecord(precursor-WATER, '[M+H-H2O]+', 'Common'),
            FragmentRecord(precursor, '[M+H]+', 'Precursor Ion'),
        ]
    if any(f.mz <= 0 for f in fragments):
        raise ValueError(f'Invalid fragment mass: {name}')
    return precursor, formula, sorted(fragments, key=lambda f: f.mz)


def normalize_protonated_single_chain_record(record: LibraryRecord) -> LibraryRecord:
    if not is_protonated_single_chain_record(record):
        return record
    _, _, fragments = protonated_single_chain_fragments(
        record.compound_class, record.lipid_chain_name, precursor_mz=record.precursor_mz,
    )
    return replace(record, fragments=fragments)


def single_chain_rule(record: LibraryRecord, rule):
    """One NAGly HG suffices; no legacy two-chain signature count is required."""
    if not is_protonated_single_chain_nagly(record):
        return rule
    return replace(rule, allow_hg_only_if_no_fah=True,
                   positive_hg_min_matches_if_no_fah=1,
                   positive_hg_min_matches=1, positive_hg_min_fraction=0.0,
                   positive_hg_complete_can_resolve_chain=True,
                   positive_signature_min_matches_if_no_fah=0)
