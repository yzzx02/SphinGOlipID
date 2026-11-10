"""Single-chain LPC sodium: mandatory NL59 and four ordinary fragments."""
from __future__ import annotations

from dataclasses import replace
import re

from .esterified_ceramide import C, H, N, O, ELECTRON
from .models import FragmentRecord
from .positive_pc_sodium import P, SODIUM_ION, TMA, SODIUM_PHOSPHOCHOLINE


def canonical_lpc_name(name):
    match = re.fullmatch(r"LPC(?:\s+|\()(\d+):(\d+)(?:[_/](\d+):(\d+))?\)?", str(name).strip())
    if match is None:
        raise ValueError(f"Not an unmodified LPC: {name}")
    chains = [tuple(map(int, match.groups()[:2]))]
    if match.group(3) is not None:
        chains.append(tuple(map(int, match.groups()[2:])))
    nonzero = [chain for chain in chains if chain != (0, 0)]
    if len(nonzero) != 1:
        raise ValueError(f"LPC requires one acyl chain: {name}")
    carbon, double_bonds = nonzero[0]
    if carbon < 2 or not 0 <= double_bonds < carbon:
        raise ValueError(f"Invalid LPC chain: {name}")
    return f"LPC({carbon}:{double_bonds})"


def sodium_lpc_fragments(name, *, precursor_mz=None):
    name = canonical_lpc_name(name)
    carbon, double_bonds = map(int, re.fullmatch(r"LPC\((\d+):(\d+)\)", name).groups())
    nc, nh = carbon + 8, 2 * carbon - 2 * double_bonds + 18
    formula = f"C{nc}H{nh}NO7P"
    mass = nc * C + nh * H + N + 7 * O + P + SODIUM_ION
    precursor = mass if precursor_mz is None else float(precursor_mz)
    fragments = [
        FragmentRecord(precursor - TMA, "[M+Na-C3H9N]+", "Diagnostic_HG", required_group="hg"),
        FragmentRecord(precursor - SODIUM_PHOSPHOCHOLINE, "[M+Na-205.0480]+", "Common"),
        FragmentRecord(5 * C + 14 * H + N + O - ELECTRON, "[C5H14NO]+", "Common"),
        FragmentRecord(5 * C + 12 * H + N - ELECTRON, "[C5H12N]+", "Common"),
        FragmentRecord(precursor, "[M+Na]+", "Precursor Ion"),
    ]
    return precursor, formula, sorted(fragments, key=lambda fragment: fragment.mz)


def normalize_lpc_sodium(record):
    if record.compound_class != "LPC" or record.adduct != "[M+Na]+":
        return record
    try:
        name = canonical_lpc_name(record.lipid_chain_name)
    except ValueError:
        return record
    _, formula, fragments = sodium_lpc_fragments(name, precursor_mz=record.precursor_mz)
    return replace(record, lipid_name=name, lipid_chain_name=name, formula=formula, fragments=fragments)
