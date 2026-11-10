"""Curated negative HexCer HG/FAH/LCB evidence with physical-peak aliases."""
from __future__ import annotations

from dataclasses import replace
import re
from typing import Iterable

from .models import FragmentRecord

NEGATIVE_ADDUCTS = {"[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"}
LCB_TYPE = "LCB碎片"
WATER = 18.01056468403


def is_negative_hexcer(compound_class: str, adduct: str) -> bool:
    return compound_class == "HexCer" and adduct in NEGATIVE_ADDUCTS


def logical_fragment_types(fragment: FragmentRecord) -> set[str]:
    types = {fragment.fragment_type}
    # An LCB P/R ion can be exactly isobaric with RCOO- or V. Store/match
    # that physical peak once, while retaining both explicit interpretations.
    if " | LCB fragment " in fragment.name:
        types.add(LCB_TYPE)
    return types


def normalize_negative_hexcer_fragments(
    name: str, precursor_mz: float, adduct: str, fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    if adduct not in NEGATIVE_ADDUCTS:
        return list(fragments)
    acyl_match = re.fullmatch(r"HexCer\([mdt]\d+:\d+/(?P<fa>h?\d+:\d+)\)", name)
    if acyl_match is None:
        raise ValueError(f"Unsupported negative HexCer name: {name}")
    acyl = acyl_match['fa']
    source = list(fragments)

    def select(names: tuple[str, ...], new_name: str, role: str) -> FragmentRecord:
        candidates = [f for f in source if any(n in f.name.split(' | ') for n in names)]
        if not candidates:
            raise ValueError(f"{name} lacks negative HexCer evidence: {names}")
        original = candidates[0]
        return replace(original, name=new_name, fragment_type=role,
                       required_group='hg' if role == 'Diagnostic_HG' else 'fah' if role == 'Diagnostic_FA' else None)

    rcoo = next((f for f in source if f.name.startswith('[RCOO]-(')), None)
    if rcoo is None:
        raise ValueError(f"{name} lacks RCOO evidence")
    rcoo_name = rcoo.name.split(' | ')[0]
    v_source = next((f for f in source if f.name.split(' | ')[0] == f'V ion ({acyl})'), None)
    v_ion = replace(v_source or rcoo, mz=round(rcoo.mz - WATER, 4),
                    name=f'V ion ({acyl})', fragment_type='Diagnostic_FA', required_group='fah')
    selected = [
        select(('M-H-C6H10O5',), 'M-H-C6H10O5', 'Diagnostic_HG'),
        replace(rcoo, name=rcoo_name, fragment_type='Diagnostic_FA', required_group='fah'),
        # The old library assigned these two letter labels in reverse.
        select(('Ceramide fragment T', f'S ion ({acyl})'), f'S ion ({acyl})', 'Diagnostic_FA'),
        select(('Ceramide fragment S', f'T ion ({acyl})'), f'T ion ({acyl})', 'Diagnostic_FA'),
        v_ion,
        select(('Ceramide fragment P', 'LCB fragment P'), 'LCB fragment P', LCB_TYPE),
        select(('Ceramide fragment R', 'LCB fragment R'), 'LCB fragment R', LCB_TYPE),
    ]
    if adduct != '[M-H]-':
        # Adduct spectra have a second HG interpretation. The adduct precursor
        # and the direct [M-H]- precursor are not part of the scoring spectrum.
        selected.append(select(('[M-H]-',), '[M-H]-', 'Diagnostic_HG'))

    physical: dict[float, FragmentRecord] = {}
    for fragment in selected:
        key = round(fragment.mz, 4)
        previous = physical.get(key)
        if previous is None:
            physical[key] = fragment
        elif previous.fragment_type == 'Diagnostic_FA' and fragment.fragment_type == LCB_TYPE:
            physical[key] = replace(previous, name=previous.name + ' | ' + fragment.name,
                                    intensity=max(previous.intensity, fragment.intensity))
        else:
            raise ValueError(f"Unexpected HexCer physical peak collision: {name} {key}")
    return sorted(physical.values(), key=lambda f: f.mz)
