"""Hsu 2016 [M-H]- product-ion sets for selected non-esterified Cer families."""
from __future__ import annotations

from dataclasses import replace
import re
from typing import Iterable

from .models import FragmentRecord

C = 12.0
H = 1.00782503223
O = 15.99491461957
WATER = 2 * H + O
FORMALDEHYDE = C + 2 * H + O
CARBON_MONOXIDE = C + O
H2_PLUS_FORMALDEHYDE = 4 * H + C + O
C3H4O = 3 * C + 4 * H + O

NEGATIVE_ADDUCTS = frozenset({"[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"})
CER_RE = re.compile(
    r"^Cer\((?P<series>[dt])(?P<lcb_c>\d+):(?P<lcb_db>\d+)/"
    r"(?P<hydroxy>h?)(?P<fa_c>\d+):(?P<fa_db>\d+)\)$"
)


def is_simple_negative_cer(name: str, adduct: str) -> bool:
    return adduct in NEGATIVE_ADDUCTS and CER_RE.fullmatch(str(name or '').strip()) is not None


def selected_family(name: str, adduct: str) -> str | None:
    """Return one of the three families that require an Hsu-2016 correction.

    d:nFA with unsaturated LCB (for example d18:1/24:1) already agrees with
    the paper and is deliberately not rewritten. A bare ``h`` is treated as
    alpha-hydroxy for this pass; beta/omega-specific names are out of scope.
    """
    if adduct not in NEGATIVE_ADDUCTS:
        return None
    match = CER_RE.fullmatch(str(name or "").strip())
    if match is None:
        return None
    series = match['series']
    lcb_db = int(match['lcb_db'])
    hydroxy = bool(match['hydroxy'])
    if series == 'd' and lcb_db == 0 and not hydroxy:
        return 'd0_nfa'
    if series == 'd' and lcb_db == 0 and hydroxy:
        return 'd0_alpha_hfa'
    if series == 't' and lcb_db == 0 and not hydroxy:
        return 't0_nfa'
    return None


def logical_negative_cer_fragment_types(fragment: FragmentRecord) -> set[str]:
    """Expose both interpretations of one exactly isobaric physical ion."""
    types = {fragment.fragment_type}
    for alias in fragment.name.split(' | '):
        if alias.startswith('LCB-') or alias in {
            'M-H-H2O-RCONH', 'M-H-HCHO-H2-RCONH',
        }:
            types.add('LCB碎片')
        elif alias.startswith(('[RCOO', '[RCONH')):
            types.add('FA类碎片')
        elif alias.startswith('[NAE'):
            types.add('NAE碎片')
    return types


def merge_negative_cer_isobars(
    name: str, adduct: str, fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Represent every exact-m/z product ion once, retaining all assignments."""
    source = list(fragments)
    if not is_simple_negative_cer(name, adduct):
        return source
    grouped: dict[float, list[FragmentRecord]] = {}
    for fragment in source:
        grouped.setdefault(round(float(fragment.mz), 4), []).append(fragment)
    result = []
    for mz, physical_group in grouped.items():
        if len(physical_group) == 1:
            result.append(physical_group[0])
            continue
        aliases = []
        for fragment in physical_group:
            for alias in fragment.name.split(' | '):
                if alias not in aliases:
                    aliases.append(alias)
        if len(aliases) != len(physical_group):
            raise ValueError(f'{name} has duplicate assignments at {mz:.4f}')
        primary = physical_group[0]
        result.append(replace(
            primary,
            name=' | '.join(aliases),
            intensity=max(fragment.intensity for fragment in physical_group),
            weight=max(fragment.weight for fragment in physical_group),
            required_group=None,
        ))
    return sorted(result, key=lambda fragment: (fragment.mz, fragment.name))


def normalize_negative_cer_hsu2016_fragments(
    name: str,
    precursor_mz: float,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    family = selected_family(name, adduct)
    source = list(fragments)
    if family is None:
        return merge_negative_cer_isobars(name, adduct, source)

    by_name: dict[str, FragmentRecord] = {}
    for fragment in source:
        if fragment.name in by_name:
            raise ValueError(f"Duplicate Cer fragment name: {name} {fragment.name}")
        by_name[fragment.name] = fragment

    def require(fragment_name: str) -> FragmentRecord:
        fragment = by_name.get(fragment_name)
        if fragment is None:
            fragment = next(
                (
                    candidate for candidate in by_name.values()
                    if fragment_name in candidate.name.split(' | ')
                ),
                None,
            )
        if fragment is None:
            raise ValueError(f"{name} lacks required Hsu-2016 source ion: {fragment_name}")
        return fragment

    def add_from(fragment: FragmentRecord, shift: float, new_name: str, role: str) -> None:
        expected_mz = round(float(fragment.mz) + shift, 4)
        aliased = [
            other for other in by_name.values()
            if new_name in other.name.split(' | ')
        ]
        if aliased:
            if len(aliased) != 1 or abs(float(aliased[0].mz) - expected_mz) > 0.00011:
                raise ValueError(f"{name} has conflicting alias {new_name}")
            return
        existing = by_name.get(new_name)
        if existing is not None:
            if abs(float(existing.mz) - expected_mz) > 0.00011 or existing.fragment_type != role:
                raise ValueError(f"{name} has conflicting {new_name}")
            return
        collisions = [
            (old_name, other) for old_name, other in by_name.items()
            if abs(float(other.mz) - expected_mz) <= 0.00005
        ]
        if collisions:
            if len(collisions) != 1:
                raise ValueError(f"{name} has multiple collisions at {expected_mz:.4f}")
            old_name, other = collisions[0]
            # Store the physical peak once. Search-time logical typing lets it
            # support either interpretation without duplicating a peak match.
            by_name.pop(old_name)
            aliased_name = f'{other.name} | {new_name}'
            by_name[aliased_name] = replace(other, name=aliased_name)
            return
        by_name[new_name] = replace(
            fragment,
            mz=expected_mz,
            name=new_name,
            fragment_type=role,
            intensity=100.0,
            weight=1.0,
            required_group=None,
        )

    if family == 'd0_nfa':
        # Hsu explicitly reports [M-H-HCHO]- as absent for d18:0/nFA and
        # Table 1 does not report the analogous b2/b4 routes in this family.
        for fragment_name in ('M-H-HCHO', 'LCB-H-HCHO', 'M-H-H2O-RCONH'):
            by_name.pop(fragment_name, None)

    elif family == 'd0_alpha_hfa':
        by_name.pop('LCB-H-HCHO', None)
        mh = require('[M-H]-')
        lcb = require('LCB-H')
        add_from(mh, -(2 * WATER + FORMALDEHYDE), 'M-H-2H2O-HCHO', 'C类碎片')
        add_from(lcb, CARBON_MONOXIDE, 'LCB-H+CO', 'LCB碎片')

    elif family == 't0_nfa':
        match = CER_RE.fullmatch(name)
        assert match is not None
        fa = f"{int(match['fa_c'])}:{int(match['fa_db'])}"
        rconh = require(f'[RCONH]-({fa})')
        lcb = require('LCB-H')
        add_from(rconh, C3H4O, f'[RCONH+C3H4O]-({fa})', 'FA类碎片')
        add_from(lcb, -H2_PLUS_FORMALDEHYDE, 'LCB-H-H2-HCHO', 'LCB碎片')

    return merge_negative_cer_isobars(name, adduct, by_name.values())
