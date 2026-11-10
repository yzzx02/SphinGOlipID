"""Theoretical nonhydroxy d-series GM3 [M+H]+ spectra."""
from __future__ import annotations

from dataclasses import dataclass, replace
import re

from .esterified_ceramide import C, H, N, O, PROTON, WATER
from .models import FragmentRecord

HEX = 6 * C + 10 * H + 5 * O
NEU5AC_RESIDUE = 11 * C + 17 * H + N + 8 * O
HG_LOSSES = {
    "M+H-C11H17NO8": NEU5AC_RESIDUE,
    "M+H-C11H17NO8-H2O": NEU5AC_RESIDUE + WATER,
    "M+H-C11H17NO8-1Hex": NEU5AC_RESIDUE + HEX,
    "M+H-C11H17NO8-1Hex-H2O": NEU5AC_RESIDUE + HEX + WATER,
    "M+H-C11H17NO8-2Hex": NEU5AC_RESIDUE + 2 * HEX,
    "M+H-C11H17NO8-2Hex-H2O": NEU5AC_RESIDUE + 2 * HEX + WATER,
}
B1_NAME = "[Neu5Ac+H-H2O]+"
B_FRAGMENT_TYPE = "Diagnostic_B"
LCB_NAME = "LCB-2H2O"
B_IONS = {
    B1_NAME: NEU5AC_RESIDUE + PROTON,
    "[Neu5Ac+Hex+H-H2O]+": NEU5AC_RESIDUE + HEX + PROTON,
    "[Neu5Ac+2Hex+H-H2O]+": NEU5AC_RESIDUE + 2 * HEX + PROTON,
    "[Neu5Ac+H-2H2O]+": NEU5AC_RESIDUE + PROTON - WATER,
}
GM3_HG = frozenset(HG_LOSSES)
_LEGACY_SPECTRUM_NAMES = frozenset({
    "[M+H]+", B1_NAME, "[Neu5Ac+H-2H2O]+", "LCB-H2O", LCB_NAME,
    "M+H-Neu5Ac", "M+H-Neu5Ac-Hex", "M+H-Neu5Ac-2Hex", "M+H-Neu5Ac-2Hex-H2O",
})


@dataclass(frozen=True)
class PositiveGm3:
    lcb_c: int
    lcb_db: int
    fa_c: int
    fa_db: int

    @property
    def name(self) -> str:
        return f'GM3(d{self.lcb_c}:{self.lcb_db}/{self.fa_c}:{self.fa_db})'

    @property
    def species_name(self) -> str:
        return f'GM3(d{self.lcb_c+self.fa_c}:{self.lcb_db+self.fa_db})'

    @property
    def formula_counts(self) -> tuple[int, int, int, int]:
        carbon = self.lcb_c + self.fa_c
        db = self.lcb_db + self.fa_db
        return carbon + 23, 2 * carbon + 38 - 2 * db, 2, 21

    @property
    def formula(self) -> str:
        c, h, n, o = self.formula_counts
        return f'C{c}H{h}N{n}O{o}'

    @property
    def precursor_mz(self) -> float:
        c, h, n, o = self.formula_counts
        return c * C + h * H + n * N + o * O + PROTON

    @property
    def in_bundled_chain_grid(self) -> bool:
        return 14 <= self.lcb_c <= 24 and 0 <= self.lcb_db <= 2 and 0 <= self.fa_db <= 2

    @property
    def lcb_2h2o_mz(self) -> float:
        lcb = self.lcb_c * C + (2 * self.lcb_c + 3 - 2 * self.lcb_db) * H + N + 2 * O + PROTON
        return lcb - 2 * WATER

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        return sorted([
            FragmentRecord(mh, '[M+H]+', 'Precursor Ion'),
            FragmentRecord(mh - WATER, 'M+H-H2O', 'Common'),
            FragmentRecord(self.lcb_2h2o_mz, LCB_NAME, 'Common'),
            FragmentRecord(self.lcb_2h2o_mz + WATER, 'LCB-H2O', 'Common'),
            *(FragmentRecord(mh - loss, name, 'Diagnostic_HG') for name, loss in HG_LOSSES.items()),
            *(FragmentRecord(mz, name, B_FRAGMENT_TYPE) for name, mz in B_IONS.items()),
        ], key=lambda f:f.mz)


def from_d_hexcer_name(name: str) -> PositiveGm3 | None:
    match = re.fullmatch(r'HexCer\(d(\d+):(\d+)/(\d+):(\d+)\)', name)
    if match is None:
        return None
    return PositiveGm3(*(int(part) for part in match.groups()))


def from_gm3_name(name: str) -> PositiveGm3:
    match = re.fullmatch(r'GM3\(d(\d+):(\d+)/(\d+):(\d+)\)', name)
    if match is None:
        raise ValueError(f'Unsupported positive GM3 (requires nonhydroxy d-series): {name}')
    return PositiveGm3(*(int(part) for part in match.groups()))


def upgrade_legacy_spectrum(
    name: str, precursor_mz: float, fragments: list[FragmentRecord],
) -> list[FragmentRecord] | None:
    """Upgrade the earlier generated nine-ion template on MSP import."""
    if len(fragments) != len(_LEGACY_SPECTRUM_NAMES) or {f.name for f in fragments} != _LEGACY_SPECTRUM_NAMES:
        return None
    try:
        model = from_gm3_name(name)
    except ValueError:
        return None
    if abs(model.precursor_mz - precursor_mz) > 0.0001:
        return None
    current = model.fragments()
    masses = {f.name: f.mz for f in current}
    expected = {
        '[M+H]+': model.precursor_mz, B1_NAME: masses[B1_NAME],
        '[Neu5Ac+H-2H2O]+': masses['[Neu5Ac+H-2H2O]+'],
        'LCB-H2O': masses[LCB_NAME] + WATER, LCB_NAME: masses[LCB_NAME],
        'M+H-Neu5Ac': masses['M+H-C11H17NO8-H2O'],
        'M+H-Neu5Ac-Hex': masses['M+H-C11H17NO8-1Hex-H2O'],
        'M+H-Neu5Ac-2Hex': masses['M+H-C11H17NO8-2Hex-H2O'],
        'M+H-Neu5Ac-2Hex-H2O': masses['M+H-C11H17NO8-2Hex-H2O'] - WATER,
    }
    if any(f.intensity != 100 or f.weight != 1 or abs(f.mz - expected[f.name]) > 0.0001 for f in fragments):
        return None
    return [replace(f, mz=round(f.mz, 4)) for f in current]
