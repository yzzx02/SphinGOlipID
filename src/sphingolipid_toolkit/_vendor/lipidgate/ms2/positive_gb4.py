"""Theoretical nonhydroxy d-series Gb4 [M+H]+ spectra."""
from __future__ import annotations

from dataclasses import dataclass
import re

from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,
    HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,
    OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,
    WATER_MONOISOTOPIC_MASS as WATER,
)
from .models import FragmentRecord


HEX = 6 * C + 10 * H + 5 * O
HEXNAC = 8 * C + 13 * H + N + 5 * O
B2_NAME = "[HexNAc+Hex+H]+"
B2_MZ = HEXNAC + HEX + PROTON
HG_LOSSES = {
    "M+H-HexNAc": HEXNAC,
    "M+H-HexNAc-H2O": HEXNAC + WATER,
    "M+H-HexNAc-1Hex": HEXNAC + HEX,
    "M+H-HexNAc-1Hex-H2O": HEXNAC + HEX + WATER,
    "M+H-HexNAc-2Hex": HEXNAC + 2 * HEX,
    "M+H-HexNAc-2Hex-H2O": HEXNAC + 2 * HEX + WATER,
    "M+H-HexNAc-3Hex": HEXNAC + 3 * HEX,
    "M+H-HexNAc-3Hex-H2O": HEXNAC + 3 * HEX + WATER,
    "M+H-HexNAc-3Hex-2H2O": HEXNAC + 3 * HEX + 2 * WATER,
}
GB4_HG = frozenset(HG_LOSSES)
_HEX3CER_CHAIN = re.compile(r"Hex3Cer\(d(\d+):(\d+)/(\d+):(\d+)\)")


@dataclass(frozen=True)
class PositiveGb4:
    lcb_c: int
    lcb_db: int
    fa_c: int
    fa_db: int

    @property
    def name(self) -> str:
        return f"Gb4(d{self.lcb_c}:{self.lcb_db}/{self.fa_c}:{self.fa_db})"

    @property
    def species_name(self) -> str:
        return f"Gb4(d{self.lcb_c + self.fa_c}:{self.lcb_db + self.fa_db})"

    @property
    def formula_counts(self) -> tuple[int, int, int, int]:
        carbon = self.lcb_c + self.fa_c
        double_bonds = self.lcb_db + self.fa_db
        return carbon + 26, 2 * carbon + 44 - 2 * double_bonds, 2, 23

    @property
    def formula(self) -> str:
        c, h, n, o = self.formula_counts
        return f"C{c}H{h}N{n}O{o}"

    @property
    def precursor_mz(self) -> float:
        c, h, n, o = self.formula_counts
        return c * C + h * H + n * N + o * O + PROTON

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        lcb = self.lcb_c * C + (2 * self.lcb_c + 3 - 2 * self.lcb_db) * H + N + 2 * O + PROTON
        return sorted([
            FragmentRecord(mh, "[M+H]+", "Precursor Ion"),
            FragmentRecord(mh - WATER, "M+H-H2O", "Common"),
            FragmentRecord(lcb - 2 * WATER, "LCB-2H2O", "Common"),
            FragmentRecord(lcb - WATER, "LCB-H2O", "Common"),
            FragmentRecord(B2_MZ, B2_NAME, "Common"),
            *(FragmentRecord(mh - loss, name, "Diagnostic_HG") for name, loss in HG_LOSSES.items()),
        ], key=lambda fragment: fragment.mz)


def from_hex3cer_name(name: str) -> PositiveGb4 | None:
    match = _HEX3CER_CHAIN.fullmatch(name)
    return PositiveGb4(*(int(part) for part in match.groups())) if match else None
