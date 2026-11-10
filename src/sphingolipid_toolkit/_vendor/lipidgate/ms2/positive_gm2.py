"""Theoretical nonhydroxy d-series GM2 [M+H]+ spectra."""
from __future__ import annotations

from dataclasses import dataclass

from .esterified_ceramide import C, H, N, O, PROTON, WATER
from .models import FragmentRecord
from .positive_gm3 import B_FRAGMENT_TYPE, HEX, LCB_NAME, NEU5AC_RESIDUE, PositiveGm3

HEXNAC = 8 * C + 13 * H + N + 5 * O
HG_LOSSES = {
    "Y0": NEU5AC_RESIDUE + HEXNAC + 2 * HEX,
    "Z0": NEU5AC_RESIDUE + HEXNAC + 2 * HEX + WATER,
    "Y1": NEU5AC_RESIDUE + HEXNAC + HEX,
    "Z1": NEU5AC_RESIDUE + HEXNAC + HEX + WATER,
    "Y2alphaY2beta": NEU5AC_RESIDUE + HEXNAC,
    "Z2alphaZ2beta": NEU5AC_RESIDUE + HEXNAC + WATER,
    "Y2alpha": NEU5AC_RESIDUE,
    "Z2alpha": NEU5AC_RESIDUE + WATER,
    "Y2beta": HEXNAC,
    "Z2beta": HEXNAC + WATER,
}
B_IONS = {
    "B1alpha": NEU5AC_RESIDUE + PROTON,
    "B2": NEU5AC_RESIDUE + HEXNAC + HEX + PROTON,
    "B3": NEU5AC_RESIDUE + HEXNAC + 2 * HEX + PROTON,
    "B2Y2alpha": HEXNAC + HEX + PROTON,
    "B2Y2beta": NEU5AC_RESIDUE + HEX + PROTON,
    "B3Y2alpha": HEXNAC + 2 * HEX + PROTON,
    "B3Y2beta": NEU5AC_RESIDUE + 2 * HEX + PROTON,
}


@dataclass(frozen=True)
class PositiveGm2:
    gm3: PositiveGm3

    @property
    def name(self) -> str:
        return "GM2" + self.gm3.name.removeprefix("GM3")

    @property
    def species_name(self) -> str:
        return "GM2" + self.gm3.species_name.removeprefix("GM3")

    @property
    def formula(self) -> str:
        c, h, n, o = self.gm3.formula_counts
        return f"C{c + 8}H{h + 13}N{n + 1}O{o + 5}"

    @property
    def precursor_mz(self) -> float:
        return self.gm3.precursor_mz + HEXNAC

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        return sorted([
            FragmentRecord(mh, "[M+H]+", "Precursor Ion"),
            FragmentRecord(mh - WATER, "M+H-H2O", "Common"),
            FragmentRecord(self.gm3.lcb_2h2o_mz, LCB_NAME, "Common"),
            FragmentRecord(self.gm3.lcb_2h2o_mz + WATER, "LCB-H2O", "Common"),
            *(FragmentRecord(mh - loss, name, "Diagnostic_HG") for name, loss in HG_LOSSES.items()),
            *(FragmentRecord(mz, name, B_FRAGMENT_TYPE) for name, mz in B_IONS.items()),
        ], key=lambda fragment: fragment.mz)
