"""Theoretical protonated GD2 spectra on the d-series ganglioside chain grid."""
from __future__ import annotations

from dataclasses import dataclass

from .esterified_ceramide import PROTON, WATER
from .models import FragmentRecord
from .positive_gm2 import HEXNAC
from .positive_gm3 import B_FRAGMENT_TYPE, HEX, LCB_NAME, NEU5AC_RESIDUE, PositiveGm3

DISIALYL = 2 * NEU5AC_RESIDUE
GLYCAN = DISIALYL + HEXNAC + 2 * HEX
B_IONS = {
    "B4": GLYCAN + PROTON,
    "B3": GLYCAN - HEX + PROTON,
    "B2alpha": DISIALYL + PROTON,
    "B4-H2O": GLYCAN + PROTON - WATER,
    "B1alpha": NEU5AC_RESIDUE + PROTON,
    "B3Y2beta": DISIALYL + HEX + PROTON,
    "B4Y1beta": DISIALYL + 2 * HEX + PROTON,
    "B4Y2alpha": HEXNAC + 2 * HEX + PROTON,
    "B3Y2alphaY2beta": NEU5AC_RESIDUE + HEX + PROTON,
    "B3Y2alpha": HEXNAC + HEX + PROTON,
}
YZ_LOSSES = {
    "Y2beta": HEXNAC,
    "Y2alphaY2beta": DISIALYL + HEXNAC,
    "Y2alpha": DISIALYL,
    "Z2beta": HEXNAC + WATER,
    "Y1": DISIALYL + HEXNAC + HEX,
    "Z0": GLYCAN + WATER,
}


@dataclass(frozen=True)
class PositiveGd2:
    gm3: PositiveGm3

    @property
    def name(self) -> str:
        return "GD2" + self.gm3.name.removeprefix("GM3")

    @property
    def species_name(self) -> str:
        return "GD2" + self.gm3.species_name.removeprefix("GM3")

    @property
    def formula(self) -> str:
        c, h, n, o = self.gm3.formula_counts
        return f"C{c + 19}H{h + 30}N{n + 2}O{o + 13}"

    @property
    def precursor_mz(self) -> float:
        return self.gm3.precursor_mz + HEXNAC + NEU5AC_RESIDUE

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        return sorted([
            FragmentRecord(mh, "[M+H]+", "Precursor Ion"),
            FragmentRecord(mh - WATER, "M+H-H2O", "Common"),
            FragmentRecord(self.gm3.lcb_2h2o_mz, LCB_NAME, "Common"),
            FragmentRecord(self.gm3.lcb_2h2o_mz + WATER, "LCB-H2O", "Common"),
            FragmentRecord(NEU5AC_RESIDUE + PROTON - WATER, "B2alphaZ3alpha", "Common"),
            FragmentRecord(mh - GLYCAN, "Y0", "Common"),
            FragmentRecord(B_IONS["B3"] - WATER, "B3-H2O", "Common"),
            *(FragmentRecord(mh - loss, name, "Diagnostic_HG") for name, loss in YZ_LOSSES.items()),
            *(FragmentRecord(mz, name, B_FRAGMENT_TYPE) for name, mz in B_IONS.items()),
        ], key=lambda fragment: fragment.mz)
