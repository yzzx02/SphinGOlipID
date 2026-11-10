"""Nonhydroxy d-series PE-Cer [M+H]+; masses derive from composition."""
from __future__ import annotations

from dataclasses import dataclass
import re

from .esterified_ceramide import C, H, N, O, PROTON, ELECTRON, WATER
from .models import FragmentRecord

P = 30.97376199842
SODIUM = 22.9897692820 - ELECTRON
PHOSPHOETHANOLAMINE = 2 * C + 8 * H + N + 4 * O + P


def is_positive_pe_cer(compound_class: str, adduct: str) -> bool:
    return compound_class == "PE-Cer" and adduct == "[M+H]+"


@dataclass(frozen=True)
class PositivePECer:
    lcb_c: int
    lcb_db: int
    fa_c: int
    fa_db: int

    @property
    def name(self) -> str:
        return f"PE-Cer(d{self.lcb_c}:{self.lcb_db}/{self.fa_c}:{self.fa_db})"

    @property
    def species_name(self) -> str:
        return f"PE-Cer(d{self.lcb_c + self.fa_c}:{self.lcb_db + self.fa_db})"

    @property
    def formula_counts(self) -> tuple[int, int]:
        c = self.lcb_c + self.fa_c + 2
        return c, 2 * c + 3 - 2 * (self.lcb_db + self.fa_db)

    @property
    def formula(self) -> str:
        c, h = self.formula_counts
        return f"C{c}H{h}N2O6P"

    @property
    def precursor_mz(self) -> float:
        c, h = self.formula_counts
        return c * C + h * H + 2 * N + 6 * O + P + PROTON

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        lcb = self.lcb_c * C + (2 * self.lcb_c + 3 - 2 * self.lcb_db) * H + N + 2 * O + PROTON
        # Acyl-derived amide ion: C(n+2)H(2n+2-2db)NO+, e.g. C18H34NO+.
        amide = (self.fa_c + 2) * C + (2 * self.fa_c + 2 - 2 * self.fa_db) * H + N + O - ELECTRON
        return sorted([
            FragmentRecord(mh, "[M+H]+", "Common"),
            FragmentRecord(mh - PHOSPHOETHANOLAMINE, "M+H-141", "Diagnostic_HG"),
            FragmentRecord(lcb - 2 * WATER, "LCB-2H2O", "LCB碎片"),
            FragmentRecord(mh - PHOSPHOETHANOLAMINE - WATER, "M+H-141-H2O", "Common"),
            FragmentRecord(mh - PROTON - PHOSPHOETHANOLAMINE + SODIUM, "M-141+Na", "Common"),
            FragmentRecord(amide, "FA-amide", "Common"),
        ], key=lambda f: f.mz)


def from_name(name: str, lipid_class: str = "PE-Cer") -> PositivePECer | None:
    match = re.fullmatch(rf"{re.escape(lipid_class)}\(d(\d+):(\d+)/(\d+):(\d+)\)", name)
    return PositivePECer(*(int(x) for x in match.groups())) if match else None
