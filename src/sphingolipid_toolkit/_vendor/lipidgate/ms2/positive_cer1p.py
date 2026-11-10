"""Positive Cer1P phosphate-loss evidence and sodium total compositions."""
from __future__ import annotations

from dataclasses import dataclass

from .esterified_ceramide import C, H, N, O, WATER
from .models import FragmentRecord
from .positive_hexosylceramide import normalize_fragments
from .positive_pc_sodium import P, NA, SODIUM_ION

PHOSPHORIC_ACID = 3*H + P + 4*O
SODIUM_PHOSPHATE = 2*H + NA + P + 4*O
GATE_NAMES = {"[M+H]+": "M+H-H3PO4", "[M+Na]+": "M+Na-H3PO4"}
POOL_WEIGHTS = {"fah": 0.0, "hg": 60.0, "other": 40.0}


def is_positive_cer1p(record):
    return record.compound_class == "Cer1P" and record.adduct in GATE_NAMES


def normalize_h_fragments(chain_name, precursor_mz, fragments):
    result = [f for f in fragments if f.name in {GATE_NAMES["[M+H]+"], "M+H-H2O"}
              or f.name == "LCB" or f.name.startswith("LCB-")]
    by_name = {f.name: f for f in result}
    for name, loss, role in [(GATE_NAMES["[M+H]+"], PHOSPHORIC_ACID, "Diagnostic_HG"),
                             ("M+H-H2O", WATER, "Common")]:
        if name not in by_name:
            result.append(FragmentRecord(round(precursor_mz-loss,4), name, role))
    return normalize_fragments("Cer1P", chain_name, precursor_mz, result)


@dataclass(frozen=True, order=True)
class Cer1PSodiumComposition:
    carbon: int
    double_bonds: int
    hydroxyls: int

    def __post_init__(self):
        if self.carbon < 2 or not 0 <= self.double_bonds < self.carbon or self.hydroxyls not in {1,2,3}:
            raise ValueError("Outside the existing m/d/t Cer1P total-composition grid")

    @property
    def name(self):
        prefix={1:'m',2:'d',3:'t'}[self.hydroxyls]
        return f"Cer1P({prefix}{self.carbon}:{self.double_bonds})"

    @property
    def formula(self):
        return f"C{self.carbon}H{2*self.carbon+2-2*self.double_bonds}NO{self.hydroxyls+4}P"

    @property
    def precursor_mz(self):
        return (self.carbon*C + (2*self.carbon+2-2*self.double_bonds)*H
                + N + (self.hydroxyls+4)*O + P + SODIUM_ION)

    def fragments(self):
        p=self.precursor_mz
        return sorted([
            FragmentRecord(p-PHOSPHORIC_ACID, GATE_NAMES["[M+Na]+"], "Diagnostic_HG"),
            FragmentRecord(p-WATER, "M+Na-H2O", "Common"),
            FragmentRecord(p-PHOSPHORIC_ACID-WATER, "M+Na-H3PO4-H2O", "Common"),
            FragmentRecord(p-SODIUM_PHOSPHATE, "M+Na-H2NaPO4", "Common"),
            FragmentRecord(p, "[M+Na]+", "Precursor Ion"),
        ], key=lambda f:f.mz)
