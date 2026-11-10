"""Sodiated ordinary Cer: total composition and precursor-derived losses only."""
from __future__ import annotations

from dataclasses import dataclass

from .esterified_ceramide import C, H, N, O, WATER, FORMALDEHYDE, is_esterified_ceramide
from .models import FragmentRecord
from .positive_pc_sodium import SODIUM_ION

GATE_NAME = "M+Na-H2O"
OTHER_NAMES = frozenset({"M+Na-2H2O", "M+Na-CH2O-H2O", "[M+Na]+"})
POOL_WEIGHTS = {"fah": 0.0, "hg": 60.0, "other": 40.0}


def is_positive_cer_sodium(record) -> bool:
    return (record.compound_class == "Cer" and record.adduct == "[M+Na]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name))


@dataclass(frozen=True, order=True)
class CerSodiumComposition:
    carbon: int
    double_bonds: int
    hydroxyls: int

    def __post_init__(self):
        if self.carbon < 2 or not 0 <= self.double_bonds < self.carbon or self.hydroxyls not in {1, 2, 3, 4}:
            raise ValueError("Invalid ordinary Cer total composition")

    @property
    def name(self):
        # Existing library convention: prefix means summed hydroxylation here,
        # never an assignment of OH to the LCB or the fatty-acyl chain.
        prefix = {1: "m", 2: "d", 3: "t", 4: "t"}[self.hydroxyls]
        return f"Cer({prefix}{self.carbon}:{self.double_bonds})" + ("(OH)" if self.hydroxyls == 4 else "")

    @property
    def formula(self):
        return f"C{self.carbon}H{2*self.carbon+1-2*self.double_bonds}NO{self.hydroxyls+1}"

    @property
    def precursor_mz(self):
        return (self.carbon*C + (2*self.carbon+1-2*self.double_bonds)*H
                + N + (self.hydroxyls+1)*O + SODIUM_ION)

    def fragments(self):
        p = self.precursor_mz
        return sorted([
            FragmentRecord(p-WATER, GATE_NAME, "C类碎片"),
            FragmentRecord(p-2*WATER, "M+Na-2H2O", "Common"),
            FragmentRecord(p-WATER-FORMALDEHYDE, "M+Na-CH2O-H2O", "Common"),
            FragmentRecord(p, "[M+Na]+", "Precursor Ion"),
        ], key=lambda f: f.mz)
