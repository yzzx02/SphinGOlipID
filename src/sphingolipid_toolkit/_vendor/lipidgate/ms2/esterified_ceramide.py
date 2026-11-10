"""Explicit esterified-Cer compositions and diagnostic ions."""
from __future__ import annotations

from dataclasses import dataclass
import re

from .models import FragmentRecord

C = 12.0
H = 1.00782503223
N = 14.00307400443
O = 15.99491461957
PROTON = 1.007276466621
ELECTRON = H - PROTON
WATER = 2 * H + O
FORMALDEHYDE = C + 2 * H + O
def is_esterified_ceramide(compound_class: str, name: str = "") -> bool:
    return compound_class == "Cer" and re.fullmatch(
        r"Cer\s*d\d+:\d+/\d+:\d+\(O-\d+:\d+\)", name,
    ) is not None


@dataclass(frozen=True)
class EsterifiedCeramide:
    lipid_class: str
    lcb_c: int
    lcb_db: int
    acyl_c: int
    acyl_db: int
    outer_c: int
    outer_db: int

    @property
    def name(self) -> str:
        return f"Cerd{self.lcb_c}:{self.lcb_db}/{self.acyl_c}:{self.acyl_db}(O-{self.outer_c}:{self.outer_db})"

    @property
    def formula_counts(self) -> tuple[int, int, int, int]:
        carbons = self.lcb_c + self.acyl_c + self.outer_c
        db = self.lcb_db + self.acyl_db + self.outer_db
        return carbons, 2 * carbons + 3 - 2 * db - 4, 1, 5

    @property
    def formula(self) -> str:
        c, h, n, o = self.formula_counts
        return f"C{c}H{h}NO{o}"

    @property
    def neutral_mass(self) -> float:
        c, h, n, o = self.formula_counts
        return c * C + h * H + n * N + o * O

    def precursor(self, adduct: str) -> float:
        shifts = {
            "[M+H]+": PROTON,
            "[M-H]-": -PROTON,
            "[M+HCOO]-": C + 2 * H + 2 * O - PROTON,
            "[M+CH3COO]-": 2 * C + 4 * H + 2 * O - PROTON,
        }
        return self.neutral_mass + shifts[adduct]

    def fragments(self, adduct: str) -> list[FragmentRecord]:
        fa = self.outer_c * C + (2 * self.outer_c - 2 * self.outer_db) * H + 2 * O
        outer = f"{self.outer_c}:{self.outer_db}"
        if adduct == "[M+H]+":
            mh = self.precursor(adduct)
            lcb = self.lcb_c * C + (2 * self.lcb_c + 3 - 2 * self.lcb_db) * H + N + 2 * O + PROTON
            result = [
                FragmentRecord(mh, "[M+H]+", "Precursor Ion"),
                FragmentRecord(mh - WATER, "M+H-H2O", "Common"),
                FragmentRecord(mh - WATER - fa, f"M+H-H2O-RCOOH({outer})", "Diagnostic_FA_Loss"),
                FragmentRecord(mh - 2 * WATER - fa, f"M+H-2H2O-RCOOH({outer})", "Diagnostic_FA_Loss"),
                FragmentRecord(lcb - WATER, "LCB-H2O", "LCB碎片"),
                FragmentRecord(lcb - 2 * WATER, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(lcb - WATER - FORMALDEHYDE, "LCB-CH2O-H2O", "LCB碎片"),
            ]
        else:
            mh = self.precursor("[M-H]-")
            # The T ion retains two LCB carbons and the omega-hydroxy N-acyl
            # residue: C(n+2)H(2n+2-2DB)NO2-, not the outer ester-linked FA.
            t_ion = (self.acyl_c + 2) * C + (2 * self.acyl_c + 2 - 2 * self.acyl_db) * H + N + 2 * O + ELECTRON
            result = [
                FragmentRecord(fa - PROTON, f"[RCOO]-({outer})", "Diagnostic_FA"),
                FragmentRecord(mh - fa, f"M-H-RCOOH({outer})", "Diagnostic_FA_Loss"),
                FragmentRecord(mh - fa + WATER, f"M-H-ketene({outer})", "Diagnostic_FA_Loss"),
                FragmentRecord(t_ion, f"T ion (omega-hydroxy {self.acyl_c}:{self.acyl_db})", "Diagnostic_FA"),
                FragmentRecord(mh, "[M-H]-", "Precursor Ion" if adduct == "[M-H]-" else "Diagnostic_HG"),
            ]
            if adduct != "[M-H]-":
                result.append(FragmentRecord(self.precursor(adduct), adduct, "Precursor Ion"))
        return sorted(result, key=lambda fragment: fragment.mz)


def parse_esterified_ceramide(name: str) -> EsterifiedCeramide:
    from .sphingolipid_naming import canonicalize_multichain_sphingolipid_name
    name = canonicalize_multichain_sphingolipid_name(name)
    matched = re.fullmatch(
        r"Cerd(?P<lc>\d+):(?P<ld>\d+)/"
        r"(?P<ac>\d+):(?P<ad>\d+)\(O-(?P<oc>\d+):(?P<od>\d+)\)", name,
    )
    if matched is None:
        raise ValueError(f"Unsupported esterified ceramide composition: {name}")
    return EsterifiedCeramide('Cer', *(int(matched[k]) for k in ('lc', 'ld', 'ac', 'ad', 'oc', 'od')))
