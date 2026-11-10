"""Oxidized fragments retain the chemistry and scoring pools of their parent class."""

import re
from dataclasses import replace

from .chain_utils import FRAGMENT_CHAIN_TOKEN_RE, canonical_chain_token
from .esterified_ceramide import ELECTRON, C, H, N, O
from .models import FragmentRecord
from .oxidized_chain_evidence import oxygen_count

WATER = 2 * H + O
PE_HEADGROUP = 2 * C + 8 * H + N + 4 * O + 30.97376199842


def parent_class_key(compound_class):
    key = re.sub(r"[^A-Z0-9]", "", str(compound_class).upper())
    return key.removeprefix("OX")


def is_oxidized_precursor_water_loss(fragment):
    # A chain label would instead describe dehydration of a chain fragment.
    name = re.sub(r"[\s\[\]]", "", fragment.name).upper()
    return re.fullmatch(r"M(?:\+H)?-(?:[123])?H2O(?:[+-])?", name) is not None


def normalize_oxidized_fragments(compound_class, chain_name, precursor, adduct, fragments):
    if not str(compound_class).upper().startswith("OX"):
        return fragments
    result = [
        replace(f, fragment_type="Common", required_group=None)
        if is_oxidized_precursor_water_loss(f) else f
        for f in fragments
    ]
    if adduct != "[M+H]+":
        return result
    tokens = [canonical_chain_token(m) for m in FRAGMENT_CHAIN_TOKEN_RE.finditer(chain_name)]
    cls = parent_class_key(compound_class)
    if cls not in {"PC", "PE"} or len(tokens) != 2 or any(t.startswith(("O-", "P-")) for t in tokens):
        return result
    if not any(oxygen_count(t) for t in tokens):
        return result
    # Replace legacy chain peaks (including duplicate oxygen spelling aliases).
    result = [f for f in result if f.fragment_type not in {"Diagnostic_FA", "Diagnostic_FA_Loss", "FA_Frag"}]
    for token in sorted(set(tokens)):
        carbon, db = map(int, token.split(";")[0].split(":"))
        acid = carbon * C + (2 * carbon - 2 * db) * H + (2 + oxygen_count(token)) * O
        ketene = acid - WATER
        if cls == "PC":
            result.extend([
                FragmentRecord(precursor - acid, f"[M-(ROOH)+H]+({token})", "Diagnostic_FA_Loss"),
                FragmentRecord(precursor - ketene, f"[M-(R=O)+H]+({token})", "Diagnostic_FA_Loss"),
            ])
        else:
            result.extend([
                FragmentRecord(precursor - PE_HEADGROUP - ketene, f"[M-R=O-C2H8O4NP+H]+({token})", "Diagnostic_FA_Loss"),
                FragmentRecord(acid - O - H - ELECTRON, f"(R=O)+({token})", "FA_Frag"),
            ])
    for count in (1, 2):
        mz = precursor - count * WATER
        if not any(abs(f.mz - mz) < 0.002 and is_oxidized_precursor_water_loss(f) for f in result):
            suffix = "H2O" if count == 1 else "2H2O"
            result.append(FragmentRecord(mz, f"[M+H-{suffix}]+", "Common"))
    if cls == "PE":
        mz = precursor - PE_HEADGROUP - WATER
        result = [f for f in result if abs(f.mz - mz) >= 0.002]
        result.append(FragmentRecord(mz, "[M+H-HG-H2O]+", "Common"))
    return sorted(result, key=lambda f: f.mz)
