"""Negative ether LPE fragment definitions and unresolved O/P report names.

The source identity stays a single O composition. Reporting an isomer pair does
not create a second chain or claim a vinyl-ether bond or an sn position.
"""
from __future__ import annotations

from dataclasses import replace
import re

from .models import FragmentRecord, LibraryRecord

ETHANOLAMINE = 61.05276384961  # C2H7NO, neutral monoisotopic mass.
HG140 = 140.0118
HG_NAMES = frozenset({"[C2H7NO4P]-", "[M-H-C2H7NO]-"})


def composition(name: object) -> tuple[int, int] | None:
    matched = re.fullmatch(
        r"(?:LPE(?:-O)?|PE)(?:\(|\s+)(?:O-)?(\d+):(\d+)(?:[/_]0:0)?\)?",
        str(name or "").strip(), flags=re.IGNORECASE,
    )
    if matched is None:
        return None
    carbon, double_bonds = map(int, matched.groups())
    return (carbon, double_bonds) if carbon > 0 and 0 <= double_bonds < carbon else None


def is_negative_lpe(record: LibraryRecord) -> bool:
    return (str(record.compound_class).strip().upper() == "LPE-O"
            and str(record.adduct).strip() == "[M-H]-"
            and composition(record.lipid_chain_name or record.lipid_name) is not None)


def negative_lpe_fragments(name: str, precursor_mz: float, fragments=()) -> list[FragmentRecord]:
    chain = composition(name)
    if chain is None:
        raise ValueError(f"Not a single ether LPE composition: {name}")
    carbon, double_bonds = chain
    source = list(fragments)

    def ion(mz, annotation, role, group=None):
        # Keep an existing source mass/intensity when it represents this ion.
        # In particular, do not move a four-decimal source ion by re-rounding
        # a loss from its independently rounded precursor.
        existing = next((f for f in source if abs(float(f.mz)-mz) <= .002), None)
        if existing is not None:
            return replace(existing, name=annotation, fragment_type=role, required_group=group)
        return FragmentRecord(round(mz, 4), annotation, role, required_group=group)

    ro = 12.0*carbon + 1.00782503223*(2*carbon+1-2*double_bonds) + 15.99491462 + .00054858
    result = [
        ion(78.9591, "[PO3]-", "Common"),
        ion(96.9696, "[H2PO4]-", "Common"),
        ion(HG140, "[C2H7NO4P]-", "Diagnostic_HG", "hg"),
        ion(196.0380, "[C5H11NO5P]-", "Common"),
        ion(ro, f"[R-O]-(O-{carbon}:{double_bonds})", "Common"),
        ion(float(precursor_mz)-ETHANOLAMINE, "[M-H-C2H7NO]-", "Diagnostic_HG", "hg"),
        replace(ion(float(precursor_mz), "[M-H]-", "Precursor Ion", "precursor"), mz=float(precursor_mz)),
    ]
    return sorted(result, key=lambda f: f.mz)


def normalize_negative_lpe_record(record: LibraryRecord) -> LibraryRecord:
    if not is_negative_lpe(record):
        return record
    fragments = negative_lpe_fragments(record.lipid_chain_name or record.lipid_name,
                                       record.precursor_mz, record.fragments)
    return record if fragments == record.fragments else replace(record, fragments=fragments)


def reported_name(record: LibraryRecord) -> str | None:
    if not is_negative_lpe(record):
        return None
    carbon, double_bonds = composition(record.lipid_chain_name or record.lipid_name)
    name = f"LPE(O-{carbon}:{double_bonds})"
    return name if double_bonds == 0 else f"{name}/LPE(P-{carbon}:{double_bonds-1})"
