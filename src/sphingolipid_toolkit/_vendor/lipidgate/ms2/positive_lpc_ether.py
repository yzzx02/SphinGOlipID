"""Protonated ether LPC: one184 HG gate, optional support and O/P reports."""
from __future__ import annotations

from dataclasses import replace
import re

from .models import FragmentRecord, LibraryRecord

HG_NAME = "[C5H15NO4P]+"


def composition(record: LibraryRecord) -> tuple[int, int] | None:
    """Return one O-equivalent composition, preserving the source identity."""
    cls = str(record.compound_class).strip().upper()
    if cls not in {"LPC-O", "LPC-P"}:
        return None
    link = "P" if cls == "LPC-P" else "O"
    name = str(record.lipid_chain_name or record.lipid_name or "").strip()
    match = re.fullmatch(
        rf"(?:LPC(?:-{link})?|PC)(?:\(|\s+)(?:{link}-)?(\d+):(\d+)(?:[/_]0:0)?\)?",
        name, flags=re.IGNORECASE,
    )
    if match is None:
        return None
    carbon, db = map(int, match.groups())
    equivalent = db + int(link == "P")
    return (carbon, equivalent) if carbon > 0 and 0 <= db < carbon and equivalent < carbon else None


def is_positive_ether_lpc(record: LibraryRecord) -> bool:
    return str(record.adduct).strip() == "[M+H]+" and composition(record) is not None


def positive_ether_lpc_fragments(precursor_mz, fragments=()):
    source = list(fragments)

    def ion(mz, name, role="Common", group=None):
        old = next((f for f in source if abs(float(f.mz)-mz) <= .002), None)
        return (replace(old, name=name, fragment_type=role, required_group=group) if old is not None
                else FragmentRecord(round(mz, 4), name, role, required_group=group))

    return sorted([
        ion(86.0964, "[C5H12N]+"),
        ion(104.1070, "[C5H14NO]+"),
        ion(124.9998, "[C2H6O4P]+"),
        ion(184.0733, HG_NAME, "Diagnostic_HG", "hg"),
        replace(ion(float(precursor_mz), "[M+H]+", "Precursor Ion"), mz=float(precursor_mz)),
    ], key=lambda f: f.mz)


def normalize_positive_ether_lpc(record: LibraryRecord) -> LibraryRecord:
    if not is_positive_ether_lpc(record):
        return record
    fragments = positive_ether_lpc_fragments(record.precursor_mz, record.fragments)
    return record if fragments == record.fragments else replace(record, fragments=fragments)


def reported_name(record: LibraryRecord) -> str | None:
    if not is_positive_ether_lpc(record):
        return None
    carbon, db = composition(record)
    name = f"LPC(O-{carbon}:{db})"
    return name if db == 0 else f"{name}/LPC(P-{carbon}:{db-1})"
