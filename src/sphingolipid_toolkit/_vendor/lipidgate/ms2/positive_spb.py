"""Protonated free sphingoid-base LCB roles, without generating new peaks."""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable

from .models import FragmentRecord, LibraryRecord


LCB_NAMES = frozenset({"M+H-H2O", "M+H-2H2O"})
LCB_TYPE = "LCB碎片"


def is_positive_spb(record: LibraryRecord) -> bool:
    return record.compound_class == "SPB" and str(record.adduct).strip() == "[M+H]+"


def fragment_type(fragment: FragmentRecord) -> str:
    if str(fragment.name).strip() in LCB_NAMES:
        return LCB_TYPE
    return "Precursor Ion" if fragment.fragment_type == "Precursor Ion" else "Common"


def normalize_fragments(fragments: Iterable[FragmentRecord]) -> list[FragmentRecord]:
    """Preserve masses, occurrences, order and provenance; only revise roles."""
    return [replace(fragment, fragment_type=fragment_type(fragment), required_group=None)
            for fragment in fragments]
