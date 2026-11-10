"""Negative GM3: mandatory Neu5Ac HG, optional HexCer P/R chain evidence."""
from __future__ import annotations

from dataclasses import replace
import re
from typing import Iterable, Sequence

from .models import FragmentMatch, FragmentRecord, LibraryRecord

LCB_NAMES = frozenset({"LCB fragment P", "LCB fragment R"})
LCB_TYPE = "LCB碎片"


def is_negative_gm3(compound_class: str, adduct: str) -> bool:
    return compound_class == "GM3" and adduct == "[M-H]-"


def lcb_fragments_from_hexcer(fragments: Iterable[FragmentRecord]) -> list[FragmentRecord]:
    """Copy P/R masses, including a P ion stored as a HexCer V/FA alias."""
    result: dict[str, FragmentRecord] = {}
    for fragment in fragments:
        for name in LCB_NAMES.intersection(fragment.name.split(" | ")):
            if name in result:
                raise ValueError(f"Duplicate HexCer {name}")
            result[name] = replace(fragment, name=name, fragment_type=LCB_TYPE,
                                   required_group=None)
    if set(result) != LCB_NAMES:
        raise ValueError("Negative GM3 requires both library P/R ions from HexCer")
    return sorted(result.values(), key=lambda fragment: fragment.mz)


def has_negative_gm3_chain_evidence(
    record: LibraryRecord, matches: Sequence[FragmentMatch],
) -> bool:
    """One of the two library P/R ions permits a chain-composition report.

    This is a resolution decision, never an identification/rejection gate.
    Sum-only legacy records must not acquire a fabricated chain identity.
    """
    if not is_negative_gm3(record.compound_class, record.adduct):
        return False
    if re.fullmatch(r"GM3\(d\d+:\d+/\d+:\d+\)", record.lipid_chain_name) is None:
        return False
    library_names = {
        fragment.name for fragment in record.fragments
        if fragment.fragment_type == LCB_TYPE and fragment.name in LCB_NAMES
    }
    if library_names != LCB_NAMES:
        return False
    return any(
        match.fragment.fragment_type == LCB_TYPE and match.fragment.name in LCB_NAMES
        for match in matches
    )
