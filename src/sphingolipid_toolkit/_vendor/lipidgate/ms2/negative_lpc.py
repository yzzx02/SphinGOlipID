"""Negative ether LPC: required demethylated ion, ordinary support and O/P names."""
from __future__ import annotations

from dataclasses import replace
import re

from .models import FragmentRecord, LibraryRecord

C, H, N, O, P = 12.0, 1.00782503223, 14.00307400443, 15.99491461957, 30.97376199842
ELECTRON = .000548579909
NL89 = 4*C + 11*H + N + O  # Net C4H11NO loss from [M-CH3]-, not from the adduct.
HG_NAME = "[M-CH3]-"
ADDUCTS = frozenset({"[M+HCOO]-", "[M+CH3COO]-", "[M+C2H3O2]-"})


def composition(name: object) -> tuple[int, int] | None:
    match = re.fullmatch(
        r"(?:LPC(?:-O)?|PC)(?:\(|\s+)(?:O-)?(\d+):(\d+)(?:[/_]0:0)?\)?",
        str(name or "").strip(), flags=re.IGNORECASE,
    )
    if match is None:
        return None
    carbon, db = map(int, match.groups())
    return (carbon, db) if carbon > 0 and 0 <= db < carbon else None


def is_negative_lpc(record: LibraryRecord) -> bool:
    return (str(record.compound_class).strip().upper() == "LPC-O"
            and str(record.adduct).strip() in ADDUCTS
            and composition(record.lipid_chain_name or record.lipid_name) is not None)


def negative_lpc_fragments(name, precursor_mz, adduct, fragments=()):
    chain = composition(name)
    if chain is None or adduct not in ADDUCTS:
        raise ValueError(f"Not a negative single ether LPC: {name} {adduct}")
    carbon, db = chain
    neutral = (carbon+8)*C + (2*carbon-2*db+20)*H + N + 6*O + P
    demethylated = neutral - (C+3*H) + ELECTRON
    ro = carbon*C + (2*carbon+1-2*db)*H + O + ELECTRON
    source = list(fragments)

    def ion(mz, annotation, role="Common", group=None):
        old = next((f for f in source if abs(float(f.mz)-mz) <= .002), None)
        if old is not None:
            return replace(old, name=annotation, fragment_type=role, required_group=group)
        return FragmentRecord(round(mz, 4), annotation, role, required_group=group)

    ions = [
        ion(P+3*O+ELECTRON, "[PO3]-"),
        ion(4*C+11*H+N+4*O+P+ELECTRON, "[C4H11NO4P]-"),
        ion(7*C+15*H+N+5*O+P+ELECTRON, "[C7H15NO5P]-"),
        ion(demethylated-NL89, "[M-CH3-C4H11NO]-"),
        ion(ro, f"[R-O]-(O-{carbon}:{db})"),
        # The audited RO-H2 mass is ordinary support, without a mechanism,
        # chain-position or vinyl-ether identity claim.
        ion(ro-2*H, f"[R-O-H2]-(O-{carbon}:{db})"),
        ion(demethylated, HG_NAME, "Diagnostic_HG", "hg"),
        replace(ion(float(precursor_mz), adduct, "Precursor Ion"), mz=float(precursor_mz)),
    ]
    return sorted(ions, key=lambda f: f.mz)


def normalize_negative_lpc_record(record: LibraryRecord) -> LibraryRecord:
    if not is_negative_lpc(record):
        return record
    fragments = negative_lpc_fragments(record.lipid_chain_name or record.lipid_name,
                                       record.precursor_mz, record.adduct, record.fragments)
    return record if fragments == record.fragments else replace(record, fragments=fragments)


def reported_name(record: LibraryRecord) -> str | None:
    if not is_negative_lpc(record):
        return None
    carbon, db = composition(record.lipid_chain_name or record.lipid_name)
    name = f"LPC(O-{carbon}:{db})"
    return name if db == 0 else f"{name}/LPC(P-{carbon}:{db-1})"
