"""Auditable fragment metadata; no masses or fragmentation rules are generated here."""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict, dataclass
from enum import Enum
from functools import wraps
import re

import pandas as pd


class EvidenceType(str, Enum):
    PRECURSOR = "Precursor"
    HG = "HG"
    LCB = "LCB"
    NL = "NL"
    COMMON = "common"


class EvidenceRole(str, Enum):
    PRECURSOR = "precursor"
    CLASS_DIAGNOSTIC = "class_diagnostic"
    CHAIN_SPECIFIC = "chain_specific"
    STRUCTURE_INFORMATIVE = "structure_informative"
    SUPPORTING = "supporting"


@dataclass(frozen=True)
class EvidenceRecord:
    fragment_type: EvidenceType
    evidence_role: EvidenceRole
    scoring_eligible: bool
    diagnostic_strength: str
    gate_eligible: bool = False
    classification_status: str = "classified"
    classification_source: str = "legacy_label"
    eligibility_reason: str = "existing discrete rule"
    lcb_identity: str = ""

    def to_dict(self):
        return {k: v.value if isinstance(v, Enum) else v for k, v in asdict(self).items()}


@dataclass(frozen=True)
class FragmentOrigin:
    function: str
    representative: bool | None = None
    lcb_identity: str = ""


# Literal labels from the existing functions. No mass-only classification:
# an isobar without the correct origin/class is not a diagnostic ion.
HG_LABELS = {
    "SM": ("[phosphocholine+H]+",),
    "PE_cer": ("[phosphoethanolamine+H]+",),
    "PI_cer": ("[phosphate+H]+",),
    "PG_cer": ("[phosphoglycerol+H]+",),
    "CAEP": ("[C2H8NO3P+H]+",),
    "N_CAEP": ("[C3H10NO3P+H]+]",),  # legacy spelling retained
    "CerP": ("H4PO4",), "S1P": ("H4PO4",),
    "Lyso_SM": ("[phosphocholine+H]+",),
    "Lyso_sulfo": ("[HSO4+H]+",),
    "GSL_fragments": ("NeuAc+H", "NeuAc+H-H2O"),
}
SUPPORT_LABELS = frozenset({"[C2H5NO+H]+", "[phosphocholine-H3PO4+H]+",
    "P-O断裂", "[C3H8O3+H]+"})
LOSS_TOKENS = {
    "SM": ("C5H14NO4P", "C3H9N"), "PE_cer": ("phosphoethanolamine", "C2H5N"),
    "PI_cer": ("phosphate", "C6H10O5"), "PG_cer": ("phosphoglycerol", "C3H8O3"),
    "CAEP": ("C2H8NO3P",), "N_CAEP": ("C3H10NO3P",),
    "CerP": ("H3PO4",), "S1P": ("H3PO4",),
    "Lyso_SM": ("C5H14NO4P", "C3H9N"), "Lyso_sulfo": ("C6H10O8S", "SO3"),
    "Glu_So": ("C6H10O5", "C6H10O5H2O"),
    "Gb3_So": ("C6H10O5", "2*C6H10O5", "3*C6H10O5"),
    "FMC_1": ("Gal", "HOAc"), "FMC_3": ("Gal_OAc", "HOAc"),
    "FMC_5": ("Gal_4OAc", "HOAc"), "EO_cer": ("E 18:2",),
    "EO_Glccer": ("E 18:2", "Glc"),
}
SOURCE_FUNCTIONS = tuple(dict.fromkeys((*HG_LABELS, *LOSS_TOKENS, "So", "O_FA_cer", "OCer")))


def _record(kind, *, eligible=True, diagnostic=False, source="legacy_label", reason="existing discrete rule", lcb="", status="classified"):
    role = {EvidenceType.PRECURSOR: EvidenceRole.PRECURSOR, EvidenceType.HG: EvidenceRole.CLASS_DIAGNOSTIC,
        EvidenceType.LCB: EvidenceRole.CHAIN_SPECIFIC, EvidenceType.NL: EvidenceRole.STRUCTURE_INFORMATIVE,
        EvidenceType.COMMON: EvidenceRole.SUPPORTING}[kind]
    strength = "diagnostic" if diagnostic else "informative" if kind == EvidenceType.NL else "supporting"
    return EvidenceRecord(kind, role, eligible, strength,
        diagnostic and eligible and kind in {EvidenceType.HG, EvidenceType.LCB, EvidenceType.NL},
        status, source, reason, lcb)


def classify_fragment(lipid_class, fragment_name, fragment_origin=None, *, registry=None) -> EvidenceRecord:
    """Prefer captured generator provenance; exact legacy labels are a fallback.

    Unknown labels remain common/UNSPECIFIED, excluded from gates and score.
    Bare dehydration is common; a glycan/headgroup loss plus water stays NL,
    with the water satellite excluded from the default scoring denominator.
    """
    from .evidence_config import get_class_rule
    from . import ms2_legacy_core as core

    label = str(fragment_name)
    origin = fragment_origin
    if isinstance(origin, str):
        origin = FragmentOrigin(origin)
    rule = get_class_rule(lipid_class, registry)
    sources = (origin.function,) if origin else rule.generator_functions
    provenance = "generator:" + origin.function if origin else "legacy_label"
    if origin and origin.function == "MS1":
        return _record(EvidenceType.PRECURSOR, eligible=False, source=provenance, reason="MS1 prerequisite only")
    # Require the label to exist in the unchanged LCB dictionary, not a regex
    # claiming that an arbitrary chain-like label is experimental evidence.
    lcb_key = re.match(r"^[mdtq]\d+:\d+", label)
    if lcb_key and label in core.LCB_fragment.get(lcb_key.group(), ()):
        return _record(EvidenceType.LCB, diagnostic=True, source=provenance, lcb=lcb_key.group())
    if any(label in HG_LABELS.get(source, ()) for source in sources):
        return _record(EvidenceType.HG, diagnostic=True, source=provenance)
    if rule.lipid_class == "Cer" and re.fullmatch(r"(?:Cer\([^)]*\)|M\+H|\[M\+H\]\+?)-(?:1?H2O|2H2O)", label):
        return _record(EvidenceType.NL, diagnostic=True, source="author:Cer_dehydration",
            reason="Author designated Cer mono/didehydration as diagnostic NL; specificity requires validation")
    # Remove the intact lipid annotation before checking a loss suffix.
    suffix = label.rsplit(")", 1)[-1] if ")" in label else label[label.find("-"):]
    suffix = suffix.split(" [", 1)[0].strip()
    loss = suffix.lstrip("-")
    water = "H2O" in loss
    tokens = [token.strip() for token in loss.split("-")]
    if label in SUPPORT_LABELS or (tokens and all(re.fullmatch(r"\d*H2O", token) for token in tokens)):
        return _record(EvidenceType.COMMON, source=provenance)
    if "GSL_fragments" in sources:
        from .glycan_encoding import RESIDUES
        residues = [token for token in tokens if token and not re.fullmatch(r"\d*H2O", token)]
        if residues and all(token in RESIDUES for token in residues):
            representative = origin is not None and origin.representative is True and not water
            return _record(EvidenceType.NL, eligible=representative, diagnostic=True, source=provenance,
                reason="first/complete unhydrated glycan loss" if representative else "glycan explanation; representative set unspecified or satellite")
    for source in sources:
        allowed = LOSS_TOKENS.get(source, ())
        meaningful = [token for token in tokens if token and not re.fullmatch(r"\d*H2O", token)]
        if meaningful and all(token in allowed for token in meaningful):
            # Fixed class-specific losses are informative; small amine loss
            # alone is not a subclass gate. No gate is ever based on water.
            diagnostic = any(token not in {"C3H9N", "C2H5N", "C3H8O3", "HOAc"} for token in meaningful)
            return _record(EvidenceType.NL, eligible=not water, diagnostic=diagnostic,
                source=provenance, reason="unhydrated class-specific loss" if not water else "water satellite; explanation only")
    if "OCer" in sources and re.fullmatch(r"\d+:\d+", loss):
        return _record(EvidenceType.NL, eligible=False, source=provenance, status="UNSPECIFIED",
            reason="OCer fixed acyl-loss enumeration; assignment/representatives TODO")
    return _record(EvidenceType.COMMON, eligible=False, source=provenance, status="UNSPECIFIED",
        reason="unverified legacy label; no diagnostic inference")


# Capture source only while an explicit adapter requests metadata. Decorated
# legacy calls otherwise return immediately through their original function.
_capture = ContextVar("fragment_evidence_capture", default=None)


@contextmanager
def capture_fragment_origins():
    records = {}
    token = _capture.set(records)
    try:
        yield records
    finally:
        _capture.reset(token)


def record_lcb_origin(identity, labels, masses):
    records = _capture.get()
    if records is not None:
        for label, mass in zip(labels, masses):
            records[(label, float(mass))] = FragmentOrigin("LCB_fragment", lcb_identity=identity)


def fragment_source(function):
    """Observe existing output without replacing a fragment label or mass."""
    @wraps(function)
    def wrapped(*args, **kwargs):
        records = _capture.get()
        if records is None:
            return function(*args, **kwargs)
        namespace = function.__globals__
        start = len(namespace["frag_id1"])
        result = function(*args, **kwargs)
        labels = namespace["frag_id1"][start:]
        masses = namespace["frag_mz1"][start:]
        representatives = set()
        if function.__name__ == "GSL_fragments":
            # Select two topological anchors BEFORE seeing any observed peak.
            # These are a provisional transparent eligibility policy, not a
            # claim of experimentally calibrated observability.
            from .glycan_encoding import parse_glycan_encoding
            import inspect
            bound = inspect.signature(function).bind(*args, **kwargs).arguments
            parsed = parse_glycan_encoding(bound["input_str"], bound["location"])
            name = bound["name"]
            first = list(parsed.main_chain[:1]) + [r for b in parsed.branches if b.position == 0 for r in b.residues]
            complete = list(parsed.main_chain) + [r for b in parsed.branches for r in b.residues]
            # Loss order in the legacy encoding can differ; multiset equality
            # identifies the same anchor without calculating a new mass.
            for label in labels:
                if not label.startswith(name + "-") or "H2O" in label:
                    continue
                losses = label[len(name)+1:].split(" [", 1)[0].split("-")
                if sorted(losses) in (sorted(first), sorted(complete)):
                    representatives.add(label)
        for label, mass in zip(labels, masses):
            records[(label, float(mass))] = FragmentOrigin(function.__name__, label in representatives)
        return result
    return wrapped


EVIDENCE_PRIORITY = {EvidenceType.HG.value: 0, EvidenceType.LCB.value: 1,
    EvidenceType.NL.value: 2, EvidenceType.COMMON.value: 3, EvidenceType.PRECURSOR.value: 4}


def annotate_theoretical_fragments(table, lipid_class, origins=None, eligibility_overrides=None, *, registry=None):
    """Add metadata and canonicalize exact-mass evidence; preserve legacy columns.

    All labels are retained. Equal masses count once with HG > LCB > NL >
    common, irrespective of which label happens to be first. Explicit eligibility
    overrides are keyed by exact fragment label and must be documented by caller.
    """
    data = table.copy()
    name_col = next((c for c in ("fragment_name", "idf", "id") if c in data), None)
    mass_col = next((c for c in ("theoretical_mz", "mz") if c in data), None)
    if name_col is None or mass_col is None:
        raise KeyError("Theory needs fragment_name/idf/id and theoretical_mz/mz")
    rows = []
    # CSV float parsing can alter the last binary bit. Provenance may fall back
    # to an exact label only when every captured use has identical metadata.
    by_label = {}
    for (label, _mass), origin in (origins or {}).items():
        by_label.setdefault(label, set()).add(origin)
    for _, row in data.iterrows():
        mass = float(row[mass_col])
        for label in str(row[name_col]).split(" | "):
            origin = (origins or {}).get((label, mass))
            if origin is None and len(by_label.get(label, ())) == 1:
                origin = next(iter(by_label[label]))
            if origin is None and "fragment_origin" in row and pd.notna(row["fragment_origin"]):
                origin = row["fragment_origin"]
            metadata = classify_fragment(lipid_class, label, origin, registry=registry).to_dict()
            if eligibility_overrides and label in eligibility_overrides:
                metadata["scoring_eligible"] = bool(eligibility_overrides[label])
                metadata["gate_eligible"] = metadata["scoring_eligible"] and metadata["diagnostic_strength"] == "diagnostic"
                metadata["eligibility_reason"] = "explicit caller eligibility override"
            rows.append({**row.to_dict(), **metadata, "fragment_name":label,"theoretical_mz":mass})
    if not rows:
        return pd.DataFrame(columns=list(data.columns) + ["fragment_name", "theoretical_mz"] + list(_record(EvidenceType.COMMON).to_dict()))
    expanded = pd.DataFrame(rows)
    result = []
    groups = ["anno", "theoretical_mz"] if "anno" in expanded else ["theoretical_mz"]
    for _, group in expanded.groupby(groups, sort=True, dropna=False):
        ordered = group.assign(_priority=group.fragment_type.map(EVIDENCE_PRIORITY)).sort_values(
            ["_priority", "scoring_eligible", "fragment_name"], ascending=[True,False,True], kind="stable")
        chosen = ordered.iloc[0].drop(labels="_priority").to_dict()
        chosen["selected_evidence_label"] = chosen["fragment_name"]
        chosen["fragment_name"] = " | ".join(sorted(set(group.fragment_name)))
        chosen["all_evidence_types"] = " | ".join(sorted(set(group.fragment_type), key=EVIDENCE_PRIORITY.get))
        # Preserve the original label field, merging instead of losing labels.
        chosen[name_col] = chosen["fragment_name"]
        result.append(chosen)
    return pd.DataFrame(result)
