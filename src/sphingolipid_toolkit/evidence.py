"""Archived evidence data vocabulary; new annotation uses the current engine.

Literal labels and containers remain readable for historical exports. The old
label classifier and origin-capture adapters are deliberately non-executable.
"""
from __future__ import annotations
from dataclasses import asdict,dataclass
from enum import Enum
from .backend_policy import removed_backend

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
    DIAGNOSTIC_NL = "diagnostic_nl"
    SUPPORTING_NL = "supporting_nl"
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


EVIDENCE_PRIORITY = {EvidenceType.HG.value: 0, EvidenceType.LCB.value: 1,
    EvidenceType.NL.value: 2, EvidenceType.COMMON.value: 3, EvidenceType.PRECURSOR.value: 4}


class ScoringPool(str, Enum):
    PRIMARY = "primary"
    SECONDARY = "secondary"
    SUPPORT = "support"


POOL_PRIORITY = {"primary":0, "secondary":1, "support":2}


classify_fragment=removed_backend
capture_fragment_origins=removed_backend
record_lcb_origin=removed_backend
fragment_source=removed_backend
scoring_pool_for_evidence=removed_backend
annotate_theoretical_fragments=removed_backend
