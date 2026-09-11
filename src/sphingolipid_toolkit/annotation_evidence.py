"""Independent evidence vocabulary; no legacy matching or scoring behavior.

Roles describe potential structural use, never sufficiency. Only explicitly
configured subclass rules may establish an annotation level.
"""

from dataclasses import dataclass
from enum import Enum


class EvidenceType(str, Enum):
    PRECURSOR = "PRECURSOR"
    HG = "HG"
    GLYCAN = "GLYCAN"
    LCB = "LCB"
    FA = "FA"
    NEUTRAL_LOSS = "NEUTRAL_LOSS"
    SUPPORT = "SUPPORT"


class EvidenceRole(str, Enum):
    SPECIES_EVIDENCE = "SPECIES_EVIDENCE"
    MOLECULAR_SPECIES_EVIDENCE = "MOLECULAR_SPECIES_EVIDENCE"
    SUPPORTING_EVIDENCE = "SUPPORTING_EVIDENCE"


def potential_role(evidence_type: EvidenceType | None) -> EvidenceRole:
    """Chemical provenance suggests a role, not a structural decision."""
    if evidence_type in (EvidenceType.HG, EvidenceType.GLYCAN):
        return EvidenceRole.SPECIES_EVIDENCE
    if evidence_type in (EvidenceType.LCB, EvidenceType.FA):
        return EvidenceRole.MOLECULAR_SPECIES_EVIDENCE
    return EvidenceRole.SUPPORTING_EVIDENCE


@dataclass(frozen=True)
class FragmentMetadata:
    fragment_id: str
    source: str
    evidence_type: EvidenceType | None = None  # None means UNSPECIFIED.
    role: EvidenceRole = EvidenceRole.SUPPORTING_EVIDENCE
    chain_composition: str | None = None
    loss_origin: EvidenceType | None = None
    theoretical_mz: float | None = None
    note: str = "TODO/UNSPECIFIED: diagnostic sufficiency requires subclass rules."

    def __post_init__(self) -> None:
        if self.role == EvidenceRole.MOLECULAR_SPECIES_EVIDENCE:
            chain_source = self.evidence_type in (EvidenceType.LCB, EvidenceType.FA)
            chain_loss = (
                self.evidence_type == EvidenceType.NEUTRAL_LOSS
                and self.loss_origin in (EvidenceType.LCB, EvidenceType.FA)
            )
            if not (chain_source or chain_loss):
                raise ValueError("Molecular-species evidence requires a chain-specific source")
        if self.role == EvidenceRole.SPECIES_EVIDENCE and self.evidence_type in (
            None, EvidenceType.PRECURSOR, EvidenceType.SUPPORT,
        ):
            raise ValueError("Unknown, precursor, and supporting fragments cannot establish Species level")


@dataclass(frozen=True)
class MatchedEvidence:
    metadata: FragmentMetadata
    matched: bool = False
    observed_mz: float | None = None
    intensity: float | None = None


@dataclass(frozen=True)
class EvidenceSummary:
    """Evidence for one candidate in one MS/MS context, never pooled by default."""

    candidate_id: str
    subclass: str
    spectrum_id: str
    evidence: tuple[MatchedEvidence, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "evidence", tuple(self.evidence))

    def matched_for(self, role: EvidenceRole) -> tuple[MatchedEvidence, ...]:
        return tuple(e for e in self.evidence if e.matched and e.metadata.role == role)
