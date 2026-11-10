"""Archived result containers; experimental scoring has been retired.

New identification and complete provenance come from the current engine. RT
metadata can still be attached to an archived result without changing its score.
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from math import isfinite
import numpy as np
import pandas as pd
from .backend_policy import removed_backend

class AnnotationLevel(str, Enum):
    CANDIDATE = "candidate"
    SPECIES = "species"
    MOLECULAR_SPECIES = "molecular_species"

    @property
    def display_name(self):
        return {self.CANDIDATE:"MS1 candidate", self.SPECIES:"Species level",
            self.MOLECULAR_SPECIES:"Molecular species level"}[self]


@dataclass(frozen=True)
class PrecursorEvidence:
    precursor_pass: bool
    precursor_error_ppm: float | None
    precursor_adduct: str | None
    source: str = "existing_ms1_match"
    observed_mz: float | None = None

    def __post_init__(self):
        if not isinstance(self.precursor_pass, (bool, np.bool_)):
            raise ValueError("precursor_pass must be an explicit boolean")
        if self.observed_mz is not None and (not isfinite(self.observed_mz) or self.observed_mz <= 0):
            raise ValueError("Observed precursor m/z must be finite and positive")
        if self.precursor_pass and (self.precursor_error_ppm is None or
            not isfinite(self.precursor_error_ppm) or not self.precursor_adduct):
            raise ValueError("Passing precursor evidence requires actual error ppm and adduct")

    @classmethod
    def from_existing_match(cls, passed, observed_mz, theoretical_mz, adduct):
        """Record the existing MS1 decision; introduce no new mass tolerance."""
        if not isfinite(observed_mz) or not isfinite(theoretical_mz) or theoretical_mz <= 0:
            raise ValueError("Finite precursor masses and positive theoretical mass required")
        return cls(passed, (observed_mz-theoretical_mz)/theoretical_mz*1e6, adduct, observed_mz=observed_mz)


@dataclass
class CandidateEvidenceResult:
    summary: dict
    fragments: pd.DataFrame


def add_rt_validation(result: CandidateEvidenceResult, status, *, error=None):
    """Attach orthogonal RT validation; never modify gates, score or level."""
    return CandidateEvidenceResult({**result.summary,"rt_validation_status":status,"rt_error":error},result.fragments.copy())


evaluate_annotation_gates=removed_backend
structural_quality=removed_backend
support_count_score=removed_backend
evaluate_candidate_evidence=removed_backend
rank_evidence_candidates=removed_backend
