"""Opt-in sidecar interfaces. Never imports or invokes the legacy core.

Metadata lookup requires provenance, not just an m/z or a substring in a lipid
name. Existing fragment IDs, dictionaries, tables and scores are not modified.
"""

from dataclasses import dataclass
from typing import Mapping, Sequence, TypeVar

from .annotation_evidence import EvidenceSummary, EvidenceType, FragmentMetadata, potential_role
from .annotation_resolution import AnnotationResolution, RuleRegistry, resolve_annotation


@dataclass(frozen=True)
class AnnotationResolutionConfig:
    annotation_resolution_enabled: bool = False


_HG_IDS = {
    "SM": {"[phosphocholine+H]+", "[phosphocholine-H3PO4+H]+"},
    "Lyso_SM": {"[phosphocholine+H]+", "[phosphocholine-H3PO4+H]+"},
    "PE_cer": {"[phosphoethanolamine+H]+"},
    "PG_cer": {"[phosphoglycerol+H]+", "[C3H8O3+H]+"},
    "CAEP": {"[C2H8NO3P+H]+"},
    "N_CAEP": {"[C3H10NO3P+H]+]"},
    "CerP": {"H4PO4"},
    "S1P": {"H4PO4"},
    "Lyso_sulfo": {"[HSO4+H]+"},
}


def lookup_legacy_metadata(fragment_id: str, source: str) -> FragmentMetadata:
    evidence_type = None
    if fragment_id in _HG_IDS.get(source, set()):
        evidence_type = EvidenceType.HG
    elif source == "GSL_fragments" and fragment_id in {"NeuAc+H", "NeuAc+H-H2O"}:
        evidence_type = EvidenceType.GLYCAN
    elif source in {"ms2_pipeline", "targeted_mzml_pipeline"} and fragment_id == "[C2H5NO+H]+":
        evidence_type = EvidenceType.SUPPORT
    return FragmentMetadata(fragment_id, source, evidence_type, potential_role(evidence_type))


def lcb_metadata(composition: str, fragment_ids: Mapping[str, Sequence[str]],
                 fragment_mz: Mapping[str, Sequence[float]]) -> tuple[FragmentMetadata, ...]:
    """Read caller-supplied LCB_fragment/LCB_mz; never generate or match peaks.

    Each entry is potential molecular-species evidence only. Isobaric chain
    ambiguity and diagnostic sufficiency remain TODO for subclass rules.
    """
    ids, masses = fragment_ids[composition], fragment_mz[composition]
    if len(ids) != len(masses):
        raise ValueError("LCB fragment IDs and masses must be aligned")
    return tuple(FragmentMetadata(
        fragment_id, "LCB_fragment/LCB_mz", EvidenceType.LCB,
        potential_role(EvidenceType.LCB), chain_composition=composition,
        theoretical_mz=mz,
    ) for fragment_id, mz in zip(ids, masses))


def neutral_loss_metadata(fragment_id: str, source: str,
                          loss_origin: EvidenceType | None = None) -> FragmentMetadata:
    """Explicit provenance hook; all legacy losses remain supporting by default."""
    return FragmentMetadata(fragment_id, source, EvidenceType.NEUTRAL_LOSS,
                            loss_origin=loss_origin)


T = TypeVar("T")


def annotate_sidecar(legacy_result: T, summary: EvidenceSummary, registry: RuleRegistry,
                     config: AnnotationResolutionConfig = AnnotationResolutionConfig(),
                     ) -> tuple[T, AnnotationResolution | None]:
    """Return the identical legacy object and an optional independent result.

    No columns, filtering, sorting, scoring, copying, I/O or RT processing. This
    function is not called by either existing pipeline, even when enabled here.
    """
    if not config.annotation_resolution_enabled:
        return legacy_result, None
    return legacy_result, resolve_annotation(summary, registry)
