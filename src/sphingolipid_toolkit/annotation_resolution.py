"""Figure 4 hierarchy, deliberately disconnected from production pipelines."""

from dataclasses import dataclass, replace
from enum import Enum
from typing import Callable

from .annotation_evidence import EvidenceRole, EvidenceSummary, MatchedEvidence


FIGURE4_WORKFLOW = (
    "MS1 Candidates", "MS/MS Targets", "Species Level",
    "Molecular Species", "RT Validation",
)


class AnnotationLevel(str, Enum):
    SPECIES = "Species level"
    MOLECULAR_SPECIES = "Molecular species level"


class EvaluationStatus(str, Enum):
    SATISFIED = "SATISFIED"
    NOT_SATISFIED = "NOT_SATISFIED"
    UNSPECIFIED = "UNSPECIFIED"


@dataclass(frozen=True)
class EvidenceEvaluation:
    status: EvaluationStatus
    reason: str


# A rule sees only matched evidence of the applicable role. RT is never input.
EvidenceRule = Callable[[tuple[MatchedEvidence, ...]], EvidenceEvaluation]


@dataclass(frozen=True)
class SubclassRules:
    species_rule: EvidenceRule | None = None
    molecular_species_rule: EvidenceRule | None = None
    provenance: str = "TODO/UNSPECIFIED: subclass requirements are not defined."


class RuleRegistry:
    """Exact subclass lookup; no universal fallback or inferred thresholds."""

    def __init__(self) -> None:
        self._rules: dict[str, SubclassRules] = {}

    def register(self, subclass: str, rules: SubclassRules) -> None:
        if not subclass:
            raise ValueError("An explicit subclass is required")
        if subclass in self._rules:
            raise ValueError(f"Rules already registered for {subclass}")
        self._rules[subclass] = rules

    def get(self, subclass: str) -> SubclassRules:
        return self._rules.get(subclass, SubclassRules())


def _evaluate(rule: EvidenceRule | None, evidence: tuple[MatchedEvidence, ...]) -> EvidenceEvaluation:
    if rule is None:
        return EvidenceEvaluation(EvaluationStatus.UNSPECIFIED, "TODO/UNSPECIFIED: subclass rule required.")
    if not evidence:
        return EvidenceEvaluation(EvaluationStatus.NOT_SATISFIED, "No eligible matched MS/MS evidence.")
    result = rule(evidence)
    if not isinstance(result, EvidenceEvaluation) or not isinstance(result.status, EvaluationStatus):
        raise TypeError("Rules must return EvidenceEvaluation with EvaluationStatus")
    return result


def evaluate_species_evidence(summary: EvidenceSummary, registry: RuleRegistry) -> EvidenceEvaluation:
    return _evaluate(registry.get(summary.subclass).species_rule,
                     summary.matched_for(EvidenceRole.SPECIES_EVIDENCE))


def _evaluate_molecular(summary: EvidenceSummary, registry: RuleRegistry,
                        species: EvidenceEvaluation) -> EvidenceEvaluation:
    if species.status != EvaluationStatus.SATISFIED:
        return EvidenceEvaluation(species.status, "Species level prerequisite has not been established.")
    evidence = tuple(
        e for e in summary.matched_for(EvidenceRole.MOLECULAR_SPECIES_EVIDENCE)
        if e.metadata.chain_composition
    )
    return _evaluate(registry.get(summary.subclass).molecular_species_rule, evidence)


def evaluate_molecular_species_evidence(summary: EvidenceSummary, registry: RuleRegistry) -> EvidenceEvaluation:
    """Re-evaluate the Species prerequisite for this same candidate/spectrum."""
    return _evaluate_molecular(summary, registry, evaluate_species_evidence(summary, registry))


@dataclass(frozen=True)
class RTValidation:
    """Externally supplied orthogonal validation; performs no RT computation."""

    status: str = "UNSPECIFIED"
    note: str = ""


@dataclass(frozen=True)
class AnnotationResolution:
    summary: EvidenceSummary
    level: AnnotationLevel | None
    species: EvidenceEvaluation
    molecular_species: EvidenceEvaluation
    rt_validation: RTValidation | None = None


def resolve_annotation(summary: EvidenceSummary, registry: RuleRegistry) -> AnnotationResolution:
    species = evaluate_species_evidence(summary, registry)
    molecular = _evaluate_molecular(summary, registry, species)
    level = None
    if species.status == EvaluationStatus.SATISFIED:
        level = AnnotationLevel.SPECIES
        if molecular.status == EvaluationStatus.SATISFIED:
            level = AnnotationLevel.MOLECULAR_SPECIES
    return AnnotationResolution(summary, level, species, molecular)


def attach_rt_validation(resolution: AnnotationResolution, validation: RTValidation) -> AnnotationResolution:
    """Attach a result without promoting, demoting, or re-evaluating MS/MS level."""
    return replace(resolution, rt_validation=validation)
