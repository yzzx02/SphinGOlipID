"""Synthetic rules below are test fixtures, not scientific subclass rules."""

import pytest

from sphingolipid_toolkit.annotation_evidence import (
    EvidenceRole, EvidenceSummary, EvidenceType, FragmentMetadata, MatchedEvidence, potential_role,
)
from sphingolipid_toolkit.annotation_resolution import (
    AnnotationLevel, EvidenceEvaluation, EvaluationStatus, FIGURE4_WORKFLOW,
    RTValidation, RuleRegistry, SubclassRules, attach_rt_validation,
    evaluate_molecular_species_evidence, resolve_annotation,
)


def evidence(kind, matched=True, chain=None, role=None):
    return MatchedEvidence(FragmentMetadata(
        kind.value, "synthetic", kind, role or potential_role(kind), chain_composition=chain,
    ), matched=matched)


def summary(*items, subclass="TEST_ONLY"):
    return EvidenceSummary("synthetic-candidate", subclass, "synthetic-spectrum", items)


def accepts_fixture(items):
    return EvidenceEvaluation(EvaluationStatus.SATISFIED, "Synthetic fixture accepted")


def registry(species=accepts_fixture, molecular=accepts_fixture):
    result = RuleRegistry()
    result.register("TEST_ONLY", SubclassRules(species, molecular, "Synthetic tests only"))
    return result


def test_workflow_and_official_level_names():
    assert FIGURE4_WORKFLOW == (
        "MS1 Candidates", "MS/MS Targets", "Species Level", "Molecular Species", "RT Validation",
    )
    assert {level.value for level in AnnotationLevel} == {"Species level", "Molecular species level"}


def test_no_real_subclass_defaults_or_single_lcb_upgrade():
    data = summary(evidence(EvidenceType.HG), evidence(EvidenceType.LCB, chain="d18:1"))
    result = resolve_annotation(data, RuleRegistry())
    assert result.level is None
    assert result.species.status == EvaluationStatus.UNSPECIFIED
    result = resolve_annotation(data, registry(molecular=None))
    assert result.level == AnnotationLevel.SPECIES
    assert result.molecular_species.status == EvaluationStatus.UNSPECIFIED
    assert resolve_annotation(summary(*data.evidence, subclass="OTHER"), registry()).level is None


def test_molecular_species_requires_species_for_same_summary():
    def forbidden(items):
        pytest.fail("Molecular rule must not run without Species evidence")

    data = summary(evidence(EvidenceType.LCB, chain="d18:1"))
    assert resolve_annotation(data, registry(molecular=forbidden)).level is None
    assert evaluate_molecular_species_evidence(data, registry(molecular=forbidden)).status == EvaluationStatus.NOT_SATISFIED


def test_explicit_synthetic_subclass_rule_controls_chain_sufficiency():
    def requires_fixture_pair(items):
        status = (EvaluationStatus.SATISFIED if {e.metadata.chain_composition for e in items}
                  == {"fixture-A", "fixture-B"} else EvaluationStatus.NOT_SATISFIED)
        return EvidenceEvaluation(status, "Artificial fixture combination, not a production rule")

    rules = registry(molecular=requires_fixture_pair)
    first = summary(evidence(EvidenceType.HG), evidence(EvidenceType.LCB, chain="fixture-A"))
    assert resolve_annotation(first, rules).level == AnnotationLevel.SPECIES
    second = summary(*first.evidence, evidence(EvidenceType.LCB, chain="fixture-B"))
    assert resolve_annotation(second, rules).level == AnnotationLevel.MOLECULAR_SPECIES


@pytest.mark.parametrize("item", [
    evidence(EvidenceType.SUPPORT), evidence(EvidenceType.PRECURSOR),
    evidence(EvidenceType.NEUTRAL_LOSS), evidence(EvidenceType.LCB),
    evidence(EvidenceType.LCB, matched=False, chain="d18:1"),
    evidence(EvidenceType.LCB, chain="d18:1", role=EvidenceRole.SUPPORTING_EVIDENCE),
])
def test_support_unmatched_and_composition_free_evidence_cannot_upgrade(item):
    result = resolve_annotation(summary(evidence(EvidenceType.HG), item), registry())
    assert result.level == AnnotationLevel.SPECIES


def test_unmatched_headgroup_cannot_establish_species():
    result = resolve_annotation(summary(evidence(EvidenceType.HG, matched=False)), registry())
    assert result.level is None


@pytest.mark.parametrize("status", ["pass", "outlier", "UNSPECIFIED"])
@pytest.mark.parametrize("items", [(), (evidence(EvidenceType.HG),),
    (evidence(EvidenceType.HG), evidence(EvidenceType.LCB, chain="d18:1"))])
def test_rt_never_changes_annotation_level(status, items):
    original = resolve_annotation(summary(*items), registry())
    validated = attach_rt_validation(original, RTValidation(status))
    assert validated.level == original.level
    assert validated.species is original.species
    assert validated.molecular_species is original.molecular_species
    assert original.rt_validation is None


def test_registry_rejects_accidental_rule_replacement():
    with pytest.raises(ValueError, match="already registered"):
        registry().register("TEST_ONLY", SubclassRules())
