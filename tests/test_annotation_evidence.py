import pytest

from sphingolipid_toolkit.annotation_evidence import (
    EvidenceRole, EvidenceType, FragmentMetadata, potential_role,
)
from sphingolipid_toolkit.annotation_legacy_adapter import (
    lcb_metadata, lookup_legacy_metadata, neutral_loss_metadata,
)


@pytest.mark.parametrize("kind,role", [
    (EvidenceType.HG, EvidenceRole.SPECIES_EVIDENCE),
    (EvidenceType.GLYCAN, EvidenceRole.SPECIES_EVIDENCE),
    (EvidenceType.LCB, EvidenceRole.MOLECULAR_SPECIES_EVIDENCE),
    (EvidenceType.FA, EvidenceRole.MOLECULAR_SPECIES_EVIDENCE),
    (EvidenceType.PRECURSOR, EvidenceRole.SUPPORTING_EVIDENCE),
    (EvidenceType.SUPPORT, EvidenceRole.SUPPORTING_EVIDENCE),
    (EvidenceType.NEUTRAL_LOSS, EvidenceRole.SUPPORTING_EVIDENCE),
    (None, EvidenceRole.SUPPORTING_EVIDENCE),
])
def test_potential_role_is_separate_from_chemical_type(kind, role):
    assert potential_role(kind) == role


def test_headgroup_and_glycan_metadata_require_source():
    hg = lookup_legacy_metadata("[phosphocholine+H]+", "SM")
    assert hg.evidence_type == EvidenceType.HG
    assert hg.role == EvidenceRole.SPECIES_EVIDENCE
    assert lookup_legacy_metadata(hg.fragment_id, "unknown").evidence_type is None
    assert lookup_legacy_metadata("NeuAc+H", "GSL_fragments").evidence_type == EvidenceType.GLYCAN


@pytest.mark.parametrize("fragment,source", [
    ("P-O断裂", "SM"), ("[phosphate+H]+", "PI_cer"),
    ("mock-18:2", "EO_cer"), ("mock-16:0", "OCer"),
    ("mock", "GSL_fragments"), ("unknown", "unknown"),
])
def test_ambiguous_ids_do_not_invent_fa_or_structural_evidence(fragment, source):
    meta = lookup_legacy_metadata(fragment, source)
    assert meta.evidence_type is None
    assert meta.role == EvidenceRole.SUPPORTING_EVIDENCE
    assert "UNSPECIFIED" in meta.note


def test_lcb_interface_reads_paired_metadata_without_mutation():
    ids, masses = {"d18:1": ["d18:1", "d18:1+OH"]}, {"d18:1": [264.2686, 282.2792]}
    metadata = lcb_metadata("d18:1", ids, masses)
    assert metadata[0].chain_composition == "d18:1"
    assert metadata[0].role == EvidenceRole.MOLECULAR_SPECIES_EVIDENCE
    assert metadata[1].theoretical_mz == 282.2792
    assert ids == {"d18:1": ["d18:1", "d18:1+OH"]}
    assert masses == {"d18:1": [264.2686, 282.2792]}
    with pytest.raises(ValueError, match="aligned"):
        lcb_metadata("d18:1", ids, {"d18:1": [264.2686]})


def test_support_cannot_be_marked_as_chain_evidence():
    with pytest.raises(ValueError, match="chain-specific"):
        FragmentMetadata("support", "test", EvidenceType.SUPPORT,
                         EvidenceRole.MOLECULAR_SPECIES_EVIDENCE, chain_composition="d18:1")


def test_neutral_loss_retains_origin_without_promoting_role():
    meta = neutral_loss_metadata("mock-Glc", "GSL_fragments", EvidenceType.GLYCAN)
    assert meta.evidence_type == EvidenceType.NEUTRAL_LOSS
    assert meta.loss_origin == EvidenceType.GLYCAN
    assert meta.role == EvidenceRole.SUPPORTING_EVIDENCE
