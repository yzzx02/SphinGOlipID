"""Archived data are readable; their classifier cannot annotate new spectra."""
import pytest
from sphingolipid_toolkit import evidence
from sphingolipid_toolkit.backend_policy import RemovedBackendError

@pytest.mark.parametrize('name',[
    'classify_fragment','capture_fragment_origins','record_lcb_origin',
    'fragment_source','scoring_pool_for_evidence','annotate_theoretical_fragments'])
def test_archived_classification_routes_are_retired(name):
    with pytest.raises(RemovedBackendError):getattr(evidence,name)()

def test_archived_evidence_metadata_remains_serializable():
    record=evidence.EvidenceRecord(evidence.EvidenceType.HG,
        evidence.EvidenceRole.CLASS_DIAGNOSTIC,True,'diagnostic')
    assert record.to_dict()['fragment_type']=='HG'
    assert record.to_dict()['evidence_role']=='class_diagnostic'
