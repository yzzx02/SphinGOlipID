import pytest

from scripts.build_sphingolipid_library import build
from sphingolipid_toolkit.engine import IdentificationConfig,SphingolipidIdentifier
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import write_standard_msp
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord,LibraryRecord


def test_standalone_fresh_library_builder_is_searchable_and_roundtrips(tmp_path):
    record=LibraryRecord(0,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',538.5199,'[M+H]+',polarity='+',
        fragments=[FragmentRecord(264.2686,'LCB-2H2O','LCB碎片'),FragmentRecord(282.2791,'LCB-H2O','LCB碎片')])
    source=write_standard_msp([record],tmp_path/'curated.msp')
    output=tmp_path/'new.sqlite'
    audit=build(source,output)
    assert audit['fresh_normalization'] and audit['blockwise_roundtrip_verified']
    with SphingolipidIdentifier(IdentificationConfig(library_path=output)) as identifier:
        rows=identifier.identify(538.5199,[(264.2686,100.)])
        assert rows and rows[0]['matched_name']=='Cer(d18:1/16:0)'


def test_standalone_library_builder_preserves_existing_destinations(tmp_path):
    output=tmp_path/'existing.sqlite';output.write_bytes(b'previous checkpoint')
    with pytest.raises(ValueError):build(tmp_path/'not_read.msp',output)
    assert output.read_bytes()==b'previous checkpoint'
