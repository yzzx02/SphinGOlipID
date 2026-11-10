"""Reference-engine behavior independent of a LipidGate installation."""
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from sphingolipid_toolkit.engine import IdentificationConfig,SphingolipidIdentifier
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord,LibraryRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint


@pytest.fixture
def cer_index(tmp_path):
    record=LibraryRecord(7,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',538.5199,'[M+H]+',
        'C34H67NO3','+',[
        FragmentRecord(264.2686,'LCB-2H2O','LCB碎片'),
        FragmentRecord(282.2791,'LCB-H2O','LCB碎片'),
        FragmentRecord(520.5093,'M+H-H2O','Common'),
        FragmentRecord(538.5199,'[M+H]+','Precursor Ion'),
    ])
    path=tmp_path/'cer.sqlite'
    import sphingolipid_toolkit._vendor.lipidgate.ms2 as ms2
    write_library_index(path,[record],dict(version=LIBRARY_CACHE_VERSION,
        source_sha256='synthetic',rules_sha256=code_fingerprint(Path(ms2.__file__).parent)))
    return path


def test_single_lcb_admission_and_ppm_setting_reach_reference_core(cer_index):
    cfg=IdentificationConfig(library_path=cer_index,fragment_ppm=20.)
    with SphingolipidIdentifier(cfg) as identifier:
        rows=identifier.identify(538.5199,[(264.2686*(1+15e-6),100.)])
        assert len(rows)==1
        assert rows[0]['passed_required_gates']
        assert rows[0]['matched_fragment_count']==1
        assert rows[0]['matched_name']=='Cer(d18:1/16:0)'
    with SphingolipidIdentifier(replace(cfg,fragment_ppm=10.)) as identifier:
        assert identifier.identify(538.5199,[(264.2686*(1+15e-6),100.)])==[]


def test_auxiliary_lcb_support_retains_raw_audit_and_pool_scores(cer_index):
    import json
    cfg=IdentificationConfig(library_path=cer_index)
    with SphingolipidIdentifier(cfg) as identifier:
        one=identifier.identify(538.5199,[(264.2686,100.),(520.5093,10.)])[0]
        two=identifier.identify(538.5199,[(264.2686,100.),(282.2791,10.),(520.5093,10.)])[0]
    assert two['lcb_score']==one['lcb_score']
    assert two['other_score']>one['other_score']
    audit=json.loads(two['ms2_evidence_json'])
    auxiliary=[f for f in audit['fragments'] if f.get('auxiliary_lcb_support')]
    assert len(auxiliary)==1 and auxiliary[0]['name']=='LCB-H2O'


def test_text_production_uses_reference_engine_without_legacy_filters(cer_index,tmp_path):
    from sphingolipid_toolkit.config import SphinGOlipIDConfig
    from sphingolipid_toolkit.ms2_pipeline import run_batch
    raw=tmp_path/'raw';raw.mkdir()
    (raw/'hilic-msms-1.txt').write_text('\n'.join([
        'spectrum:','index:1','scan start time,1.25,min','target m/z,538.5199, m/z',
        'binaryDataArray:','binary: 0 264.2686','binaryDataArray:','binary: 0 100.0',
    ]),encoding='utf-8')
    pd.DataFrame(dict(target=[538.5199],下限=[538.51],上限=[538.53],RT=[1.25],Abund=[1000.])).to_excel(raw/'hilic-msms-1.xlsx',index=False)
    config=SphinGOlipIDConfig(raw,raw,None,tmp_path/'out',file_indices=(1,),ms2_library_path=cer_index)
    result=pd.read_excel(run_batch(config)[0])
    assert len(result)==1
    assert result.annotation_engine.iloc[0]=='SphinGOlipID'
    assert result.matched_fragment_count.iloc[0]==1
    assert result.total_score.iloc[0]>=50
    assert 'ms2_evidence_json' in result and 'lcb_score' in result


@pytest.mark.parametrize('kwargs',[{'fragment_ppm':0},{'fragment_da':.3},{'min_total_score':101},{'top_n':0}])
def test_reference_configuration_rejects_ambiguous_or_invalid_settings(kwargs):
    with pytest.raises(ValueError):IdentificationConfig(**kwargs)


def test_bundled_library_rejects_a_changed_payload(tmp_path,monkeypatch):
    import hashlib,json
    import sphingolipid_toolkit.engine as engine
    root=tmp_path/'data/ms2';root.mkdir(parents=True)
    bank=root/'sphingolipid_positive.sqlite';bank.write_bytes(b'reference fixture')
    (root.parent/'library_revision.json').write_text(json.dumps({'positive':{
        'index_sha256':hashlib.sha256(bank.read_bytes()).hexdigest()}}),encoding='utf-8')
    monkeypatch.setattr(engine,'files',lambda _:tmp_path)
    assert engine.library_path('+')==bank
    bank.write_bytes(b'changed fixture')
    with pytest.raises(ValueError,match='fingerprint mismatch'):engine.library_path('+')


def write_native_mzml(path,with_ms1=True):
    import base64,struct
    def array(values,name,accession):
        binary=base64.b64encode(struct.pack('<'+'d'*len(values),*values)).decode()
        return (f'<binaryDataArray encodedLength="{len(binary)}">'
            '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
            '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>'
            f'<cvParam cvRef="MS" accession="{accession}" name="{name}"/>'
            f'<binary>{binary}</binary></binaryDataArray>')
    def scan(index,level,rt,masses,intensities,precursor=''):
        return (f'<spectrum index="{index-1}" id="scan={index}" defaultArrayLength="{len(masses)}">'
            f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>'
            '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>'
            '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>'
            '<scanList count="1"><scan>'
            f'<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt}" '
            'unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>'
            '</scan></scanList>'+precursor+'<binaryDataArrayList count="2">'
            +array(masses,'m/z array','MS:1000514')+array(intensities,'intensity array','MS:1000515')
            +'</binaryDataArrayList></spectrum>')
    raw=538.5199*(1+12e-6);measured=538.5199*(1+5e-6)
    precursor=('<precursorList count="1"><precursor spectrumRef="scan=1">'
        '<selectedIonList count="1"><selectedIon>'
        f'<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{raw}"/>'
        '<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="1"/>'
        '</selectedIon></selectedIonList><activation>'
        '<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="60"/>'
        '</activation></precursor></precursorList>')
    spectra=[]
    if with_ms1:spectra.append(scan(1,1,60.,[measured],[10000.]))
    spectra.append(scan(2,2,61.,[264.2686],[1000.],precursor))
    if with_ms1:spectra.append(scan(3,1,62.,[measured],[9000.]))
    path.write_text('<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">'
        f'<run id="synthetic"><spectrumList count="{len(spectra)}">'+''.join(spectra)
        +'</spectrumList></run></mzML>',encoding='utf-8')
    return path


def test_native_mzml_cli_uses_measured_ms1_and_preserves_scan_and_energy(cer_index,tmp_path):
    from sphingolipid_toolkit.identify_cli import main
    path=write_native_mzml(tmp_path/'reference.mzML')
    output=tmp_path/'identified.xlsx'
    assert main([str(path),'--output',str(output),'--library',str(cer_index),'--top-n','1'])==0
    result=pd.read_excel(output)
    assert len(result)==1
    row=result.iloc[0]
    assert row.precursor_mz_source=='MS1_REFERENCE_CENTROID'
    assert row.precursor_mz==pytest.approx(538.5199*(1+5e-6))
    assert row.precursor_mz_raw==pytest.approx(538.5199*(1+12e-6))
    assert row.native_scan_id=='scan=2' and row.collision_energy==60
    assert row.matched_name=='Cer(d18:1/16:0)' and row.total_score==60


def test_native_mzml_does_not_replace_missing_ms1_with_theoretical_mass(cer_index,tmp_path):
    path=write_native_mzml(tmp_path/'missing_ms1.mzML',with_ms1=False)
    with SphingolipidIdentifier(IdentificationConfig(library_path=cer_index)) as identifier:
        assert identifier.identify_mzml(path).empty
