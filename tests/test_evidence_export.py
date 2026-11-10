import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from sphingolipid_toolkit import result_writer


def test_excel_overflow_retains_exact_evidence_and_existing_references(tmp_path):
    long=json.dumps({'spectrum':[[123.456789,100.123456]]*3000},separators=(',',':'))
    old='{"external_evidence_sha256":"old"}'
    data=pd.DataFrame({'ms2_evidence_json':[long,old],'evidence_sha256':[None,'old'],
        'evidence_file':[None,'prior.json.gz'],'total_score':[88.1234,76.5],'rank':[1,2]})
    path=tmp_path/'result.xlsx'
    result_writer.write_identifications(data,path)
    exported=pd.read_excel(path)
    assert exported.loc[1,'evidence_file']=='prior.json.gz'
    key=hashlib.sha256(long.encode()).hexdigest()
    assert exported.loc[0,'evidence_sha256']==key
    complete=json.loads(gzip.decompress((tmp_path/exported.loc[0,'evidence_file']).read_bytes()))
    assert complete[key]==long
    assert exported.total_score.tolist()==[88.1234,76.5]
    assert data.loc[0,'ms2_evidence_json']==long


def test_failed_excel_export_keeps_previous_result_and_cleans_new_sidecar(tmp_path,monkeypatch):
    path=tmp_path/'result.xlsx';path.write_bytes(b'previous result')
    data=pd.DataFrame({'ms2_evidence_json':[json.dumps({'spectrum':'x'*40000})]})
    def fail(*args,**kwargs):raise PermissionError('workbook open')
    monkeypatch.setattr(result_writer,'write_excel_atomic',fail)
    with pytest.raises(PermissionError):result_writer.write_identifications(data,path)
    assert path.read_bytes()==b'previous result'
    assert not list(tmp_path.glob('ms2-evidence-*.json.gz'))


def test_save_bundle_preserves_nested_references_in_all_worksheets(tmp_path):
    source = tmp_path / 'source' / 'result.xlsx'
    source.parent.mkdir()
    paths = ['a/evidence.json.gz', 'b/evidence.json.gz']
    for index, name in enumerate(paths):
        sidecar = source.parent / name
        sidecar.parent.mkdir()
        sidecar.write_bytes(gzip.compress(json.dumps({'key': str(index)}).encode()))
    with pd.ExcelWriter(source) as writer:
        for index, name in enumerate(paths):
            pd.DataFrame({'evidence_file': [name], 'total_score': [91.2345],
                          'rank': [index + 1]}).to_excel(writer, sheet_name=f'part{index}', index=False)
    destination = tmp_path / 'saved' / 'result.xlsx'
    result_writer.copy_identification_bundle(source, destination)
    assert destination.read_bytes() == source.read_bytes()
    for name in paths:
        assert (destination.parent / name).read_bytes() == (source.parent / name).read_bytes()


@pytest.mark.parametrize('absolute', [True, False])
def test_save_bundle_localizes_absolute_or_parent_references(tmp_path, absolute):
    evidence = tmp_path / 'shared.json.gz'
    evidence.write_bytes(gzip.compress(b'{"key":"complete evidence"}'))
    source = tmp_path / 'source' / 'result.xlsx'
    source.parent.mkdir()
    reference = str(evidence) if absolute else '../shared.json.gz'
    frame = pd.DataFrame({'evidence_file': [reference], 'evidence_sha256': ['preserve-me'],
                          'total_score': [88.7654], 'rank': [2]})
    with pd.ExcelWriter(source) as writer:
        frame.to_excel(writer, sheet_name='annotations', index=False)
        pd.DataFrame({'parameter': ['keep me']}).to_excel(writer, sheet_name='metadata', index=False)
    original = source.read_bytes()
    destination = tmp_path / 'saved' / 'result.xlsx'
    result_writer.copy_identification_bundle(source, destination)
    result = pd.read_excel(destination)
    localized = Path(result.loc[0, 'evidence_file'])
    assert not localized.is_absolute() and '..' not in localized.parts
    assert (destination.parent / localized).read_bytes() == evidence.read_bytes()
    pd.testing.assert_frame_equal(result.drop(columns='evidence_file'), frame.drop(columns='evidence_file'))
    assert pd.read_excel(destination, sheet_name='metadata').parameter.tolist() == ['keep me']
    assert source.read_bytes() == original


def test_bundle_publish_failure_retains_previous_result_and_removes_new_sidecar(tmp_path, monkeypatch):
    source = tmp_path / 'source' / 'result.xlsx'
    source.parent.mkdir()
    result_writer.write_identifications(pd.DataFrame({'ms2_evidence_json': ['x' * 40000]}), source)
    destination = tmp_path / 'saved' / 'result.xlsx'
    destination.parent.mkdir()
    destination.write_bytes(b'previous result')
    replace = Path.replace
    def fail_publish(path, target):
        if Path(target) == destination:
            raise PermissionError('destination workbook is open')
        return replace(path, target)
    monkeypatch.setattr(Path, 'replace', fail_publish)
    with pytest.raises(PermissionError, match='workbook is open'):
        result_writer.copy_identification_bundle(source, destination)
    assert destination.read_bytes() == b'previous result'
    assert list(destination.parent.iterdir()) == [destination]


def test_repeated_evidence_references_are_hashed_once_per_source(tmp_path, monkeypatch):
    source = tmp_path / 'source' / 'result.xlsx'
    source.parent.mkdir()
    sidecar = source.parent / 'evidence.json.gz'
    sidecar.write_bytes(gzip.compress(b'{"key":"complete evidence"}'))
    pd.DataFrame({'evidence_file': [sidecar.name] * 100}).to_excel(source, index=False)
    calls = []
    digest = result_writer._file_digest
    def count(path):
        calls.append(path)
        return digest(path)
    monkeypatch.setattr(result_writer, '_file_digest', count)
    result_writer.copy_identification_bundle(source, tmp_path / 'saved' / 'result.xlsx')
    assert calls.count(sidecar.resolve()) == 1


def test_multisheet_export_preserves_complete_distinct_evidence_and_numeric_scores(tmp_path):
    first=json.dumps({'spectrum':'x'*40000},separators=(',',':'))
    second=json.dumps({'spectrum':'y'*40000},separators=(',',':'))
    tables={'annotations':pd.DataFrame({'ms2_evidence_json':[first],'total_score':[88.1234],'rank':[1]}),
        'spectrum_matches':pd.DataFrame({'ms2_evidence_json':[first,second],'total_score':[88.1234,77.6],'rank':[1,2]}),
        'features':pd.DataFrame({'xcms_peak_id':['CP1'],'peak_area':[123456.789]})}
    path=tmp_path/'result.xlsx'
    result_writer.write_identification_workbook(tables,path)
    selected=pd.read_excel(path,sheet_name='annotations')
    all_matches=pd.read_excel(path,sheet_name='spectrum_matches')
    assert selected.evidence_file.iloc[0]==all_matches.evidence_file.iloc[1]
    payload=json.loads(gzip.decompress((path.parent/selected.evidence_file.iloc[0]).read_bytes()))
    assert set(payload.values())=={first,second}
    assert all_matches.total_score.tolist()==[88.1234,77.6]
    assert all_matches['rank'].tolist()==[1,2]
    assert tables['annotations'].ms2_evidence_json.iloc[0]==first


def test_failed_multisheet_export_keeps_previous_output_and_cleans_artifacts(tmp_path,monkeypatch):
    path=tmp_path/'result.xlsx'
    path.write_bytes(b'previous result')
    def fail(*args,**kwargs):raise PermissionError('workbook open')
    monkeypatch.setattr(result_writer.pd,'ExcelWriter',fail)
    with pytest.raises(PermissionError):
        result_writer.write_identification_workbook({'annotations':pd.DataFrame({'ms2_evidence_json':['x'*40000]})},path)
    assert path.read_bytes()==b'previous result'
    assert list(tmp_path.iterdir())==[path]
