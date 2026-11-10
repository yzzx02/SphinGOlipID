import gzip
import hashlib
import json

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
