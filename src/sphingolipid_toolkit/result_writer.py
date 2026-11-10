"""Preserve complete MS2 evidence when a JSON cell exceeds Excel's limit."""
import gzip
import hashlib
import json
from pathlib import Path

from .io_utils import write_excel_atomic


def write_identifications(data,path):
    path=Path(path)
    output=data.copy()
    evidence={}
    external_indices=[]
    if 'ms2_evidence_json' in output:
        for index,value in output['ms2_evidence_json'].items():
            if isinstance(value,str) and len(value)>32767:
                key=hashlib.sha256(value.encode('utf-8')).hexdigest()
                evidence[key]=value
                output.loc[index,'ms2_evidence_json']=json.dumps({'external_evidence_sha256':key},separators=(',',':'))
                output.loc[index,'evidence_sha256']=key
                external_indices.append(index)
    created=None
    if evidence:
        payload=json.dumps(evidence,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode('utf-8')
        digest=hashlib.sha256(payload).hexdigest()
        sidecar=path.parent/f'ms2-evidence-{digest[:16]}.json.gz'
        path.parent.mkdir(parents=True,exist_ok=True)
        if sidecar.exists():
            if gzip.decompress(sidecar.read_bytes())!=payload:
                raise ValueError('An existing evidence file has different contents')
        else:
            with sidecar.open('xb') as handle:handle.write(gzip.compress(payload,mtime=0))
            created=sidecar
        mask=output.index.isin(external_indices)
        output.loc[mask,'evidence_file']=sidecar.name
    try:
        write_excel_atomic(output,path,index=False)
    except BaseException:
        if created is not None and created.resolve().parent==path.resolve().parent:
            created.unlink(missing_ok=True)
        raise
    return path
