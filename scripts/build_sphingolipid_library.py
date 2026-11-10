"""Freshly normalize a curated MSP into a separate, verified standalone index."""
from dataclasses import replace
import argparse
import json
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import sphingolipid_toolkit._vendor.lipidgate.ms2 as ms2
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_streaming_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import (
    LIBRARY_CACHE_VERSION,HYDROXY_FA_SPHINGOLIPID_CLASSES,iter_standard_msp_records,
    _without_d18_1_300_lcb_marker,_expand_hydroxy_fa_sphingolipids,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import CLASS_KEYS
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint,sha256


def build(source,destination):
    source,destination=map(Path,(source,destination))
    if destination.exists() or destination.suffix!='.sqlite':
        raise ValueError('Use a new, separate .sqlite destination for review')
    rules=code_fingerprint(Path(ms2.__file__).parent)
    before=sha256(source)
    allowed={re.sub('[^A-Z0-9]','',value.upper()) for value in CLASS_KEYS}|{'HEXIPC'}
    def records():
        seeds=[];next_id=0
        for record in iter_standard_msp_records(source):
            next_id=max(next_id,record.record_id+1)
            if re.sub('[^A-Z0-9]','',record.compound_class.upper()) not in allowed:continue
            record=_without_d18_1_300_lcb_marker(record)
            if record.compound_class in HYDROXY_FA_SPHINGOLIPID_CLASSES and record.adduct=='[M+H]+':
                seeds.append(record)
            yield record
        derived=_expand_hydroxy_fa_sphingolipids(seeds)[len(seeds):]
        for offset,record in enumerate(derived):yield replace(record,record_id=next_id+offset)
    metadata=write_streaming_library_index(destination,records(),dict(version=LIBRARY_CACHE_VERSION,
        source_sha256=before,rules_sha256=rules))
    if sha256(source)!=before or code_fingerprint(Path(ms2.__file__).parent)!=rules:
        raise ValueError('Source or rules changed during index verification')
    return dict(metadata=metadata,index_sha256=sha256(destination),index_bytes=destination.stat().st_size,
        fresh_normalization=True,blockwise_roundtrip_verified=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--audit',type=Path,required=True)
    args=parser.parse_args()
    if args.audit.exists():parser.error('Use a new audit filename')
    result=build(args.source,args.destination)
    args.audit.parent.mkdir(parents=True,exist_ok=True)
    args.audit.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Verified records:',result['metadata']['record_count'])
