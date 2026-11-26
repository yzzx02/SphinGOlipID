"""Build and check the positive HG8 library against the preceding curated release."""
import argparse
from collections import Counter
from dataclasses import replace
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3

from extend_verified_modified_sphingoid_library import rebuild_verified_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.headgroup_core import CURATED, GM2_B, REVISION, curate_headgroup
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_ganglioside import NH3
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

PEAK = re.compile(rb'([^\s]+)\s+([^\s]+)\s+"([^"]+)"\s+"([^"]+)"$')
UPDATED = {c.encode() for c in CURATED}

@lru_cache(maxsize=65536)
def peak_role(cls, adduct, name, kind):
    if cls=='GM2' and name=='B1alpha':
        return name,kind
    f=FragmentRecord(100.,name,kind)
    result,=curate_headgroup(cls,adduct,[f])
    return result.name,result.fragment_type

def curate_block(raw,cls,adduct):
    header, marker, peaks=raw.partition(b'Num Peaks: ')
    if not marker: raise ValueError("Missing MSP peak count")
    n,_,peakdata=peaks.partition(b'\n')
    lines=peakdata.splitlines()
    if len(lines)!=int(n): raise ValueError("Invalid MSP peak count")
    revised=[]
    for line in lines:
        m=PEAK.fullmatch(line)
        if not m:raise ValueError(line)
        mz,intensity,name,kind=m.groups()
        new_name,new_kind=peak_role(cls,adduct,name.decode(),kind.decode())
        if (new_name.encode(),new_kind.encode())!=(name,kind):
            line=mz+b' '+intensity+b' "'+new_name.encode()+b'" "'+new_kind.encode()+b'"'
        revised.append(line)
    if cls=='GM2':
        revised.append(b'274.0921 100.00 "B1alpha-H2O" "Diagnostic_B"')
    return header+marker+str(len(revised)).encode()+b'\n'+b'\n'.join(revised)

def expected(record):
    if record.compound_class not in CURATED or record.polarity!='+':return record
    return replace(record,fragments=curate_headgroup(record.compound_class,record.adduct,record.fragments))

def build(source,old_index,directory):
    with sqlite3.connect(old_index.resolve().as_uri()+'?mode=ro',uri=True) as connection:
        prior=json.loads(connection.execute("SELECT value FROM metadata WHERE key='library'").fetchone()[0])
    if prior['revision']==REVISION or prior['source_sha256']!=sha256(source):
        raise ValueError('Supply the matching MSP and index from the preceding library release.')
    directory.mkdir(exist_ok=True,parents=True)
    output=directory/'sphingolipid_positive.msp.gz'
    residual=b'';counts=Counter();added=None;content=hashlib.sha256()
    with gzip.open(source,'rb') as reader,output.open('xb') as target, gzip.GzipFile(filename='',fileobj=target,mode='wb',mtime=0,compresslevel=5) as writer:
        def emit(data):
            writer.write(data);content.update(data)
        for chunk in iter(lambda:reader.read(4<<20),b''):
            blocks=(residual+chunk).split(b'\n\n');residual=blocks.pop()
            for block in blocks:
                cls=block.partition(b'CompoundClass: ')[2].partition(b'\n')[0]
                if cls in UPDATED:
                    adduct=block.partition(b'PrecursorType: ')[2].partition(b'\n')[0].decode()
                    block=curate_block(block,cls.decode(),adduct)
                    if cls==b'GM2' and block.startswith(b'Name: GM2(d18:1/18:0)\n') and adduct=='[M+H]+':
                        if added is not None:raise ValueError("Duplicate source template")
                        _,_,peaks=block.partition(b'Num Peaks: ')
                        products=[l for l in peaks.splitlines()[1:] if b'"LCB-' not in l]
                        mh=float(block.partition(b'PrecursorMZ: ')[2].partition(b'\n')[0])
                        products=[l.replace(b'"[M+H]+" "Precursor Ion"',b'"M+NH4-NH3" "NL"') for l in products]
                        mass=mh+NH3
                        products.append(f'{mass:.8f} 100.00 "[M+NH4]+" "Precursor Ion"'.encode())
                        hdr=f'Name: GM2(d36:1)\nPrecursorMZ: {mass:.8f}\nPrecursorType: [M+NH4]+\nCompoundClass: GM2\nFormula: C67H121N3O26\nComment: MS1_name=GM2(d36:1);polarity=+\nLibraryRevision: {REVISION}\nAnnotationScope: sum composition only\nNum Peaks: {len(products)}\n'
                        added=hdr.encode()+b'\n'.join(products)+b'\n\n'
                    counts[cls.decode()]+=1
                emit(block+b'\n\n')
        if residual.strip() or added is None:raise ValueError("Incomplete source library")
        emit(added)
    print('Curated MSP',counts,flush=True)
    audit=rebuild_verified_index(output,old_index,directory/'positive.sqlite',expected_transform=expected,
        revision=REVISION,additional_identities={('GM2','GM2(d36:1)','C67H121N3O26','[M+NH4]+')})
    audit.update(updated_classes=counts,raw_content_sha256=content.hexdigest())
    (directory/'positive.audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')

if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ('source','previous_index','destination'):cli.add_argument(key,type=Path)
    args=cli.parse_args()
    build(args.source,args.previous_index,args.destination)
