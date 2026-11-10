"""Append the native-CID DMS reference and freshly rebuild with stable source IDs.

The previous index supplies identities for unchanged records, never cached
normalization or search plans. Every source block and every derived record is
freshly normalized, compared in full, then reindexed with the current rules.
"""
from __future__ import annotations

from dataclasses import replace
from contextlib import closing
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import re
import sqlite3
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import sphingolipid_toolkit._vendor.lipidgate.ms2 as ms2
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import IndexedLibrary,write_streaming_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import (
    LIBRARY_CACHE_VERSION,HYDROXY_FA_SPHINGOLIPID_CLASSES,iter_standard_msp_records,
    _without_d18_1_300_lcb_marker,_expand_hydroxy_fa_sphingolipids,_encode_fragment_fields,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_legacy_sphingolipid import fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint,sha256

REVISION='sphingolipid-literature-20261110-r1'
DMS_CLASS='N,N-dimethylSo'
DMS_NAME='N,N-dimethylSo(d18:1)'
DMS_FORMULA='C20H41NO2'
PARENT_PATTERN=re.compile(r'^[^()]+\(d\d+:\d+/\d+:\d+\)$',re.I)


def scientific_payload(record):
    return [record.compound_class,record.lipid_name,record.lipid_chain_name,float(record.precursor_mz),
        record.adduct,record.formula,record.polarity,
        [[float(f.mz),f.name,f.fragment_type,float(f.intensity),float(f.weight),f.required_group]
         for f in record.fragments],record.metadata]


def signature(record):
    return hashlib.sha256(json.dumps(scientific_payload(record),sort_keys=True,
        separators=(',',':'),ensure_ascii=False).encode('utf-8')).digest()


def composition_key(record):
    return json.dumps([record.compound_class,record.lipid_chain_name,
        round(float(record.precursor_mz),4),record.adduct],separators=(',',':'))


def append_dms(source,baseline_index,destination):
    source,baseline_index,destination=map(Path,(source,baseline_index,destination))
    if destination.exists() or destination.resolve() in {source.resolve(),baseline_index.resolve()}:
        raise ValueError('Use a new, separate MSP destination')
    parent=None
    expected_parent=18*C+37*H+N+2*O+PROTON
    with closing(IndexedLibrary(baseline_index,allowed_class_keys={'SPB','NNDIMETHYLSO'})) as library:
        for record in library:
            if record.compound_class==DMS_CLASS:
                raise ValueError('DMS reference already present; do not append twice')
            if (record.compound_class=='SPB' and record.lipid_chain_name=='SPB(d18:1)'
                    and record.adduct=='[M+H]+' and record.formula in {'','C18H37NO2'}
                    and abs(record.precursor_mz-expected_parent)/expected_parent*1e6<=10.):
                parent=record
    if parent is None:raise ValueError('Missing eligible d18:1 SPB parent')
    mass=20*C+41*H+N+2*O+PROTON
    products=fragments(DMS_CLASS,mass,2)
    lines=[f'Name: {DMS_NAME}',f'PrecursorMZ: {mass:.8f}','PrecursorType: [M+H]+',
        f'CompoundClass: {DMS_CLASS}',f'Formula: {DMS_FORMULA}',
        f'Comment: MS1_name={DMS_NAME};polarity=+',f'LibraryRevision: {REVISION}',
        'Source: native positive CID standard; doi:10.1038/nprot.2013.004; Fig.7b',
        'SpectrumOrigin: literature-derived elemental masses; theoretical equal intensities; not measured spectrum',
        f'EligibleParent: {parent.lipid_chain_name}; previous_record_id={parent.record_id}',
        f'Num Peaks: {len(products)}']
    lines.extend(f'{f.mz:.4f} {f.intensity:.2f} {_encode_fragment_fields(f)}' for f in products)
    addition=('\n'.join(lines)+'\n\n').encode('utf-8')
    before=sha256(source);prefix=hashlib.sha256();count=0;tail=b''
    destination.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(source,'rb') as old,destination.open('wb') as output:
        with gzip.GzipFile(filename='',fileobj=output,mode='wb',compresslevel=9,mtime=0) as new:
            for chunk in iter(lambda:old.read(1024*1024),b''):
                new.write(chunk);prefix.update(chunk);count+=len(chunk);tail=chunk[-2:]
            if tail!=b'\n\n':new.write(b'\n\n')
            new.write(addition)
    check=hashlib.sha256()
    with gzip.open(destination,'rb') as handle:
        remaining=count
        while remaining:
            chunk=handle.read(min(remaining,1024*1024));assert chunk
            check.update(chunk);remaining-=len(chunk)
    if sha256(source)!=before or check.digest()!=prefix.digest():
        raise ValueError('An existing source byte changed')
    return dict(added_records=1,added_name=DMS_NAME,source_sha256=before,
        preserved_decompressed_prefix_bytes=count,preserved_decompressed_prefix_sha256=prefix.hexdigest(),
        destination_sha256=sha256(destination),eligible_parent_record_id=parent.record_id)


def rebuild_verified_index(source,baseline_index,destination,*,allow_dms=False):
    source,baseline_index,destination=map(Path,(source,baseline_index,destination))
    if destination.exists() or destination.resolve() in {source.resolve(),baseline_index.resolve()}:
        raise ValueError('Use a new, separate SQLite destination')
    before=sha256(source);baseline_hash=sha256(baseline_index)
    rules=code_fingerprint(Path(ms2.__file__).parent)
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent,prefix='fresh_scientific_') as temporary, \
            closing(sqlite3.connect(Path(temporary)/'verification.sqlite')) as stage:
        stage.execute('PRAGMA journal_mode=OFF');stage.execute('PRAGMA synchronous=OFF')
        stage.execute('PRAGMA temp_store=FILE');stage.execute('PRAGMA cache_size=-8192')
        stage.execute('CREATE TABLE old (ordinal INTEGER PRIMARY KEY,sig BLOB,id INTEGER,name TEXT,payload BLOB)')
        stage.execute('CREATE INDEX signature_lookup ON old(sig,ordinal) WHERE payload IS NULL')
        stage.execute('CREATE TABLE keys (value TEXT PRIMARY KEY)')
        stage.execute('CREATE TABLE parents (ordinal INTEGER PRIMARY KEY,payload BLOB)')
        old_count=0;maximum_id=-1
        with closing(IndexedLibrary(baseline_index)) as library:
            for position,record in enumerate(library):
                stage.execute('INSERT INTO old(ordinal,sig,id,name) VALUES (?,?,?,?)',
                    (position,signature(record),record.record_id,record.lipid_chain_name))
                maximum_id=max(maximum_id,record.record_id);old_count+=1
                if old_count%25000==0:stage.commit()
        stage.commit();print('Previous identities verified:',old_count,flush=True)
        matched=0;raw=0;derived=0;additions=[]

        def assign(record):
            nonlocal matched
            existing=stage.execute('SELECT ordinal,id FROM old WHERE sig=? AND payload IS NULL ORDER BY ordinal LIMIT 1',
                (signature(record),)).fetchone()
            if existing is None:
                if (allow_dms and record.compound_class==DMS_CLASS and record.lipid_chain_name==DMS_NAME
                        and record.formula==DMS_FORMULA and record.adduct=='[M+H]+' and not additions):
                    additions.append(replace(record,record_id=maximum_id+1));return
                raise ValueError('Fresh normalization differs from previous content: '+record.lipid_chain_name)
            position,identity=existing
            fresh=replace(record,record_id=identity)
            stage.execute('UPDATE old SET payload=? WHERE ordinal=?',
                (pickle.dumps(fresh,protocol=5),position));matched+=1

        for record in iter_standard_msp_records(source):
            record=_without_d18_1_300_lcb_marker(record)
            assign(record);raw+=1
            stage.execute('INSERT OR IGNORE INTO keys VALUES (?)',(composition_key(record),))
            if (record.compound_class in HYDROXY_FA_SPHINGOLIPID_CLASSES and record.adduct=='[M+H]+'
                    and PARENT_PATTERN.fullmatch(record.lipid_chain_name)):
                stage.execute('INSERT INTO parents VALUES (?,?)',(raw,pickle.dumps(record,protocol=5)))
            if raw%25000==0:
                stage.commit();print('Fresh source records:',raw,flush=True)
        stage.commit()
        # Same shared expansion, evaluated one parent at a time. The complete
        # source-key set preserves the original global duplicate suppression.
        for parent_row, in stage.execute('SELECT payload FROM parents ORDER BY ordinal'):
            parent=pickle.loads(parent_row)
            for candidate in _expand_hydroxy_fa_sphingolipids([parent])[1:]:
                inserted=stage.execute('INSERT OR IGNORE INTO keys VALUES (?)',(composition_key(candidate),)).rowcount
                if not inserted:continue
                assign(candidate);derived+=1
        stage.commit()
        missing=stage.execute('SELECT id,name FROM old WHERE payload IS NULL LIMIT 3').fetchall()
        if matched!=old_count or missing:raise ValueError(f'Missing previous records after fresh normalization: {missing}')
        if len(additions)!=(1 if allow_dms else 0):raise ValueError('Unexpected appended reference count')

        def records():
            for payload, in stage.execute('SELECT payload FROM old ORDER BY ordinal'):
                yield pickle.loads(payload)
            yield from additions

        metadata=write_streaming_library_index(destination,records(),dict(version=LIBRARY_CACHE_VERSION,
            source_sha256=before,rules_sha256=rules,revision=REVISION,
            previous_index_sha256=baseline_hash,fresh_normalization=True))
        with closing(IndexedLibrary(baseline_index)) as original,closing(IndexedLibrary(destination)) as revised:
            previous=iter(original);preserved=0;new_ids={r.record_id for r in additions}
            for record in revised:
                if record.record_id in new_ids:continue
                old=next(previous)
                if old.record_id!=record.record_id or scientific_payload(old)!=scientific_payload(record):
                    raise ValueError('Source IDs, scientific fields or relative order changed')
                preserved+=1
            if next(previous,None) is not None or preserved!=old_count:raise ValueError('Incomplete index roundtrip')
    if sha256(source)!=before or sha256(baseline_index)!=baseline_hash or code_fingerprint(Path(ms2.__file__).parent)!=rules:
        raise ValueError('Source, baseline or rules changed during verification')
    return dict(metadata=metadata,index_sha256=sha256(destination),index_bytes=destination.stat().st_size,
        raw_records_freshly_normalized=raw,derived_records_freshly_generated=derived,
        preserved_original_records=old_count,preserved_ids_order_and_all_scientific_fields=True,
        rebuilt_search_plans=True,blockwise_roundtrip_verified=True,new_record_ids=[r.record_id for r in additions])


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['append-dms','rebuild'])
    for name in ['source','previous-index','destination','audit']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--allow-dms',action='store_true',help='Allow exactly the single reviewed DMS addition during rebuild.')
    args=parser.parse_args()
    if args.audit.exists():parser.error('Use a new audit path')
    result=(append_dms(args.source,args.previous_index,args.destination) if args.operation=='append-dms'
        else rebuild_verified_index(args.source,args.previous_index,args.destination,allow_dms=args.allow_dms))
    args.audit.parent.mkdir(parents=True,exist_ok=True)
    args.audit.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Verified:',args.destination)
