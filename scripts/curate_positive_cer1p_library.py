"""Curate positive Cer1P and append sodium total compositions in a separate MSP."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value,iter_blocks,parse_peaks,rebuild_block,validate_block
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_cer1p import Cer1PSodiumComposition,normalize_h_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import explicit_acyl_carbons, chain_double_bond_violations
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

REVISION="positive_cer1p_pools_20261108"
CHAIN=re.compile(r"Cer1P\(([mdt])(\d+):(\d+)/(\d+):(\d+)\)")


def is_protonated(block):
    return (header_value(block,"CompoundClass"),header_value(block,"PrecursorType"))==("Cer1P","[M+H]+")


def curate_block(block):
    if not is_protonated(block):return block
    validate_block(block)
    fragments=[FragmentRecord(p['mz'],p['name'],p['type'],intensity=p['intensity']) for p in parse_peaks(block)]
    revised=normalize_h_fragments(header_value(block,'Name'),float(header_value(block,'PrecursorMZ')),fragments)
    result=rebuild_block(block,[dict(mz=f.mz,name=f.name,type=f.fragment_type,intensity=f.intensity) for f in revised])
    validate_block(result)
    return result


def sodium_composition(block):
    if not is_protonated(block):return None
    name=header_value(block,'Name');m=CHAIN.fullmatch(name)
    if m is None or any(c<14 for c in explicit_acyl_carbons(name)) or chain_double_bond_violations(name,'Cer1P'):return None
    s,lc,ld,fc,fd=m.groups()
    return Cer1PSodiumComposition(int(lc)+int(fc),int(ld)+int(fd),{'m':1,'d':2,'t':3}[s])


def sodium_block(c):
    fragments=c.fragments()
    result=[f'Name: {c.name}',f'PrecursorMZ: {c.precursor_mz:.8f}','PrecursorType: [M+Na]+','CompoundClass: Cer1P',
            f'Formula: {c.formula}',f'Comment: MS1_name={c.name};polarity=+','Source: theoretical',
            f'LibraryRevision: {REVISION}','AnnotationScope: sum_composition_only',f'Num Peaks: {len(fragments)}']
    result += [f'{f.mz:.4f} {f.intensity:.2f} "{f.name}" "{f.fragment_type}"' for f in fragments]
    validate_block(result);return result


def curate(source,destination):
    source,destination=Path(source),Path(destination)
    if source.resolve()==destination.resolve():raise ValueError('Write to a separate destination before installing')
    destination.parent.mkdir(parents=True,exist_ok=True)
    compositions=set();existing=set();changed=Counter();preserved=0;source_count=0;digest=hashlib.sha256()
    with (destination.open('wb') as binary,
          gzip.GzipFile(filename='',fileobj=binary,mode='wb',compresslevel=9,mtime=0) as compressed,
          io.TextIOWrapper(compressed,encoding='utf-8',newline='\n') as output):
        for block in iter_blocks(source):
            source_count+=1
            if (c:=sodium_composition(block)) is not None:compositions.add(c)
            if (header_value(block,'CompoundClass'),header_value(block,'PrecursorType'))==('Cer1P','[M+Na]+'):
                existing.add(header_value(block,'Name'))
            revised=curate_block(block)
            if revised!=block:changed[CHAIN.fullmatch(header_value(block,'Name'))[1]]+=1
            else:preserved+=1;digest.update(('\n'.join(block)+'\n\n').encode('utf-8'))
            output.write('\n'.join(revised)+'\n\n')
        additions=[sodium_block(c) for c in sorted(compositions) if c.name not in existing]
        for block in additions:output.write('\n'.join(block)+'\n\n')
    revised=iter(iter_blocks(destination))
    for before in iter_blocks(source):
        after=next(revised,None)
        if after!=curate_block(before):raise ValueError('Unexpected prior record or order change')
        if is_protonated(before):
            old={p['name']:p for p in parse_peaks(before)};new={p['name']:p for p in parse_peaks(after)}
            assert all(n in new and all(new[n][k]==p[k] for k in ('mz','intensity')) for n,p in old.items())
            assert all(header_value(before,k)==header_value(after,k) for k in ('Name','PrecursorMZ','Formula','Comment'))
    for block in additions:
        if next(revised,None)!=block:raise ValueError('Unexpected sodium addition')
    if next(revised,None) is not None:raise ValueError('Unexpected trailing library record')
    return dict(revision=REVISION,changed_protonated_by_series=dict(changed),added_sodium_records=len(additions),
                added_sodium_by_total_hydroxyls=dict(Counter(c.hydroxyls for c in compositions if c.name not in existing)),
                original_source_records=source_count,preserved_records=preserved,preserved_records_sha256=digest.hexdigest(),
                source_before_sha256=sha256(source),source_sha256=sha256(destination),source_bytes=destination.stat().st_size)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--destination',type=Path,required=True)
    args=parser.parse_args();audit=curate(args.source,args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)
