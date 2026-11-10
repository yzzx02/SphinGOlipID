"""Add positive ganglioside family/adduct coverage to a separate branch MSP."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import replace
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import PEAK_RE,header_value,iter_blocks,parse_peaks,validate_block
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_ganglioside import (
    ADDUCTS,LEGACY_CLASSES as CLASSES,GangliosideComposition,from_name,normalize_fragments,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import fragment_type,DIHEXOSYL_CLASSES
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import chain_double_bond_violations,explicit_acyl_carbons
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

CHAIN=re.compile(r'^(Cer|GM3|GM1)\(([dt])(\d+):(\d+)/(h?)(\d+):(\d+)\)$')
REVISION='ganglioside_coverage_20261108'

def chain_tuple(name):
    m=CHAIN.fullmatch(name)
    if m is None:return None
    return m[2],int(m[3]),int(m[4]),bool(m[5]),int(m[6]),int(m[7])

def eligible(name,cls):
    return not chain_double_bond_violations(name,cls) and all(c>=14 for c in explicit_acyl_carbons(name))

def collect_models(source):
    gm3,cer,legacy,existing=set(),set(),set(),set()
    for block in iter_blocks(source):
        cls,adduct,name=(header_value(block,k) for k in ('CompoundClass','PrecursorType','Name'))
        if cls in CLASSES and adduct in ADDUCTS:existing.add((cls,name,adduct))
        parsed=chain_tuple(name)
        if parsed is None or not eligible(name,cls):continue
        if cls=='GM3' and adduct=='[M+H]+':gm3.add(parsed)
        if cls=='Cer' and adduct=='[M+H]+':cer.add(parsed)
        if cls=='GM1' and adduct=='[M+2NH4]2+':legacy.add(parsed)
    if not gm3:raise ValueError('No existing protonated GM3 chain grid to extend')
    common=set(gm3)
    for chain in cer:
        series,lc,ld,hydroxy,fc,fd=chain
        long_acyl=series=='d' and lc==18 and ld in (1,2) and 26<=fc<=30
        if series=='t' or chain in gm3 or long_acyl:common.add(chain)
        # This is the existing ordinary Cer d/h expansion, not an inferred
        # chain from a sum composition and not a benchmark-derived addition.
        if series=='d' and not hydroxy and (chain in gm3 or long_acyl):
            common.add((series,lc,ld,True,fc,fd))
    models={cls:tuple(GangliosideComposition(cls,*chain) for chain in sorted(common|(legacy if cls=='GM1' else set())))
            for cls in sorted(CLASSES)}
    return models,existing,dict(gm3_source_chains=len(gm3),ordinary_cer_source_chains=len(cer),
        shared_seed_chains=len(common),additional_existing_gm1_chains=len(legacy-common))

def curated_block(block):
    cls,adduct,name=(header_value(block,k) for k in ('CompoundClass','PrecursorType','Name'))
    old=parse_peaks(block)
    fragments=[FragmentRecord(p['mz'],p['name'],p['type'],p['intensity']) for p in old]
    if cls in CLASSES and adduct in ADDUCTS:
        revised=normalize_fragments(cls,name,adduct,fragments)
    elif cls in DIHEXOSYL_CLASSES and adduct=='[M+H]+':
        revised=[replace(f,fragment_type=fragment_type(cls,name,f)) for f in fragments]
    else:return block
    assert [(f.name,f.mz,f.intensity) for f in revised[:len(old)]]==[(p['name'],p['mz'],p['intensity']) for p in old]
    result=[];position=0
    for line in block:
        match=PEAK_RE.match(line)
        if match is not None:
            f=revised[position];position+=1
            line=line[:match.start('type')]+f.fragment_type+line[match.end('type'):]
        elif line.startswith('Num Peaks:'):line=f'Num Peaks: {len(revised)}'
        result.append(line)
    result += [f'{f.mz:.4f} {f.intensity:.2f} "{f.name}" "{f.fragment_type}"' for f in revised[len(old):]]
    validate_block(result)
    return result

def model_block(model,adduct):
    fragments=model.fragments(adduct)
    block=[f'Name: {model.name}',f'PrecursorMZ: {model.precursor_mz(adduct):.8f}',
        f'PrecursorType: {adduct}',f'CompoundClass: {model.compound_class}',f'Formula: {model.formula}',
        f'Comment: MS1_name={model.species_name};polarity=+',
        'Source: theoretical glycosidic-cleavage family model',f'LibraryRevision: {REVISION}',
        'AnnotationScope: glycan_composition_family;no_linkage_isomer_assignment',f'Num Peaks: {len(fragments)}']
    block += [f'{f.mz:.4f} {f.intensity:.2f} "{f.name}" "{f.fragment_type}"' for f in fragments]
    validate_block(block)
    return block

def curate(source,destination):
    source,destination=Path(source),Path(destination)
    root=Path(__file__).resolve().parents[1]
    protected={root/'libraries/ms2'/name for name in ('current_positive.msp.gz','current_negative.msp.gz')}
    if source.resolve()==destination.resolve() or destination.resolve() in {p.resolve() for p in protected}:
        raise ValueError('Write to a separate destination; never overwrite original libraries')
    before=sha256(source)
    models,existing,seeds=collect_models(source)
    additions=[(model,adduct) for cls,group in models.items() for model in group for adduct in ADDUCTS
               if (cls,model.name,adduct) not in existing]
    changes=Counter();added=Counter();prefix=hashlib.sha256();source_count=0
    destination.parent.mkdir(parents=True,exist_ok=True)
    with (destination.open('wb') as binary,gzip.GzipFile(filename='',fileobj=binary,mode='wb',compresslevel=9,mtime=0) as compressed,
          io.TextIOWrapper(compressed,encoding='utf-8',newline='\n') as out):
        for block in iter_blocks(source):
            revised=curated_block(block);source_count+=1
            if revised!=block:changes[header_value(block,'CompoundClass')+' '+header_value(block,'PrecursorType')]+=1
            prefix.update(('\n'.join(revised)+'\n\n').encode('utf-8'))
            out.write('\n'.join(revised)+'\n\n')
        for model,adduct in additions:
            if not eligible(model.name,model.compound_class):raise ValueError('Seed outside branch chain limits: '+model.name)
            out.write('\n'.join(model_block(model,adduct))+'\n\n')
            added[model.compound_class+' '+adduct]+=1
    revised=iter(iter_blocks(destination))
    for block in iter_blocks(source):
        if next(revised,None)!=curated_block(block):raise ValueError('Existing identity, mass, intensity, fragment or order changed unexpectedly')
    for model,adduct in additions:
        if next(revised,None)!=model_block(model,adduct):raise ValueError('Unexpected appended theoretical record')
    if next(revised,None) is not None or sha256(source)!=before:raise ValueError('Unexpected trailing record or source mutation')
    return dict(revision=REVISION,previous_source_records=source_count,source_records=source_count+len(additions),
        changed_existing_by_class_adduct=dict(changes),added_records=len(additions),added_by_class_adduct=dict(added),
        seed_policy=seeds,source_before_sha256=before,source_sha256=sha256(destination),source_bytes=destination.stat().st_size,
        original_prefix_record_ids_names_precursors_masses_intensities_and_order_preserved=True,
        curated_prefix_sha256=prefix.hexdigest(),prior_legacy_gd1a_gd1b_gt1b_records_unchanged=True,
        new_gd1_gt1_profiles_are_linkage_unresolved_families=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--destination',type=Path,required=True)
    args=parser.parse_args();result=curate(args.source,args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
