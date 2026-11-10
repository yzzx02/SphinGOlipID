"""Append formula-checked glycan totals backed by eligible existing Cer parents."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import _encode_fragment_fields
from sphingolipid_toolkit._vendor.lipidgate.ms2.library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C, HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N, OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,
    PHOSPHORUS_MONOISOTOPIC_MASS as P,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_glycosphingolipid import PROFILES, RESIDUE_FORMULAS, fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_phosphoceramide import HEADS as PHOSPHO_HEADS, ALIASES as PHOSPHO_ALIASES, fragments as phosphoceramide_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_legacy_sphingolipid import (HEADS as LEGACY_HEADS,ALIASES as LEGACY_ALIASES,
    FREE_BASE_CLASSES,fragments as legacy_fragments,S)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import chain_double_bond_violations, explicit_acyl_carbons
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

NAME = re.compile(r'^(?P<class>.*)\((?P<series>[mdt])(?P<c>\d+):(?P<db>\d+)\)(?:\((?P<oh>\d*)OH\))?$')
PARENT = re.compile(r'^Cer\((?P<series>[mdt])(?P<lcb_c>\d+):(?P<lcb_db>\d+)/(?P<hydroxy>h?)(?P<fa_c>\d+):(?P<fa_db>\d+)\)(?:\((?P<oh>\d*)OH\))?$')
FREE_PARENT=re.compile(r'^(?:SPB|Sph|DHSph|PhytoSph)\(([mdt])(\d+):(\d+)\)$')


def formula_atoms(value):
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', str(value))
    if ''.join(a+b for a,b in tokens) != value or set(a for a,_ in tokens)-set('CHNO'):
        return None
    values = {a:int(b or 1) for a,b in tokens}
    return tuple(values.get(a,0) for a in 'CHNO')


def read_catalog(path):
    data = Path(path).read_bytes()
    for encoding in ('utf-8-sig','gb18030'):
        try:
            text = data.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError('Catalog is neither UTF-8 nor GB18030')
    rows = list(csv.reader(text.splitlines()))
    start = next(i for i,row in enumerate(rows) if row and row[0].startswith('# Formula'))
    return [row for row in rows[start+1:] if len(row)>3]


def catalog_candidates(rows):
    """Derive only unique residue formulas; preserve uncertain names separately."""
    masses = dict(zip(PROFILES, (tuple(sum(n*r[i] for n,r in zip(comp,RESIDUE_FORMULAS)) for i in range(4))
                                for comp in PROFILES.values())))
    by_atoms = {}
    aliases = json.loads((Path(__file__).resolve().parents[1] / 'config/positive_catalog_glycosphingolipid_aliases.json').read_text(encoding='utf-8'))['aliases']
    for cls,atoms in masses.items():
        by_atoms.setdefault(atoms,[]).append(cls)
    for row in rows:
        match = NAME.fullmatch(row[3])
        atoms = formula_atoms(row[0])
        if not match or atoms is None:
            continue
        expected_class = aliases.get(match['class'])
        if expected_class is None:
            continue
        carbon,db = int(match['c']),int(match['db'])
        head_atoms=masses[expected_class]
        core=tuple(a-b for a,b in zip(atoms,head_atoms))
        # Historical (2OH) may denote a hydroxyl position. Total oxygen is
        # established from the formula and eligible parent, never that suffix.
        if core[:3] == (carbon,2*carbon+1-2*db,1):
            yield expected_class,core,row,match


def phosphoceramide_candidates(rows):
    for row in rows:
        match=NAME.fullmatch(row[3])
        if match is None or match['class'] not in PHOSPHO_ALIASES:continue
        cls=PHOSPHO_ALIASES[match['class']]
        tokens=re.findall(r'([A-Z][a-z]?)(\d*)',row[0])
        if ''.join(a+b for a,b in tokens)!=row[0] or set(a for a,b in tokens)-set('CHNOP'):continue
        counts={a:int(b or 1) for a,b in tokens}
        atoms=tuple(counts.get(a,0) for a in 'CHNOP')
        # Recover the attached head residue, never an acyl chain from a sum.
        residue=list(PHOSPHO_HEADS[cls]);residue[1]-=2;residue[3]-=1
        core5=tuple(a-b for a,b in zip(atoms,residue))
        if core5[-1]!=0:continue
        carbon,db=int(match['c']),int(match['db'])
        if core5[:3]==(carbon,2*carbon+1-2*db,1):yield cls,core5[:4],row,match


def legacy_candidates(rows):
    for row in rows:
        match=NAME.fullmatch(row[3])
        if match is None or match['class'] not in LEGACY_ALIASES:continue
        cls=LEGACY_ALIASES[match['class']]
        if cls=='KetoSPB' and (match['series'],int(match['c']),int(match['db']))!=('d',18,0):
            continue  # Primary-standard 3KDS coverage; no unverified keto grid.
        tokens=re.findall(r'([A-Z][a-z]?)(\d*)',row[0])
        if ''.join(a+b for a,b in tokens)!=row[0] or set(a for a,b in tokens)-set('CHNOPS'):continue
        values={a:int(b or 1) for a,b in tokens}
        atoms=tuple(values.get(a,0) for a in 'CHNOPS')
        core6=tuple(a-b for a,b in zip(atoms,LEGACY_HEADS[cls]))
        carbon,db=int(match['c']),int(match['db'])
        hydrogens=2*carbon+(3 if cls in FREE_BASE_CLASSES else 1)-2*db
        if core6[:3]==(carbon,hydrogens,1) and core6[4:]==(0,0):
            yield cls,core6[:4],row,match


def extend(source, destination, catalog):
    source,destination,catalog = map(Path,(source,destination,catalog))
    originals = {source.parent/name for name in ('current_positive.msp.gz','current_negative.msp.gz')}
    if source.resolve()==destination.resolve() or destination.resolve() in {p.resolve() for p in originals}:
        raise ValueError('Use a separate destination; preserve original libraries')
    before = sha256(source)
    parents = {}
    phyto_parents={}
    free_parents={}
    existing = set()
    source_count = 0
    for block in iter_blocks(source):
        source_count += 1
        cls,adduct,name = (header_value(block,key) for key in ('CompoundClass','PrecursorType','Name'))
        if cls in PROFILES or cls in PHOSPHO_HEADS or cls in LEGACY_HEADS:
            existing.add((cls,header_value(block,'Formula')))
        free=FREE_PARENT.fullmatch(name) if cls=='SPB' and adduct=='[M+H]+' else None
        if free is not None and not chain_double_bond_violations(name,'SPB'):
            series,carbon,db=free.groups();carbon,db=int(carbon),int(db)
            atoms=(carbon,2*carbon+3-2*db,1,{'m':1,'d':2,'t':3}[series])
            expected=sum(n*m for n,m in zip(atoms,(C,H,N,O)))+PROTON
            if abs(float(header_value(block,'PrecursorMZ'))-expected)<.01:free_parents.setdefault(atoms,name)
        parent = PARENT.fullmatch(name)
        if cls!='Cer' or adduct!='[M+H]+' or parent is None:
            continue
        acyls = explicit_acyl_carbons(name)
        if not acyls or min(acyls)<14 or chain_double_bond_violations(name,'Cer'):
            continue
        carbon = int(parent['lcb_c'])+int(parent['fa_c'])
        db = int(parent['lcb_db'])+int(parent['fa_db'])
        hydroxyls = (1 if parent['hydroxy'] else 0)+(int(parent['oh'] or 1) if parent['oh'] is not None else 0)
        atoms = (carbon,2*carbon+1-2*db,1,{'m':2,'d':3,'t':4}[parent['series']]+hydroxyls)
        stored_formula = header_value(block,'Formula')
        if stored_formula and formula_atoms(stored_formula)!=atoms:
            continue
        expected = sum(n*m for n,m in zip(atoms,(C,H,N,O)))+PROTON
        if abs(float(header_value(block,'PrecursorMZ'))-expected)>.01:
            continue
        parents.setdefault(atoms,name)
        if parent['series']=='t':
            phyto_parents.setdefault(atoms,name)
    rows = read_catalog(catalog)
    additions = {}
    aliases = {}
    from itertools import chain
    for cls,core,row,match in chain(catalog_candidates(rows),phosphoceramide_candidates(rows),legacy_candidates(rows)):
        eligible=free_parents if cls in FREE_BASE_CLASSES else parents
        if core not in eligible:
            continue
        aliases.setdefault(cls,set()).add(match['class'])
        # O is molecular total oxygen, not a location of a hydroxyl group.
        oxygen=core[3]
        series=({1:'m',2:'d'} if cls in FREE_BASE_CLASSES else {2:'m',3:'d'}).get(oxygen,'t')
        total = f"{cls}({series}{core[0]}:{match['db']})"
        if oxygen>4:
            total += '(OH)' if oxygen==5 else f'({oxygen-4}OH)'
        key=(cls,row[0])
        if key not in existing:additions.setdefault(key,(total,eligible[core],core))
    # The old catalog lists only hydroxylated M(IP)2C. Preserve the family on
    # the CURRENT eligible phyto-Cer grid, rather than adding unsupported OH
    # chains just to reproduce historical compositions.
    if any(row[3].startswith('M(IP)2C(') for row in rows):
        residue=LEGACY_HEADS['MIP2C']
        for core,parent in sorted(phyto_parents.items()):
            values=tuple(a+b for a,b in zip(core+(0,0),residue))
            formula=''.join(a+(str(n) if n!=1 else '') for a,n in zip('CHNOPS',values) if n)
            key=('MIP2C',formula)
            db=(2*core[0]+1-core[1])//2
            name=f'MIP2C(t{core[0]}:{db})'
            if core[3]>4:name+=f'({core[3]-4}OH)'
            if key not in existing:additions.setdefault(key,(name,parent,core))
        if phyto_parents:aliases.setdefault('MIP2C',set()).add('M(IP)2C')
    destination.parent.mkdir(parents=True,exist_ok=True)
    digest = hashlib.sha256()
    with destination.open('wb') as binary,gzip.GzipFile(fileobj=binary,filename='',mode='wb',mtime=0,compresslevel=6) as compressed,io.TextIOWrapper(compressed,encoding='utf-8',newline='\n') as handle:
        for block in iter_blocks(source):
            text = '\n'.join(block)+'\n\n'
            handle.write(text)
            digest.update(text.encode('utf-8'))
        for key,(name,parent,core) in sorted(additions.items()):
            cls,formula=key[:2]
            parts={a:int(b or 1) for a,b in re.findall(r'([A-Z][a-z]?)(\d*)',formula)}
            mass=sum(parts.get(a,0)*m for a,m in zip('CHNOPS',(C,H,N,O,P,S)))+PROTON
            products = (fragments(cls,mass) if cls in PROFILES else
                phosphoceramide_fragments(cls,mass) if cls in PHOSPHO_HEADS else legacy_fragments(cls,mass,core[3]))
            handle.write(f'Name: {name}\nPrecursorMZ: {mass:.8f}\nPrecursorType: [M+H]+\nCompoundClass: {cls}\nFormula: {formula}\n')
            handle.write(f'Comment: MS1_name={name};polarity=+\n')
            handle.write(f'Curation: catalog_glycan_totals_20261110\nEligibleParent: {parent}\n')
            handle.write(f'Num Peaks: {len(products)}\n')
            for f in products:
                handle.write(f'{f.mz:.4f} {f.intensity:.2f} {_encode_fragment_fields(f)}\n')
            handle.write('\n')
    revised = iter(iter_blocks(destination))
    for old in iter_blocks(source):
        if next(revised,None)!=old:
            raise ValueError('An existing record or its order changed')
    if sum(1 for _ in revised)!=len(additions) or sha256(source)!=before:
        raise ValueError('Unexpected source mutation or addition count')
    return dict(source_sha256=before,catalog_sha256=sha256(catalog),catalog_rows=len(rows),
        original_records=source_count,added_records=len(additions),added_by_class=dict(Counter(key[0] for key in additions)),
        reference_aliases={k:sorted(v) for k,v in sorted(aliases.items())},
        eligible_core_compositions=len(parents),eligible_free_base_compositions=len(free_parents),preserved_records_and_order=True,
        retained_original_content_sha256=digest.hexdigest(),destination_sha256=sha256(destination))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--audit',type=Path,required=True)
    args = parser.parse_args()
    result = extend(args.source,args.destination,args.catalog)
    args.audit.parent.mkdir(parents=True,exist_ok=True)
    args.audit.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='reference_aliases'},indent=2),flush=True)
