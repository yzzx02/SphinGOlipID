"""Retain historical reference types separately from qualified MS2 evidence."""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sqlite3
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sphingolipid_toolkit._vendor.lipidgate.ms2.esterified_ceramide import is_esterified_ceramide
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import IndexedLibrary
from sphingolipid_toolkit._vendor.lipidgate.ms2.library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C, HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N, OXYGEN_MONOISOTOPIC_MASS as O,
    PHOSPHORUS_MONOISOTOPIC_MASS as P, PROTON_MONOISOTOPIC_MASS as PROTON,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint, sha256
try:
    from scripts.extend_positive_catalog_glycosphingolipid_library import read_catalog, NAME
except ModuleNotFoundError:
    from extend_positive_catalog_glycosphingolipid_library import read_catalog, NAME

ELEMENTS = {'C':C, 'H':H, 'N':N, 'O':O, 'P':P, 'S':31.9720711744}
MODIFIED = {'N-methylSo', 'ketoSo', '3-keto Cer',
            '1-O-acetyl-Cer', '1-O-acetyl-ω-linoleoyloxy-Cer'}


def reference_type_name(name):
    match=NAME.fullmatch(name)
    if match:return match['class']
    if name.startswith('M(IP)2C('):return 'M(IP)2C'
    # Positional prefixes such as (6-t...) describe a modified composition;
    # they must not create a separate type for every carbon/DB combination.
    return name.split('(',1)[0]


def reviewed_reference_formula(reference_type, formula):
    """Return reviewed composition and provenance without rewriting source rows."""
    if reference_type!='LM1':return formula,''
    tokens=re.findall(r'([A-Z][a-z]?)(\d*)',formula)
    if ''.join(a+b for a,b in tokens)!=formula or any(a not in 'CHNO' for a,b in tokens):
        raise ValueError('Review unexpected LM1 reference formula before correction')
    counts={a:int(b or 1) for a,b in tokens}
    if counts.get('N')!=2:
        raise ValueError('LM1 correction requires the reviewed two-nitrogen legacy formula')
    for element,count in zip('CHNO',(11,17,1,8)):
        counts[element]=counts.get(element,0)+count
    corrected=''.join(element+(str(counts[element]) if counts[element]!=1 else '')
        for element in 'CHNO' if counts.get(element))
    return corrected,'LM1_add_missing_NeuAc_residue;doi:10.3390/ijms22168844'


def build(catalog, index, output, audit, aliases_output, types_output):
    catalog, index, output, audit, aliases_output, types_output = map(
        Path, (catalog, index, output, audit, aliases_output, types_output))
    config = json.loads((ROOT/'config/sphingolipid_reference_class_aliases.json').read_text(encoding='utf-8'))
    if sha256(catalog) != config['source_catalog_sha256']:
        raise ValueError('Review aliases before using a different historical catalog')
    aliases = config['aliases']
    masses = defaultdict(list)
    with sqlite3.connect(f'file:{index.as_posix()}?mode=ro', uri=True) as conn:
        meta = json.loads(conn.execute("SELECT value FROM metadata WHERE key='library'").fetchone()[0])
        if meta['rules_sha256'] != code_fingerprint(ROOT/'src/sphingolipid_toolkit/_vendor/lipidgate/ms2'):
            raise ValueError('Rebuild the library index using current rules first')
        adduct = meta['adducts'].index('[M+H]+')
        for cls, value in conn.execute('SELECT class_id,precursor FROM entries WHERE adduct_id=? ORDER BY precursor', (adduct,)):
            masses[meta['classes'][cls]].append(value)
    # Do not count an ordinary Cer with an isobaric precursor as evidence for
    # an esterified subtype. The subtype itself must be present in the library.
    ester_masses=[]
    library=IndexedLibrary(index, allowed_class_keys={'CER'})
    try:
        for record in library:
            if record.adduct=='[M+H]+' and is_esterified_ceramide('Cer',record.lipid_chain_name):
                ester_masses.append(record.precursor_mz)
    finally:
        library.close()
    ester_masses.sort()
    rows = read_catalog(catalog)
    records, seen = [], set()
    counts, by_type = Counter(), defaultdict(Counter)
    aliases_by_class = defaultdict(set)
    canonical_by_type = {}
    for row in rows:
        original_formula, name = row[0], row[3]
        formula,correction=reviewed_reference_formula(reference_type_name(name),original_formula)
        if (name,formula) in seen: continue
        seen.add((name,formula))
        head=reference_type_name(name)
        cls=aliases.get(head,'')
        canonical_by_type[head]=cls
        if cls: aliases_by_class[cls].add(head)
        tokens=re.findall(r'([A-Z][a-z]?)(\d*)',formula)
        valid=''.join(a+b for a,b in tokens)==formula and all(a in ELEMENTS for a,b in tokens)
        value=sum(ELEMENTS[a]*int(b or 1) for a,b in tokens)+PROTON if valid else None
        status='reference_only_unvalidated_class'
        if not valid:status='invalid_reference_formula'
        elif name.startswith('myriocin'):status='not_a_sphingolipid'
        elif head=='N,N,N-trimethylSo':
            status='reference_only_insufficient_structural_signature'
            value=value-PROTON-0.000548579909065 if value is not None else None
            correction='intrinsic_monocation_[M]+;doi:10.1371/journal.pone.0144817'
        elif head=='LM1':status='reference_only_verified_composition_unresolved_linkage'
        elif head in MODIFIED:status='reference_only_ambiguous_modification'
        elif cls:
            values=ester_masses if head in {'ω-linoleoyloxy-Cer','ω-linoleoyloxy-DHCer'} else masses[cls]
            left=bisect.bisect_left(values,value*(1-10e-6))
            covered=left<len(values) and abs(values[left]-value)/value<=10e-6
            status='ms2_composition_covered' if covered else 'outside_current_eligible_library_grid'
        counts[status]+=1
        by_type[head][status]+=1
        records.append(dict(reference_name=name,formula=formula,precursor_mz=f'{value:.8f}' if value else '',
            adduct=('[M]+' if head=='N,N,N-trimethylSo' else '[M+H]+') if value else '',
            canonical_class=cls,coverage_status=status,original_formula=original_formula if correction else '',
            reference_correction=correction))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as binary,gzip.GzipFile(filename='',fileobj=binary,mode='wb',mtime=0,compresslevel=9) as compressed,io.TextIOWrapper(compressed,encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    report=dict(source_sha256=sha256(catalog),input_rows=len(rows),unique_name_formula_rows=len(records),
        unique_names=len({r['reference_name'] for r in records}),reference_types=len(by_type),
        coverage_status=dict(counts),by_reference_class={k:dict(v) for k,v in sorted(by_type.items())},
        ms2_coverage_library_source_sha256=meta['source_sha256'],rules_sha256=meta['rules_sha256'],
        asset_sha256=sha256(output),asset_bytes=output.stat().st_size,
        coverage_definition='Class/subtype and precursor composition within 10 ppm; theoretical availability, not measured recall or unique structure confirmation.')
    audit.parent.mkdir(parents=True,exist_ok=True)
    audit.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    aliases_output.parent.mkdir(parents=True,exist_ok=True)
    aliases_output.write_text(json.dumps({k:sorted(v) for k,v in sorted(aliases_by_class.items())},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    types_output.parent.mkdir(parents=True,exist_ok=True)
    with types_output.open('w',encoding='utf-8-sig',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['reference_type','canonical_class','ms2_type_available','reference_identities','composition_covered','status'])
        writer.writeheader()
        for name,values in sorted(by_type.items()):
            cls=canonical_by_type[name]
            writer.writerow(dict(reference_type=name,canonical_class=cls,ms2_type_available=bool(masses.get(cls)),
                reference_identities=sum(values.values()),composition_covered=values['ms2_composition_covered'],
                status='; '.join(sorted(values))))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('catalog','index','output','audit','aliases-output','types-output'):
        parser.add_argument('--'+arg,type=Path,required=True)
    args=parser.parse_args()
    result=build(args.catalog,args.index,args.output,args.audit,args.aliases_output,args.types_output)
    print(json.dumps({k:v for k,v in result.items() if k!='by_reference_class'},indent=2))
