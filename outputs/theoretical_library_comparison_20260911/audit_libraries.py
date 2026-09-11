"""Read theoretical libraries only; never imports the identification pipeline."""
import csv
import hashlib
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook

OUT = Path(__file__).resolve().parent
PRIMARY = Path(r'C:\Users\Administrator\xwechat_files\wxid_uig3uqawfvzj12_127a\msg\file\2026-06\常规鞘脂一级库3.0.csv')


def convolution(a, b):
    result = Counter()
    for (c, u), n in a.items():
        for (d, v), m in b.items():
            result[c + d, u + v] += n * m
    return result


def primary():
    original = PRIMARY.read_bytes()
    rows = [r for r in csv.reader(io.StringIO(original.decode('gb18030')))
            if r and re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', r[0])]
    names = {r[3] for r in rows}
    pattern = re.compile(r'^(.*)\(([^()]*(?:m|d|t|q))(\d+):(\d+)\)(.*)$')
    single_heads = {'So', 'ketoSo', 'Sa', 'ketoSa', 'N-methylSo', 'N,N-dimethylSo',
                    'N,N,N-trimethylSo', 'So1P', 'Sa1P', 'Lyso-SM', 'Lyso-sulfo',
                    'GlcSo', 'GalSo', 'Gb3So'}
    lcb = Counter((c, u) for c in range(12, 31) for u in range(6))
    fa = Counter((c, u) for c in range(2, 43) for u in range(4))
    pair = convolution(lcb, fa)
    triple = convolution(pair, fa)
    assert sum(lcb.values()) == 114
    assert sum(fa.values()) == 164
    assert sum(pair.values()) == 18696
    assert sum(triple.values()) == 3066144
    # Independent analytical check of the two-chain distribution.
    for (c, u), n in pair.items():
        n_c = max(0, min(30, c-2) - max(12, c-42) + 1)
        n_u = max(0, min(5, u) - max(0, u-3) + 1)
        assert n == n_c * n_u
    totals, parents, templates = Counter(), Counter(), defaultdict(set)
    for name in names:
        m = pattern.fullmatch(name)
        if not m:
            assert name == 'myriocin', name
            totals['fixed_named'] += 1
            continue
        head, prefix, c, u, suffix = m.groups()
        c, u = int(c), int(u)
        kind = ('single' if head in single_heads else 'three_variable_chains'
                if head.startswith('1-O-acetyl') else 'two_variable_chains')
        # The user explicitly authorized two independently attached FAs with
        # the same ranges for the 1-O-acyl families. Fixed linoleoyl groups
        # outside the name's variable-chain total are not enumerated again.
        pool = {'single': lcb, 'two_variable_chains': pair,
                'three_variable_chains': triple}[kind]
        totals[kind] += pool[c, u]
        parents[kind] += 1
        templates[kind].add((head, prefix, suffix))
    assert len(set((r[0], r[3]) for r in rows)) == len(names)
    result = dict(source_file=PRIMARY.name, sha256=hashlib.sha256(original).hexdigest(),
                  manuscript_ms1=84240, file_records=len(rows), unique_names=len(names),
                  duplicate_records=len(rows)-len(names), expansion=sum(totals.values()),
                  expansion_by_chain_type=dict(totals), unique_parents_by_type=dict(parents),
                  templates_by_type={k: len(v) for k,v in templates.items()},
                  lcb_c=[12,30], lcb_db=[0,5], fa_c=[2,42], fa_db=[0,3],
                  chemically_filtered=False)
    return result


def lda():
    all_keys, counts = set(), {}
    for p in sorted((OUT/'sources').glob('LDA2*noSphingoOH.xlsx')):
        book = load_workbook(p, read_only=True, data_only=True)
        keys = set()
        for head in ['SM', 'Cer']:
            for row in list(book[head].values)[2:]:
                if (isinstance(row[0], (int,float)) and isinstance(row[2], (int,float))
                        and all(isinstance(x, (int,float)) for x in row[3:8])):
                    # C,H,O,P,N are neutral composition columns; exclude mass/adducts.
                    keys.add((head, tuple(row[3:8])))
        counts[p.name] = len(keys)
        all_keys.update(keys)
        book.close()
    return dict(unique_species_candidates=len(all_keys), per_file=counts,
                scope='Official recommended example mass lists; not software capacity')


def msdial():
    paths = [Path('F:/MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp'),
             OUT/'sources'/'MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp']
    expected_md5 = ['f1a28bfdadaadd5f47c1ecf04c0e88a6', '0b90719cc5f4f84ed7a656ec2c5f1fa3']
    all_keys, per_file, classes = set(), {}, Counter()
    for p, expected in zip(paths, expected_md5):
        with p.open('rb') as f:
            assert hashlib.file_digest(f, 'md5').hexdigest() == expected, f'Incomplete/unverified file: {p}'
        keys, records, rec = set(), 0, {}
        def flush():
            nonlocal records
            h = rec.get('COMPOUNDCLASS', '')
            sp = ('Cer' in h or h in {'SM','ASM','SL','SPB','Sphingosine','Sphinganine',
                  'Phytosphingosine','Sph','DHSph','PhytoSph','S1P','SPBP',
                  'GM1','GM2','GM3','GD1a','GD1b','GD2','GD3','GT1b','GQ1b'})
            if sp:
                assert rec.get('NAME') and rec.get('FORMULA')
                keys.add((h, rec['NAME'], rec['FORMULA']))
                classes[h] += 1
                records += 1
        with p.open(encoding='utf-8-sig') as f:
            for line in f:
                line=line.strip()
                if not line:
                    flush();rec={}
                elif ': ' in line:
                    k,v=line.split(': ',1)
                    if k in {'NAME','FORMULA','COMPOUNDCLASS'}:rec[k]=v
            if rec:flush()
        all_keys.update(keys)
        per_file[p.name] = dict(records=records, unique_candidates=len(keys))
    return dict(unique_candidates=len(all_keys), per_file=per_file,
                included_class_records=dict(classes), version='VS69')


def lipidin():
    all_keys, per_file = set(), {}
    for mode in ['pos', 'neg']:
        p=OUT/f'lipidin_{mode}_structures.tsv'
        with p.open(encoding='utf-8-sig') as f:
            keys={(r['structure_key'],r['Formula']) for r in csv.DictReader(f,delimiter='\t')}
        per_file[mode]=len(keys)
        all_keys.update(keys)
    return dict(unique_candidates=len(all_keys), per_file=per_file,
                publication_SP_hierarchical_records=6252820,
                version='Local LipidIN-main.zip pos_ALL.rda + neg_ALL.rda')


if __name__ == '__main__':
    data={'SphinGOlipID':primary(),'LDA2':lda()}
    (OUT/'primary_audit.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
