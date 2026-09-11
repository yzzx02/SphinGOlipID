"""Combine existing library audits only; no experimental data inputs."""
import json
import csv
from audit_libraries import OUT, lipidin

data = json.loads((OUT / 'primary_audit.json').read_text(encoding='utf-8'))
data['MS-DIAL'] = json.loads((OUT / 'msdial_audit.json').read_text(encoding='utf-8'))
data['LipidIN'] = lipidin()
sets = {}
for mode in ['pos', 'neg']:
    with (OUT / f'lipidin_{mode}_structures.tsv').open(encoding='utf-8-sig') as f:
        sets[mode] = {(r['structure_key'], r['Formula']) for r in csv.DictReader(f, delimiter='\t')}
data['LipidIN']['polarity_overlap'] = len(sets['pos'] & sets['neg'])
data['LipiDetective'] = {'unique_candidates': None, 'scope': 'Not applicable: sequence prediction model'}
(OUT / 'library_totals.json').write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(data['LipidIN'], indent=2))

# Review exact-key disagreement across polarities without silently normalizing it.
names = {}
for mode, keys in sets.items():
    names[mode] = {}
    for name, formula in keys:
        names[mode].setdefault(name, set()).add(formula)
conflicts = {name: {m: sorted(names[m][name]) for m in names}
             for name in names['pos'].keys() & names['neg'].keys()
             if names['pos'][name] != names['neg'][name]}
(OUT / 'lipidin_formula_conflicts.json').write_text(json.dumps(conflicts, indent=2), encoding='utf-8')
print('Cross-polarity name/formula conflicts:', len(conflicts))
