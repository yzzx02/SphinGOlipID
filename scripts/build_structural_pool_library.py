"""Build the positive reference from the preceding checked MSP, retaining stable records and order."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re

from extend_verified_modified_sphingoid_library import rebuild_verified_index, sha256
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import POSITIVE_CHOLINE_COMMON_FRAGMENTS, normalize_imported_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_phosphoceramide import HEADS, fragments, C, H, N, O, P, PROTON
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_complex_ceramide import CeramideCore

EVIDENCE = Path(__file__).resolve().parents[1] / 'config/positive_structural_chain_evidence_20261126.json'

def chain_block(entry, revision):
    name = entry['name']
    matched = re.fullmatch(r'([^()]+)\(([dt])(\d+):(\d+)/(h?)(\d+):(\d+)\)', name)
    if not matched:
        raise ValueError(f'Expected an explicit curated chain pair: {name}')
    cls, series, lc, db, hyd, fa, fdb = matched.groups()
    lc, db, fa, fdb = map(int, (lc, db, fa, fdb))
    if fa < 14 or db > 3 or fdb > 3 or entry['adduct'] != '[M+H]+':
        raise ValueError(name)
    core = CeramideCore(series, lc, db, bool(hyd), fa, fdb)
    cc, ch, cn, co = core.counts
    hc, hh, hn, ho, hp = HEADS[cls]
    counts = (cc + hc, ch + hh - 2, cn + hn, co + ho - 1, hp)
    formula = ''.join(e + (str(n) if n > 1 else '') for e, n in zip('CHNOP', counts) if n)
    precursor = sum(x * m for x, m in zip(counts, (C, H, N, O, P))) + PROTON
    if entry['formula'] != formula or abs(entry['observed_precursor_mz'] / precursor - 1) > 1e-5:
        raise ValueError(f'Formula or precursor mismatch for {name}')
    ions = [f for f in core.lcb_ions if f.name in entry['library_lcb_channels']]
    measured = {ion['name']: ion for ion in entry['measured_lcb_products']}
    if not ions or len(ions) != len(measured):
        raise ValueError(f'Missing evidence for {name}')
    for f in ions:
        ion = measured[f.name]
        if abs(ion['mz'] / f.mz - 1) > 2e-5 or ion.get('fraction_of_base_peak', 1.) < .01:
            raise ValueError(f'LCB mass or intensity outside the reference criteria: {name}')
    symbols = {2: 'd', 3: 't'}
    total = f'{cls}({symbols[(2 if series == "d" else 3) + bool(hyd)]}{lc+fa}:{db+fdb})'
    products = sorted([*fragments(cls, precursor), *ions], key=lambda x: x.mz)
    origin = entry.get('publication') or f"{entry['acquisition']} scan {entry['scan_id']}"
    lines = [f'Name: {name}', f'PrecursorMZ: {precursor:.8f}', 'PrecursorType: [M+H]+',
             f'CompoundClass: {cls}', f'Formula: {formula}', f'Comment: MS1_name={total};polarity=+',
             f'LibraryRevision: {revision}', f'Source: {origin}', 'SpectrumOrigin: theoretical products; recorded chain evidence',
             f'Num Peaks: {len(products)}']
    lines.extend(f'{f.mz:.4f} 100.00 "{f.name}" "{f.fragment_type}"' for f in products)
    return '\n'.join(lines).encode('utf8'), (cls, name, formula, '[M+H]+')

def update_lyso(block):
    header, _, peakpart = block.partition(b'Num Peaks: ')
    count, _, data = peakpart.partition(b'\n')
    lines = data.splitlines()
    if len(lines) != int(count):
        raise ValueError('Invalid number of Lyso-SM peaks')
    for mz, name in POSITIVE_CHOLINE_COMMON_FRAGMENTS:
        if any(f'"{name}"'.encode() in line for line in lines):
            raise ValueError('Use the preceding Lyso-SM library as input')
        role = 'Diagnostic_HG' if name == '[C2H6O4P]+' else 'Common'
        lines.append(f'{mz:.4f} 100.00 "{name}" "{role}"'.encode())
    return header + b'Num Peaks: ' + str(len(lines)).encode() + b'\n' + b'\n'.join(lines)

def expected(record):
    if record.compound_class == 'LSM' and record.adduct == '[M+H]+':
        return normalize_imported_record(record)
    return record

def build(source, previous_index, directory):
    directory.mkdir(exist_ok=True, parents=True)
    config = json.loads(EVIDENCE.read_text(encoding='utf8'))
    revision = config['revision']
    output = directory / 'sphingolipid_positive.msp.gz'
    content = hashlib.sha256()
    additions = set()
    statistics = Counter()
    unchanged = hashlib.sha256()
    residual = b''
    with gzip.open(source, 'rb') as reader, output.open('xb') as raw_out, gzip.GzipFile(filename='', fileobj=raw_out, mode='wb', mtime=0, compresslevel=5) as writer:
        def emit(block):
            data = block + b'\n\n'
            writer.write(data)
            content.update(data)
        for chunk in iter(lambda: reader.read(4 << 20), b''):
            blocks = (residual + chunk).split(b'\n\n')
            residual = blocks.pop()
            for block in blocks:
                if b'\nCompoundClass: LSM\n' in block and b'\nPrecursorType: [M+H]+\n' in block:
                    block = update_lyso(block)
                    statistics['Lyso-SM modified'] += 1
                else:
                    unchanged.update(block + b'\n\n')
                    statistics['unchanged'] += 1
                emit(block)
        if residual.strip() or statistics['Lyso-SM modified'] != 51:
            raise ValueError('Unexpected source coverage or truncated MSP')
        for entry in config['entries']:
            block, identity = chain_block(entry, revision)
            if identity in additions:raise ValueError(identity)
            additions.add(identity)
            emit(block)
    print('MSP written', statistics, 'new explicit chains', len(additions), flush=True)
    out = rebuild_verified_index(output, previous_index, directory / 'sphingolipid_positive.sqlite',
                                 revision=revision, additional_identities=additions, expected_transform=expected)
    out.update(raw_content_sha256=content.hexdigest(), untouched_blocks_sha256=unchanged.hexdigest(),
               previous_MSP_sha256=sha256(source), block_counts=dict(statistics), added_names=sorted(x[1] for x in additions))
    (directory/'positive.audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n', encoding='utf8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','previous_index','destination'):parser.add_argument(name,type=Path)
    args=parser.parse_args()
    build(args.source,args.previous_index,args.destination)

