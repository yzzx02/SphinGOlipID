"""Freshly normalize the positive triene revision or its matched negative index."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import re
import shutil

from extend_verified_modified_sphingoid_library import (
    ROOT, fragments, rebuild_verified_index, _encode_fragment_fields,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import (
    chain_double_bond_violations, explicit_chain_compositions, is_positive_sphingolipid,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

REVISION = 'sphingolipid-lcb-triene-20261115-r1'
NAME = 'SHexSPB(d22:3)'
IDENTITY = ('SHexSPB', NAME, 'C28H51NO10S', '[M+H]+')
PRE_LIMIT_SHA = 'b4b360b4016b84267113c54088e3e52d320e5dca1d253f49ac676e5c10d2cc1e'


def restore_curated_trienes(previous, output):
    """Recover explicit pre-limit references whose only exclusion was the C18 exception."""
    if sha256(previous) != PRE_LIMIT_SHA:
        raise ValueError('Use the archived positive source preceding the chain-limit curation')
    restored = set()
    with gzip.open(previous, 'rb') as src, output.open('ab') as dst, gzip.GzipFile(
        filename='', fileobj=dst, mode='wb', mtime=0
    ) as extra:
        remainder = b''
        for chunk in iter(lambda: src.read(4*1024*1024), b''):
            records = (remainder + chunk).split(b'\n\n')
            remainder = records.pop()
            for raw in records:
                first = raw.split(b'\n',1)[0]
                if not re.search(rb'[mdt]\d+:3(?!\d)', first, re.I):
                    continue
                headers = dict(line.split(': ',1) for line in raw.decode().splitlines() if ': ' in line)
                cls, name, adduct = (headers[k] for k in ('CompoundClass','Name','PrecursorType'))
                lcbs, acyls = explicit_chain_compositions(name, cls)
                if not is_positive_sphingolipid(cls, adduct) or not any(c!=18 and db==3 for c,db in lcbs):
                    continue
                if chain_double_bond_violations(name,cls) or any(c<14 for c,_ in acyls):
                    continue
                key = (cls, name, headers.get('Formula',''), adduct)
                if key in restored:
                    raise ValueError(f'Duplicate pre-limit candidate: {key}')
                restored.add(key)
                extra.write(raw+b'\n\n')
        if remainder.strip() or len(restored)!=1680:
            raise ValueError('Unexpected pre-limit library layout or restored count')
    return restored


def append_candidate(source, output, catalog):
    with gzip.open(catalog, 'rt', encoding='utf-8-sig') as handle:
        reference, = [r for r in csv.DictReader(handle) if r['reference_name'] == 'Lyso-sulfo(d22:3)']
    if (reference['canonical_class'], NAME, reference['formula'], reference['adduct']) != IDENTITY:
        raise ValueError('Unexpected archived candidate identity')
    if chain_double_bond_violations(NAME, 'SHexSPB'):
        raise ValueError('Candidate rejected by the current curation rule')
    mass = float(reference['precursor_mz'])
    products = fragments('SHexSPB', mass)
    lines = [f'Name: {NAME}', f'PrecursorMZ: {mass:.8f}', 'PrecursorType: [M+H]+',
             'CompoundClass: SHexSPB', f'Formula: {IDENTITY[2]}', f'Comment: MS1_name={NAME};polarity=+',
             f'LibraryRevision: {REVISION}', 'Source: Lyso-sulfo(d22:3); archived explicit single-base composition',
             'SpectrumOrigin: calculated class fragments; paper annotation is provisional', f'Num Peaks: {len(products)}']
    lines.extend(f'{p.mz:.4f} {p.intensity:.2f} {_encode_fragment_fields(p)}' for p in products)
    # A second gzip member preserves every original compressed and uncompressed byte.
    with source.open('rb') as old, output.open('xb') as out:
        shutil.copyfileobj(old, out)
        with gzip.GzipFile(filename='', fileobj=out, mode='wb', mtime=0) as extra:
            extra.write(('\n' + '\n'.join(lines) + '\n\n').encode('utf-8'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('polarity', choices=['positive', 'negative'])
    parser.add_argument('destination', type=Path)
    parser.add_argument('--baseline', type=Path, default=ROOT/'src/sphingolipid_toolkit/data')
    parser.add_argument('--pre-limit-source', type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    stem = 'sphingolipid_' + args.polarity
    source = args.baseline/'ms2'/(stem+'.msp.gz')
    allow = ()
    if args.polarity == 'positive':
        if args.pre_limit_source is None:
            parser.error('Positive builds require --pre-limit-source')
        output = args.destination/(stem+'.msp.gz')
        append_candidate(source, output, args.baseline/'reference_names.csv.gz')
        source = output
        allow = restore_curated_trienes(args.pre_limit_source, output) | {IDENTITY}
    report = rebuild_verified_index(source, args.baseline/'ms2'/(stem+'.sqlite'),
                                    args.destination/(stem+'.sqlite'), additional_identities=allow, revision=REVISION)
    (args.destination/(stem+'.audit.json')).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
