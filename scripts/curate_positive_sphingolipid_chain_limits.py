"""Apply branch positive sphingolipid chain DB limits in a separate MSP."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import (
    chain_double_bond_violations, impossible_two_chain_sum_double_bonds, is_positive_sphingolipid,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

try:
    from scripts.extend_positive_cer_sodium_library import seed_compositions
    from scripts.curate_positive_cer1p_library import sodium_composition
except ModuleNotFoundError:
    from extend_positive_cer_sodium_library import seed_compositions
    from curate_positive_cer1p_library import sodium_composition


def explicit_reasons(block):
    cls, adduct = header_value(block, 'CompoundClass'), header_value(block, 'PrecursorType')
    if not is_positive_sphingolipid(cls, adduct):
        return ()
    return chain_double_bond_violations(header_value(block, 'Name'), cls)


def supported_sodium_compositions(source):
    supported = {'Cer': set(), 'Cer1P': set()}
    for block in iter_blocks(source):
        if explicit_reasons(block):
            continue
        supported['Cer'].update(c.name for c in seed_compositions(block))
        composition = sodium_composition(block)
        if composition is not None:
            supported['Cer1P'].add(composition.name)
    return supported


def removal_reasons(block, supported):
    reasons = explicit_reasons(block)
    if reasons:
        return reasons
    cls, adduct, name = (header_value(block, key) for key in ('CompoundClass', 'PrecursorType', 'Name'))
    if adduct == '[M+Na]+' and cls in supported and '/' not in name:
        return () if name in supported[cls] else ('sodium_total_without_eligible_source',)
    if cls == 'SM' and adduct == '[M+Na]+' and impossible_two_chain_sum_double_bonds(name):
        return ('total_double_bonds_above_two_chain_maximum',)
    return ()


def curate(source, destination):
    source, destination = Path(source), Path(destination)
    originals = {Path(__file__).resolve().parents[1] / 'libraries/ms2' / name
                 for name in ('current_positive.msp.gz', 'current_negative.msp.gz')}
    if source.resolve() == destination.resolve() or destination.resolve() in {p.resolve() for p in originals}:
        raise ValueError('Write to a separate destination; never overwrite an original library')
    before_sha = sha256(source)
    supported = supported_sodium_compositions(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    removed, reasons, examples = Counter(), Counter(), []
    retained = total = 0
    digest = hashlib.sha256()
    with (destination.open('wb') as binary,
          gzip.GzipFile(filename='', fileobj=binary, mode='wb', compresslevel=9, mtime=0) as compressed,
          io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as output):
        for block in iter_blocks(source):
            total += 1
            why = removal_reasons(block, supported)
            if why:
                removed[header_value(block, 'CompoundClass') + ' ' + header_value(block, 'PrecursorType')] += 1
                reasons.update(why)
                if len(examples) < 20:
                    examples.append(dict(name=header_value(block, 'Name'), reasons=list(why)))
                continue
            text = '\n'.join(block) + '\n\n'
            output.write(text); digest.update(text.encode('utf-8')); retained += 1
    revised = iter(iter_blocks(destination)); verified = 0
    for block in iter_blocks(source):
        if removal_reasons(block, supported):
            continue
        if next(revised, None) != block:
            raise ValueError('A retained record, identity, fragment or order changed')
        verified += 1
    if next(revised, None) is not None or verified != retained or sha256(source) != before_sha:
        raise ValueError('Unexpected record addition or original source mutation')
    return dict(original_records=total, retained_records=retained, removed_records=total-retained,
        removed_records_by_class_adduct=dict(removed), removal_reason_counts=dict(reasons), removed_examples=examples,
        supported_sodium_totals={k: len(v) for k, v in supported.items()},
        retained_records_sha256=digest.hexdigest(), retained_records_content_and_order_unchanged=True,
        source_before_sha256=before_sha, source_sha256=sha256(destination), source_bytes=destination.stat().st_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    result = curate(args.source, args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
