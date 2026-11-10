"""Append total-composition Cer [M+Na]+ spectra to a separate positive library."""
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks, validate_block
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_cer_sodium import CerSodiumComposition
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import explicit_acyl_carbons, chain_double_bond_violations
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256

REVISION = "positive_cer_sodium_sum_20261108"
CHAIN = re.compile(r"Cer\(([mdt])(\d+):(\d+)/(h?)(\d+):(\d+)\)")


def seed_compositions(block):
    if (header_value(block, "CompoundClass"), header_value(block, "PrecursorType")) != ("Cer", "[M+H]+"):
        return ()
    name = header_value(block, "Name")
    match = CHAIN.fullmatch(name)
    if match is None or any(c < 14 for c in explicit_acyl_carbons(name)) or chain_double_bond_violations(name, 'Cer'):
        return ()
    series, lc, ld, hydroxy, fc, fd = match.groups()
    carbon, db = int(lc)+int(fc), int(ld)+int(fd)
    oh = {"m": 1, "d": 2, "t": 3}[series] + bool(hydroxy)
    # Cover both ordinary and mono-hydroxy N-acyl totals on the existing chain
    # grid. Merge by total composition; never assign the OH site in a Na result.
    hydroxyls = (oh,) if hydroxy else (oh, oh+1)
    return tuple(CerSodiumComposition(carbon, db, h) for h in hydroxyls)


def sodium_block(composition):
    fragments = composition.fragments()
    block = [f"Name: {composition.name}", f"PrecursorMZ: {composition.precursor_mz:.4f}",
             "PrecursorType: [M+Na]+", "CompoundClass: Cer", f"Formula: {composition.formula}",
             f"Comment: MS1_name={composition.name};polarity=+", "Source: theoretical",
             f"LibraryRevision: {REVISION}", "AnnotationScope: sum_composition_only",
             f"Num Peaks: {len(fragments)}"]
    block += [f'{f.mz:.4f} {f.intensity:.2f} "{f.name}" "{f.fragment_type}"' for f in fragments]
    validate_block(block)
    return block


def extend(source: Path, destination: Path):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Write to a separate destination before installing the library")
    destination.parent.mkdir(parents=True, exist_ok=True)
    compositions, existing = set(), set()
    preserved = hashlib.sha256()
    count = 0
    with (destination.open("wb") as binary,
          gzip.GzipFile(filename="", fileobj=binary, mode="wb", compresslevel=9, mtime=0) as compressed,
          io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as output):
        for block in iter_blocks(source):
            compositions.update(seed_compositions(block))
            if (header_value(block, "CompoundClass"), header_value(block, "PrecursorType")) == ("Cer", "[M+Na]+"):
                existing.add(header_value(block, "Name"))
            text = "\n".join(block) + "\n\n"
            preserved.update(text.encode("utf-8"))
            output.write(text)
            count += 1
        additions = [sodium_block(c) for c in sorted(compositions) if c.name not in existing]
        for block in additions:
            output.write("\n".join(block) + "\n\n")
    revised = iter(iter_blocks(destination))
    for before in iter_blocks(source):
        if next(revised, None) != before:
            raise ValueError("Original library content or record order changed")
    for block in additions:
        if next(revised, None) != block:
            raise ValueError("Appended sodium spectrum differs from expected content")
    if next(revised, None) is not None:
        raise ValueError("Unexpected additional library record")
    return dict(revision=REVISION, added_records=len(additions),
                added_by_hydroxyls=dict(Counter(c.hydroxyls for c in compositions if c.name not in existing)),
                added_names=[header_value(b, "Name") for b in additions],
                preserved_records=count, preserved_records_sha256=preserved.hexdigest(),
                source_before_sha256=sha256(source), source_sha256=sha256(destination),
                source_bytes=destination.stat().st_size)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    result = extend(args.source, args.destination)
    args.destination.with_suffix(".extension.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "added_names"}, indent=2))
