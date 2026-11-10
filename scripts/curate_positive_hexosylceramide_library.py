"""Revise protonated hexosylceramide fragments in a separate, verified MSP."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
from itertools import zip_longest
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks, parse_peaks, rebuild_block, validate_block
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import CLASSES, is_positive_hexosylceramide, normalize_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256


def curate_block(block, classes=CLASSES):
    compound_class = header_value(block, "CompoundClass")
    if compound_class not in classes or not is_positive_hexosylceramide(
        compound_class, header_value(block, "PrecursorType"),
    ):
        return block
    validate_block(block)
    source = [
        FragmentRecord(p["mz"], p["name"], p["type"], intensity=p["intensity"])
        for p in parse_peaks(block)
    ]
    fragments = normalize_fragments(
        compound_class, header_value(block, "Name"), float(header_value(block, "PrecursorMZ")), source,
    )
    revised = rebuild_block(block, [
        dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity) for f in fragments
    ])
    validate_block(revised)
    return revised


def curate(source: Path, destination: Path, classes=CLASSES) -> dict:
    if source.resolve() == destination.resolve():
        raise ValueError("Write to a separate destination before installing the library")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with (destination.open("wb") as binary,
          gzip.GzipFile(filename="", fileobj=binary, mode="wb", compresslevel=9, mtime=0) as compressed,
          io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as stream):
        for block in iter_blocks(source):
            stream.write("\n".join(curate_block(block, classes)) + "\n\n")

    changed = Counter()
    unchanged = 0
    untouched = hashlib.sha256()
    for before, after in zip_longest(iter_blocks(source), iter_blocks(destination)):
        if before is None or after is None or after != curate_block(before, classes):
            raise ValueError("Library record content or order differs from the expected revision")
        if before != after:
            changed[header_value(before, "CompoundClass")] += 1
        else:
            unchanged += 1
            untouched.update(("\n".join(before) + "\n\n").encode("utf-8"))
    return dict(changed_records=dict(changed), preserved_records=unchanged,
                preserved_records_sha256=untouched.hexdigest(), source_before_sha256=sha256(source),
                source_sha256=sha256(destination), source_bytes=destination.stat().st_size)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--class", dest="classes", action="append", choices=sorted(CLASSES))
    args = parser.parse_args()
    audit = json.dumps(curate(args.source, args.destination, set(args.classes or CLASSES)), indent=2) + "\n"
    args.destination.with_suffix(".curation.json").write_text(audit, encoding="utf-8")
    print(audit, end="")
