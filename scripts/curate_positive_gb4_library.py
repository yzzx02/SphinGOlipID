"""Append protonated Gb4 on the existing nonhydroxy d-Hex3Cer chain grid."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks, peak, rebuild_block, validate_block
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_gb4 import PositiveGb4, from_hex3cer_name
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import sha256


def gb4_block(model: PositiveGb4) -> list[str]:
    header = [
        f"Name: {model.name}", f"PrecursorMZ: {model.precursor_mz:.4f}",
        "PrecursorType: [M+H]+", "CompoundClass: Gb4", f"Formula: {model.formula}",
        f"Comment: source=theoretical Hex3Cer chain grid;MS1_name={model.species_name};polarity=+",
        "Num Peaks: 0",
    ]
    block = rebuild_block(header, [peak(f.mz, f.name, f.fragment_type) for f in model.fragments()])
    validate_block(block)
    return block


def curate(source: Path, destination: Path) -> dict:
    if source.resolve() == destination.resolve():
        raise ValueError("Write a separate Gb4 library and verify it before installation")
    destination.parent.mkdir(parents=True, exist_ok=True)
    models = {}
    with (destination.open("wb") as binary,
          gzip.GzipFile(filename="", fileobj=binary, mode="wb", compresslevel=9, mtime=0) as compressed,
          io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as stream):
        for block in iter_blocks(source):
            compound_class = header_value(block, "CompoundClass")
            if compound_class == "Gb4":
                raise ValueError("The source already contains Gb4")
            if compound_class == "Hex3Cer" and header_value(block, "PrecursorType") == "[M+H]+":
                model = from_hex3cer_name(header_value(block, "Name"))
                if model is not None:
                    if model.name in models:
                        raise ValueError(f"Duplicate Hex3Cer chain composition for {model.name}")
                    models[model.name] = model
            stream.write("\n".join(block) + "\n\n")
        if not models:
            raise ValueError("No protonated nonhydroxy d-Hex3Cer chain templates")
        for model in models.values():
            stream.write("\n".join(gb4_block(model)) + "\n\n")

    previous_records = 0
    previous_digest = hashlib.sha256()
    revised = iter_blocks(destination)
    for block in iter_blocks(source):
        if next(revised, None) != block:
            raise ValueError(f"Existing record content or order changed at {block[0]}")
        previous_digest.update(("\n".join(block) + "\n\n").encode("utf-8"))
        previous_records += 1
    for model in models.values():
        if next(revised, None) != gb4_block(model):
            raise ValueError(f"Unexpected Gb4 record for {model.name}")
    if next(revised, None) is not None:
        raise ValueError("Unexpected trailing Gb4 records")
    return dict(added_records=len(models), sum_compositions=len({m.species_name for m in models.values()}),
                preserved_records=previous_records, preserved_records_sha256=previous_digest.hexdigest(),
                source_before_sha256=sha256(source), source_sha256=sha256(destination),
                source_bytes=destination.stat().st_size)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    audit = json.dumps(curate(args.source, args.destination), indent=2) + "\n"
    args.destination.with_suffix(".curation.json").write_text(audit, encoding="utf-8")
    print(audit, end="")
