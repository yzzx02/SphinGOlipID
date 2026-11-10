"""Streaming MSP primitives shared by repeatable library curation tools."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

PEAK_RE = re.compile(r'^(?P<mz>\S+)\s+(?P<intensity>\S+)\s+"(?P<name>[^"]*)"\s+"(?P<type>[^"]*)"$')

def header_value(lines: list[str], key: str) -> str:
    prefix = f"{key}:"
    for line in lines:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def parse_peaks(lines: list[str]) -> list[dict[str, object]]:
    peaks = []
    for line in lines:
        match = PEAK_RE.match(line)
        if match is None:
            continue
        peaks.append(
            {
                "mz": float(match.group("mz")),
                "intensity": float(match.group("intensity")),
                "name": match.group("name"),
                "type": match.group("type"),
            }
        )
    return peaks


def replace_header(block: list[str], key: str, value: str) -> list[str]:
    prefix = f"{key}:"
    result = list(block)
    for index, line in enumerate(result):
        if line.startswith(prefix):
            result[index] = f"{key}: {value}"
            return result
    return result


def format_peak(item: dict[str, object]) -> str:
    return (
        f'{float(item["mz"]):.4f} {float(item.get("intensity", 100.0)):.2f} '
        f'"{item["name"]}" "{item["type"]}"'
    )


def rebuild_block(
    block: list[str],
    peaks: list[dict[str, object]],
    replacements: dict[str, str] | None = None,
) -> list[str]:
    header = [line for line in block if PEAK_RE.match(line) is None]
    for key, value in (replacements or {}).items():
        header = replace_header(header, key, value)
    header = replace_header(header, "Num Peaks", str(len(peaks)))
    peaks = sorted(peaks, key=lambda item: float(item["mz"]))
    return header + [format_peak(item) for item in peaks]


def iter_blocks(path: Path):
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if block:
                yield block
                block = []
        if block:
            yield block


def validate_block(block: list[str]) -> None:
    for key in ("Name", "PrecursorMZ", "PrecursorType", "CompoundClass", "Num Peaks"):
        if sum(line.startswith(key + ":") for line in block) != 1:
            raise ValueError(f"Duplicate/missing {key}: {block[0]}")
    peaks = parse_peaks(block)
    if len(peaks) != int(header_value(block, "Num Peaks")) or not peaks:
        raise ValueError(f"Peak count mismatch or empty record: {block[0]}")



def peak(mz: float, name: str, fragment_type: str, intensity: float = 100.0) -> dict[str, object]:
    return {"mz": mz, "intensity": intensity, "name": name, "type": fragment_type}


def role(source: dict[str, object], fragment_type: str) -> dict[str, object]:
    return {**source, "type": fragment_type}


