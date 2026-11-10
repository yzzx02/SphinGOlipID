from __future__ import annotations

import re


def canonicalize_single_chain_name(value: object, compound_class: object = "") -> str:
    """Remove the empty glycerol position from LCE-PE display identities."""
    text = str(value or "").strip()
    if str(compound_class or "").strip().upper() != "LCE-PE":
        return text
    match = re.fullmatch(r"LCE-PE\(([^()/]+)[/_]([^()/]+)\)", text, re.IGNORECASE)
    if match and (match[1] == "0:0") != (match[2] == "0:0"):
        return f"LCE-PE({match[2] if match[1] == '0:0' else match[1]})"
    return text


def canonicalize_n_acyl_glycerophospholipid_name(
    value: object,
    compound_class: object = "",
) -> str:
    """Keep glycerol chains before the N-acyl chain in NAPE/NAPS names."""

    text = canonicalize_single_chain_name(value, compound_class)
    class_name = str(compound_class or "").strip().upper()
    if class_name not in {"NAPE", "NAPS"}:
        return text
    matched = re.fullmatch(
        r"(?P<prefix>NAPE|NAPS)\((?P<inner>[^()]*)\)",
        text,
        flags=re.IGNORECASE,
    )
    if matched is None:
        return text
    chain = r"(?:O-|P-)?\d+:\d+"
    inner = matched.group("inner")
    glycerol_first = re.fullmatch(
        rf"(?P<g1>{chain})_(?P<g2>{chain})-N-(?P<n>{chain})",
        inner,
        flags=re.IGNORECASE,
    )
    n_acyl_first = re.fullmatch(
        rf"(?P<n>{chain})-N-(?P<g1>{chain})_(?P<g2>{chain})",
        inner,
        flags=re.IGNORECASE,
    )
    parsed = glycerol_first or n_acyl_first
    if parsed is None:
        # Total-DAG notation such as NAPS(36:1-N-16:0) is already canonical.
        return text
    glycerol_chains = sorted((parsed.group("g1"), parsed.group("g2")))
    return (
        f"{matched.group('prefix').upper()}("
        f"{'_'.join(glycerol_chains)}-N-{parsed.group('n')})"
    )
