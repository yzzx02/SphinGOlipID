from __future__ import annotations

import re


_AHEXCER_MSDIAL_POSITIVE_RE = re.compile(
    r"^AHexCer\s+\((?P<o_acyl>O-\d+:\d+)\)"
    r"(?P<lcb>\d+:\d+);(?P<lcb_oxygen>\d*)O/"
    r"(?P<n_acyl>\d+:\d+);(?P<n_acyl_oxygen>\d*)O$",
    flags=re.IGNORECASE,
)
_AHEXCER_MSDIAL_NEGATIVE_RE = re.compile(
    r"^AHexCer\((?P<o_acyl>\d+:\d+)/"
    r"(?P<lcb>\d+:\d+);(?P<lcb_oxygen>\d*)O/"
    r"(?P<n_acyl>\d+:\d+);(?P<n_acyl_oxygen>\d*)O\)$",
    flags=re.IGNORECASE,
)
_AHEXCER_CANONICAL_RE = re.compile(
    r"^AHexCer\s+[mdt]\d+:\d+\(O-\d+:\d+\)/"
    r"\d+:\d+(?:\((?:\d+)?OH\))?$",
    flags=re.IGNORECASE,
)
_AHEXCER_SPECIES_SOURCE_RE = re.compile(
    r"^AHexCer\((?P<o_acyl>\d+:\d+)/(?P<cer_total>\d+:\d+;\d*O)\)$",
    flags=re.IGNORECASE,
)
_ESTERIFIED_CER_CANONICAL_RE = re.compile(
    r"^Cer\s*"
    r"(?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)\((?P<o_acyl>O-\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_TOTAL_WITH_FA_RE = re.compile(
    r"^ASM\s+(?P<sm_total>\d+:\d+;\d*O)\(FA\s+(?P<o_acyl>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_HYDROXY_FA_SOURCE_RE = re.compile(
    r"^(?P<class>Cer|HexCer)\((?P<lcb>[mdt]\d+:\d+)/"
    r"(?P<fa>\d+:\d+)\)\(OH\)$",
    flags=re.IGNORECASE,
)
_HYDROXY_FA_CANONICAL_RE = re.compile(
    r"^(?P<class>Cer|HexCer)\((?P<lcb>[mdt]\d+:\d+)/"
    r"h(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_OXYGENATED_CONJUGATED_CER_SOURCE_RE = re.compile(
    r"^(?P<class>PE-Cer\+O|PI-Cer\+O)\("
    r"(?P<lcb>\d+:\d+);2O/(?P<fa>\d+:\d+);O\)$",
    flags=re.IGNORECASE,
)
_OXYGENATED_CONJUGATED_CER_CANONICAL_RE = re.compile(
    r"^(?P<class>PE-Cer\+O|PI-Cer\+O)\("
    r"m(?P<lcb>\d+:\d+)/(?P<fa>\d+:\d+);2O\)$",
    flags=re.IGNORECASE,
)
_SHEXCER_CANONICAL_RE = re.compile(
    r"^SHexCer\((?P<lcb>[mdt]\d+:\d+)/(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SHEXCER_REVERSED_RE = re.compile(
    r"^SHexCer\((?P<fa>\d+:\d+)/(?P<lcb>[mdt]\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SHEXCER_O_SOURCE_RE = re.compile(
    r"^SHexCer\+O\((?P<lcb>[mdt]\d+:\d+)/(?P<fa>\d+:\d+)\)\(OH\)$",
    flags=re.IGNORECASE,
)
_SHEXCER_O_MALFORMED_RE = re.compile(
    r"^SHexCer\+O\((?P<fa>\d+:\d+)\)\(OH/(?P<lcb>[mdt]\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SHEXCER_O_CANONICAL_RE = re.compile(
    r"^SHexCer\+O\((?P<lcb>[mdt]\d+:\d+)/h(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SL_CANONICAL_RE = re.compile(
    r"^SL\(m(?P<lcb>\d+:\d+)/(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SL_SOURCE_RE = re.compile(
    r"^SL\((?P<lcb>\d+:\d+);O/(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SL_REVERSED_M_RE = re.compile(
    r"^SL\((?P<fa>\d+:\d+)/m(?P<lcb>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_SL_REVERSED_O_RE = re.compile(
    r"^SL\((?P<fa>\d+:\d+)/(?P<lcb>\d+:\d+);O\)$",
    flags=re.IGNORECASE,
)
_SL_O_SOURCE_RE = re.compile(
    r"^SL\+O\((?P<lcb>\d+:\d+);O/(?P<fa>\d+:\d+);O\)$",
    flags=re.IGNORECASE,
)
_SL_O_CANONICAL_RE = re.compile(
    r"^SL\+O\(m(?P<lcb>\d+:\d+)/h(?P<fa>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_COLLAPSED_SOURCE_RE = re.compile(
    r"^ASM\s+d(?P<core_c>\d+):(?P<core_db>\d+)/"
    r"(?P<o_acyl>\d+:\d+)$",
    flags=re.IGNORECASE,
)
_ASM_COLLAPSED_CANONICAL_RE = re.compile(
    r"^ASM\s+d(?P<core_c>\d+):(?P<core_db>\d+)"
    r"\(O-(?P<o_acyl>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_COMPLETE_SOURCE_RE = re.compile(
    r"^ASM(?:\s+|\()(?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)"
    r"(?:-(?P<marked_o_acyl>O-\d+:\d+)|\(FA\s+(?P<fa>\d+:\d+)\))\)?$",
    flags=re.IGNORECASE,
)
_ASM_TWO_CHAIN_PAREN_RE = re.compile(
    r"^ASM\((?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
def _lcb_prefix(oxygen_count: str) -> str:
    return {1: "m", 2: "d", 3: "t"}.get(int(oxygen_count or "1"), "d")


def _hydroxy_suffix(oxygen_count: str) -> str:
    count = int(oxygen_count or "1")
    return "(OH)" if count == 1 else f"({count}OH)"


def canonicalize_multichain_sphingolipid_name(
    value: object,
    compound_class: object = "",
) -> str:
    """Normalize sphingolipids that carry an additional esterified fatty acid.

    The extra fatty acid is written explicitly at its attachment site. ASM is
    kept at the summed sphingomyelin-core level, ``ASM dX:Y(O-C:DB)``, rather
    than inventing a particular LCB/N-acyl pair.
    """

    text = str(value or "").strip()
    cls = str(compound_class or "").strip().upper()
    if not cls:
        cls = text.split("(", 1)[0].split(" ", 1)[0].strip().upper()

    esterified = _ESTERIFIED_CER_CANONICAL_RE.fullmatch(text)
    if esterified is not None:
        return f"Cer{esterified.group('lcb')}/{esterified.group('n_acyl')}({esterified.group('o_acyl')})"

    if cls == "AHEXCER":
        if _AHEXCER_CANONICAL_RE.fullmatch(text):
            return text
        matched = _AHEXCER_MSDIAL_POSITIVE_RE.fullmatch(text)
        if matched is None:
            matched = _AHEXCER_MSDIAL_NEGATIVE_RE.fullmatch(text)
        if matched is not None:
            o_acyl = matched.group("o_acyl")
            if not o_acyl.upper().startswith("O-"):
                o_acyl = f"O-{o_acyl}"
            return (
                f"AHexCer {_lcb_prefix(matched.group('lcb_oxygen'))}{matched.group('lcb')}"
                f"({o_acyl})/{matched.group('n_acyl')}"
                f"{_hydroxy_suffix(matched.group('n_acyl_oxygen'))}"
            )
        matched_species = _AHEXCER_SPECIES_SOURCE_RE.fullmatch(text)
        if matched_species is not None:
            return (
                f"AHexCer {matched_species.group('cer_total')}"
                f"(O-{matched_species.group('o_acyl')})"
            )
        return text

    if cls in {"CER", "HEXCER"}:
        matched = _HYDROXY_FA_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _HYDROXY_FA_SOURCE_RE.fullmatch(text)
        if matched is not None:
            class_name = "HexCer" if matched.group("class").upper() == "HEXCER" else "Cer"
            return f"{class_name}({matched.group('lcb')}/h{matched.group('fa')})"
        return text

    if cls in {"PE-CER+O", "PI-CER+O"}:
        matched = _OXYGENATED_CONJUGATED_CER_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _OXYGENATED_CONJUGATED_CER_SOURCE_RE.fullmatch(text)
        if matched is not None:
            class_name = "PI-Cer+O" if matched.group("class").upper() == "PI-CER+O" else "PE-Cer+O"
            return f"{class_name}(m{matched.group('lcb')}/{matched.group('fa')};2O)"
        return text

    if cls == "SHEXCER":
        matched = _SHEXCER_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _SHEXCER_REVERSED_RE.fullmatch(text)
        if matched is not None:
            return f"SHexCer({matched.group('lcb')}/{matched.group('fa')})"
        return text

    if cls == "SHEXCER+O":
        matched = _SHEXCER_O_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _SHEXCER_O_SOURCE_RE.fullmatch(text)
        if matched is None:
            matched = _SHEXCER_O_MALFORMED_RE.fullmatch(text)
        if matched is not None:
            return f"SHexCer+O({matched.group('lcb')}/h{matched.group('fa')})"
        return text

    if cls == "SL":
        matched = _SL_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _SL_SOURCE_RE.fullmatch(text)
        if matched is None:
            matched = _SL_REVERSED_M_RE.fullmatch(text)
        if matched is None:
            matched = _SL_REVERSED_O_RE.fullmatch(text)
        if matched is not None:
            return f"SL(m{matched.group('lcb')}/{matched.group('fa')})"
        return text

    if cls == "SL+O":
        matched = _SL_O_CANONICAL_RE.fullmatch(text)
        if matched is None:
            matched = _SL_O_SOURCE_RE.fullmatch(text)
        if matched is not None:
            return f"SL+O(m{matched.group('lcb')}/h{matched.group('fa')})"
        return text

    if cls == "ASM":
        matched = _ASM_COLLAPSED_CANONICAL_RE.fullmatch(text)
        if matched is not None:
            return (
                f"ASM d{matched.group('core_c')}:{matched.group('core_db')}"
                f"(O-{matched.group('o_acyl')})"
            )
        matched = _ASM_COLLAPSED_SOURCE_RE.fullmatch(text)
        if matched is not None:
            return (
                f"ASM d{matched.group('core_c')}:{matched.group('core_db')}"
                f"(O-{matched.group('o_acyl')})"
            )
        matched = _ASM_TOTAL_WITH_FA_RE.fullmatch(text)
        if matched is not None:
            core_total = matched.group("sm_total").split(";", 1)[0]
            return f"ASM d{core_total}(O-{matched.group('o_acyl')})"
        matched = _ASM_COMPLETE_SOURCE_RE.fullmatch(text)
        if matched is not None:
            lcb_match = re.fullmatch(r"[mdt](\d+):(\d+)", matched.group("lcb"), flags=re.IGNORECASE)
            n_acyl_match = re.fullmatch(r"(\d+):(\d+)", matched.group("n_acyl"))
            if lcb_match is not None and n_acyl_match is not None:
                core_c = int(lcb_match.group(1)) + int(n_acyl_match.group(1))
                core_db = int(lcb_match.group(2)) + int(n_acyl_match.group(2))
                o_acyl = matched.group("marked_o_acyl") or matched.group("fa")
                o_acyl = str(o_acyl).removeprefix("O-")
                return f"ASM d{core_c}:{core_db}(O-{o_acyl})"
        matched = _ASM_TWO_CHAIN_PAREN_RE.fullmatch(text)
        if matched is not None:
            return f"ASM {matched.group('lcb')}/{matched.group('n_acyl')}"
        return text

    return text


def has_complete_multichain_sphingolipid_identity(
    value: object,
    compound_class: object = "",
) -> bool:
    """Return true only when all three structural chains are explicit."""

    text = canonicalize_multichain_sphingolipid_name(value, compound_class)
    cls = str(compound_class or "").strip().upper()
    if not cls:
        cls = text.split("(", 1)[0].split(" ", 1)[0].strip().upper()
    if cls == "AHEXCER":
        return _AHEXCER_CANONICAL_RE.fullmatch(text) is not None
    if cls == "CER" or _ESTERIFIED_CER_CANONICAL_RE.fullmatch(text):
        return _ESTERIFIED_CER_CANONICAL_RE.fullmatch(text) is not None
    if cls == "ASM":
        return False
    if cls == "SHEXCER":
        return _SHEXCER_CANONICAL_RE.fullmatch(text) is not None
    if cls == "SHEXCER+O":
        return _SHEXCER_O_CANONICAL_RE.fullmatch(text) is not None
    if cls == "SL":
        return _SL_CANONICAL_RE.fullmatch(text) is not None
    if cls == "SL+O":
        return _SL_O_CANONICAL_RE.fullmatch(text) is not None
    return False
