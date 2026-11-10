from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
import sys
import pickle
import logging
import shutil
import tempfile
import re
from dataclasses import replace
from pathlib import Path
from typing import Dict, Iterable, Iterator, List

import pandas as pd

from . import safe_pickle
from .esterified_ceramide import is_esterified_ceramide
from .negative_hexcer import is_negative_hexcer, normalize_negative_hexcer_fragments
from .negative_cer_hsu2016 import normalize_negative_cer_hsu2016_fragments

from .import_msdial_sphingo_positive import annotate_positive_sl_fragments
from .lipid_naming import canonicalize_n_acyl_glycerophospholipid_name
from .negative_ptdglc import canonicalize_glc_gp_identity
from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS,
    HYDROGEN_MONOISOTOPIC_MASS,
    OXYGEN_MONOISOTOPIC_MASS,
    PROTON_MONOISOTOPIC_MASS,
    WATER_MONOISOTOPIC_MASS,
    fragment_with_role as _fragment_with_role,
    normalize_special_aminophospholipid_fragments as _normalize_special_aminophospholipid_fragments,
    normalize_tg_est_fragments as _normalize_tg_est_fragments,
    normalize_tg_o_fragments as _normalize_tg_o_fragments,
    curate_positive_lyso_pe_fragments as _curate_positive_lyso_pe_fragments,
    curate_positive_pe_p_fragments as _curate_positive_pe_p_fragments,
    without_pi_ammonium_acyl_losses as _without_pi_ammonium_acyl_losses,
    synthetic_fragment as _synthetic_fragment,
)
from .models import FragmentRecord, LibraryRecord
from .nae_precursor import normalize_nae_record
from .negative_lpe import normalize_negative_lpe_record
from .negative_lpc import normalize_negative_lpc_record
from .positive_lpc_ether import normalize_positive_ether_lpc
from .positive_single_chain_h import normalize_protonated_single_chain_record
from .positive_lpc_sodium import normalize_lpc_sodium
from .positive_spb import normalize_fragments as _normalize_positive_spb_fragments
from .positive_complex_ceramide import normalize_fragments as _normalize_complex_ceramide_fragments
from .positive_cer1p import normalize_h_fragments as _normalize_cer1p_h_fragments
from .positive_ganglioside import (
    is_positive_ganglioside as is_positive_extended_ganglioside,
    normalize_fragments as _normalize_extended_ganglioside_fragments,
)
from .positive_hexosylceramide import (
    is_positive_hexosylceramide,
    normalize_fragments as _normalize_positive_hexosyl_fragments,
)
from .chain_labels import annotate_record_chain_labels
from .provenance import sha256, code_fingerprint
from .sphingolipid_naming import (
    canonicalize_multichain_sphingolipid_name,
    has_complete_multichain_sphingolipid_identity,
)


LIBRARY_CACHE_VERSION = 43

TRIMETHYLAMINE_MONOISOTOPIC_MASS = 59.07349929
PHOSPHORIC_ACID_MONOISOTOPIC_MASS = 97.97689557
HYDROXY_FA_SPHINGOLIPID_CLASSES = {"Cer", "HexCer", "LacCer", "Hex2Cer"}
LCB_FRAGMENT_TYPE = "LCB碎片"
D18_1_T18_0_MARKER_MZ = 300.2897

POSITIVE_CHOLINE_COMMON_FRAGMENTS = (
    (86.0964, "[C5H12N]+"),
    (104.1070, "[C5H14NO]+"),
    (124.9998, "[C2H6O4P]+"),
)
POSITIVE_CHOLINE_CLASS_KEYS = {
    "PC",
    "PCO",
    "PCP",
    "LPC",
    "LPCO",
    "SM",
}
FREE_SPHINGOID_BASE_CLASSES = {"SPB"}

def _canonicalize_ahexcer_name(value: object) -> str | None:
    text = str(value or "").strip()
    canonical = canonicalize_multichain_sphingolipid_name(text, "AHexCer")
    if canonical == text and not has_complete_multichain_sphingolipid_identity(
        canonical,
        "AHexCer",
    ):
        return None
    return canonical


def _canonicalize_sphingoid_base_identity(
    compound_class: object,
    lipid_name: object,
    lipid_chain_name: object,
) -> tuple[str, str, str]:
    cls = str(compound_class or "").strip()
    name = canonicalize_multichain_sphingolipid_name(lipid_name, cls)
    chain_name = canonicalize_multichain_sphingolipid_name(lipid_chain_name, cls)
    if is_esterified_ceramide(cls, chain_name):
        return "Cer", chain_name, chain_name
    if cls == "AHexCer":
        canonical_ahexcer = _canonicalize_ahexcer_name(chain_name)
        if canonical_ahexcer is not None:
            return cls, canonical_ahexcer, canonical_ahexcer
    if cls not in FREE_SPHINGOID_BASE_CLASSES:
        return cls, name, chain_name

    def canonicalize_name(value: str) -> str:
        if "(" in value:
            return f"SPB({value.split('(', 1)[1]}"
        return value

    return "SPB", canonicalize_name(name), canonicalize_name(chain_name)


def _add_formula_oxygen(formula: str) -> str:
    text = str(formula or "").strip()
    if not text:
        return text
    match = re.search(r"O(\d*)", text)
    if match is None:
        return f"{text}O"
    current_count = int(match.group(1) or "1")
    return f"{text[:match.start()]}O{current_count + 1}{text[match.end():]}"


def _without_d18_1_300_lcb_marker(record: LibraryRecord) -> LibraryRecord:
    """Keep m/z 300.2897 as a t18:0 marker instead of a d18:1 LCB ion."""

    if re.search(r"\(d18:1/(?:h)?\d+:\d+\)", record.lipid_chain_name, flags=re.IGNORECASE) is None:
        return record
    fragments = [
        fragment
        for fragment in record.fragments
        if not (
            (fragment.fragment_type == LCB_FRAGMENT_TYPE or (
                (is_positive_hexosylceramide(record.compound_class, record.adduct)
                 or (record.compound_class == "SM" and record.adduct == "[M+H]+"))
                and fragment.name.startswith("LCB")
            ))
            and abs(float(fragment.mz) - D18_1_T18_0_MARKER_MZ) <= 0.02
        )
    ]
    if len(fragments) == len(record.fragments):
        return record
    return replace(record, fragments=fragments)


def _fatty_acid_anion_mz(carbons: int, double_bonds: int) -> float:
    neutral_hydrogens = 2 * int(carbons) - 2 * int(double_bonds)
    return (
        int(carbons) * CARBON_MONOISOTOPIC_MASS
        + neutral_hydrogens * HYDROGEN_MONOISOTOPIC_MASS
        + 2 * OXYGEN_MONOISOTOPIC_MASS
        - PROTON_MONOISOTOPIC_MASS
    )


def _normalize_targeted_positive_sphingolipid_fragments(
    compound_class: str,
    lipid_chain_name: str,
    precursor_mz: float,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> List[FragmentRecord]:
    """Apply curated evidence pools for the targeted sphingolipid classes."""

    result = list(fragments)
    adduct_text = str(adduct or "").strip()
    cls = str(compound_class or "").strip()
    chain_name = str(lipid_chain_name or "").strip()
    if is_esterified_ceramide(cls, chain_name):
        return result
    if is_positive_extended_ganglioside(cls,adduct_text):
        return _normalize_extended_ganglioside_fragments(cls,chain_name,adduct_text,result)
    if is_negative_hexcer(cls, adduct_text):
        return normalize_negative_hexcer_fragments(chain_name, precursor_mz, adduct_text, result)
    corrected_negative_cer = normalize_negative_cer_hsu2016_fragments(
        chain_name, precursor_mz, adduct_text, result,
    )
    if corrected_negative_cer != result:
        return corrected_negative_cer
    by_name = {str(fragment.name).strip(): fragment for fragment in result}

    if adduct_text in {"[M+CH3COO]-", "[M+HCOO]-"} and cls == "SM":
        fa_match = re.search(r"/(?P<carbons>\d+):(?P<double_bonds>\d+)\)", chain_name)
        phosphate = by_name.get("PO3-") or _synthetic_fragment(78.9591, "PO3-", "Common")
        phosphocholine = by_name.get("[C4H11NO4P]-") or _synthetic_fragment(
            168.0431,
            "[C4H11NO4P]-",
            "Diagnostic_HG",
        )
        methyl_loss = by_name.get("M-CH3")
        precursor = by_name.get(adduct_text) or _synthetic_fragment(
            float(precursor_mz),
            adduct_text,
            "Precursor Ion",
        )
        curated = [
            _fragment_with_role(phosphate, "Common"),
            _fragment_with_role(phosphocholine, "Diagnostic_HG"),
            _fragment_with_role(precursor, "Precursor Ion"),
        ]
        if methyl_loss is not None:
            curated.append(_fragment_with_role(methyl_loss, "Diagnostic_HG"))
        fa_loss = next(
            (
                fragment
                for fragment in result
                if fragment.fragment_type == "Diagnostic_FA_Loss"
                or str(fragment.name or "").startswith("M-CH3-(R=O)(")
            ),
            None,
        )
        if fa_loss is not None:
            curated.append(_fragment_with_role(fa_loss, "Diagnostic_FA_Loss"))
        if fa_match is not None:
            carbons = int(fa_match.group("carbons"))
            double_bonds = int(fa_match.group("double_bonds"))
            fa_token = f"{carbons}:{double_bonds}"
            rcoo_name = f"[RCOO]-({fa_token})"
            rcoo = by_name.get(rcoo_name) or _synthetic_fragment(
                _fatty_acid_anion_mz(carbons, double_bonds),
                rcoo_name,
                "Diagnostic_FA",
            )
            curated.append(_fragment_with_role(rcoo, "Diagnostic_FA"))
        return curated

    if adduct_text == "[M-H]-" and cls == "Cer1P":
        phosphate = by_name.get("PO3-") or _synthetic_fragment(
            78.9591,
            "PO3-",
            "Diagnostic_HG",
        )
        hydrogen_phosphate = by_name.get("H2PO4-") or _synthetic_fragment(
            96.9696,
            "H2PO4-",
            "Diagnostic_HG",
        )
        precursor = by_name.get("[M-H]-") or _synthetic_fragment(
            float(precursor_mz),
            "[M-H]-",
            "Precursor Ion",
        )
        water_loss = by_name.get("M-H-H2O") or _synthetic_fragment(
            float(precursor_mz) - WATER_MONOISOTOPIC_MASS,
            "M-H-H2O",
            "Common",
        )
        ketene_loss = next(
            (
                fragment
                for fragment in result
                if fragment.fragment_type == "Diagnostic_FA_Loss"
                and (
                    str(fragment.name).startswith("NL_Ketene(")
                    or "M-H-(R=O)" in str(fragment.name)
                )
            ),
            None,
        )
        curated = [
            _fragment_with_role(phosphate, "Diagnostic_HG"),
            _fragment_with_role(hydrogen_phosphate, "Diagnostic_HG"),
            _fragment_with_role(precursor, "Precursor Ion"),
            _fragment_with_role(water_loss, "Common"),
        ]
        if ketene_loss is not None:
            fa_match = re.search(r"/(?:n|h)?(?P<fa>\d+:\d+)\)", chain_name, flags=re.IGNORECASE)
            ketene_name = (
                f"NL_Ketene(n{fa_match.group('fa')})"
                if fa_match is not None
                else str(ketene_loss.name)
            )
            curated.append(
                _fragment_with_role(replace(ketene_loss, name=ketene_name), "Diagnostic_FA_Loss")
            )
            ketene_water_loss = next(
                (
                    fragment
                    for fragment in result
                    if fragment.fragment_type == "Diagnostic_FA_Loss"
                    and (
                        str(fragment.name).startswith("NL_Ketene-H2O(")
                        or "M-H-(ROOH)" in str(fragment.name)
                    )
                ),
                None,
            ) or replace(
                ketene_loss,
                mz=float(ketene_loss.mz) - WATER_MONOISOTOPIC_MASS,
            )
            ketene_water_name = (
                f"NL_Ketene-H2O(n{fa_match.group('fa')})"
                if fa_match is not None
                else str(ketene_water_loss.name)
            )
            curated.append(
                _fragment_with_role(
                    replace(ketene_water_loss, name=ketene_water_name),
                    "Diagnostic_FA_Loss",
                )
            )
        return curated

    if adduct_text != "[M+H]+":
        return result

    if is_positive_hexosylceramide(cls, adduct_text) or cls in {"SM", "Cer"}:
        return _normalize_positive_hexosyl_fragments(cls, chain_name, precursor_mz, result)

    if (cls == "AHexCer" and len(result) == 10
            and all(str(fragment.name).startswith('m/z ') or fragment.name in {'[M+H]+','M+H-H2O'} for fragment in result)):
        ordered = sorted(result, key=lambda fragment: float(fragment.mz))
        o_acyl_match = re.search(r"\((O-\d+:\d+)\)", chain_name, flags=re.IGNORECASE)
        if o_acyl_match is not None:
            o_acyl = o_acyl_match.group(1)
            names_and_types = (
                ("LCB-C2H5N", LCB_FRAGMENT_TYPE),
                ("LCB-CH4O2", LCB_FRAGMENT_TYPE),
                ("LCB-2H2O", LCB_FRAGMENT_TYPE),
                ("LCB-H2O", LCB_FRAGMENT_TYPE),
                (f"{o_acyl}-Hex+", "Diagnostic_HG"),
                (f"M+H-Acyl({o_acyl})-C6H10O5-2H2O", "Diagnostic_HG"),
                (f"M+H-Acyl({o_acyl})-C6H10O5-H2O", "Diagnostic_HG"),
                (f"M+H-Acyl({o_acyl})-C6H10O5", "Diagnostic_HG"),
                ("M+H-H2O", "Common"),
                ("[M+H]+", "Common"),
            )
            return [
                _fragment_with_role(replace(fragment, name=name), fragment_type)
                for fragment, (name, fragment_type) in zip(ordered, names_and_types)
            ]

    if cls == "SPB":
        return _normalize_positive_spb_fragments(result)

    if cls == "LSM":
        hg = by_name.get("[C5H15NO4P]+") or _synthetic_fragment(
            184.0733,
            "[C5H15NO4P]+",
            "Diagnostic_HG",
        )
        water_loss = by_name.get("M+H-H2O") or _synthetic_fragment(
            float(precursor_mz) - WATER_MONOISOTOPIC_MASS,
            "M+H-H2O",
            "Common",
        )
        trimethylamine_loss = by_name.get("M+H-trimethylamine(-59)") or _synthetic_fragment(
            float(precursor_mz) - TRIMETHYLAMINE_MONOISOTOPIC_MASS,
            "M+H-trimethylamine(-59)",
            "Common",
        )
        lcb_h2o = by_name.get("LCB-H2O")
        lcb_2h2o = by_name.get("LCB-2H2O")
        curated = [
            _fragment_with_role(hg, "Diagnostic_HG"),
            _fragment_with_role(water_loss, "Common"),
            _fragment_with_role(trimethylamine_loss, "Common"),
        ]
        if lcb_h2o is not None:
            curated.append(_fragment_with_role(lcb_h2o, LCB_FRAGMENT_TYPE))
        if lcb_2h2o is not None:
            curated.append(_fragment_with_role(lcb_2h2o, LCB_FRAGMENT_TYPE))
        return curated

    if cls == "Cer1P":
        return _normalize_cer1p_h_fragments(chain_name, float(precursor_mz), result)

    return result


def _hydroxy_fa_fragment(fragment: FragmentRecord, compound_class: str = "") -> FragmentRecord:
    # Adding the FA hydroxyl shifts precursor-derived ions and the acyl-derived
    # Ceramide-U ion by one oxygen. LCB-only ions retain their original masses.
    shifts_with_fa_oxygen = fragment.fragment_type in {
        "Precursor Ion",
        "C类碎片",
        "Diagnostic_HG",
    } or fragment.name == "Ceramide fragment U" or (
        (compound_class == "Cer" and fragment.name.startswith("M+H-"))
        or
        (compound_class == "HexCer" and fragment.name in {"M+H-H2O", "M+H-C6H10O5"})
        or (compound_class in {"Hex2Cer", "LacCer"} and fragment.name.startswith("M+H-"))
    )
    if not shifts_with_fa_oxygen:
        return fragment
    return replace(fragment, mz=float(fragment.mz) + OXYGEN_MONOISOTOPIC_MASS)


def _expand_hydroxy_fa_sphingolipids(records: Iterable[LibraryRecord]) -> List[LibraryRecord]:
    """Add positive-mode d/h molecular species and normalize three-OH totals to t."""

    normalized = [_without_d18_1_300_lcb_marker(record) for record in records]
    existing_keys = {
        (
            record.compound_class,
            record.lipid_chain_name,
            round(float(record.precursor_mz), 4),
            record.adduct,
        )
        for record in normalized
    }
    next_id = max((record.record_id for record in normalized), default=-1) + 1
    generated: List[LibraryRecord] = []
    chain_pattern = re.compile(
        r"^(?P<prefix>[^()]+)\(d(?P<base_c>\d+):(?P<base_db>\d+)/"
        r"(?P<fa_c>\d+):(?P<fa_db>\d+)\)$",
        flags=re.IGNORECASE,
    )

    for source in normalized:
        if source.compound_class not in HYDROXY_FA_SPHINGOLIPID_CLASSES:
            continue
        if source.adduct != "[M+H]+":
            continue
        match = chain_pattern.fullmatch(source.lipid_chain_name)
        if match is None:
            continue

        base_c = int(match.group("base_c"))
        base_db = int(match.group("base_db"))
        fa_c = int(match.group("fa_c"))
        fa_db = int(match.group("fa_db"))
        chain_name = (
            f"{match.group('prefix')}(d{base_c}:{base_db}/h{fa_c}:{fa_db})"
        )
        precursor_mz = float(source.precursor_mz) + OXYGEN_MONOISOTOPIC_MASS
        key = (
            source.compound_class,
            chain_name,
            round(precursor_mz, 4),
            source.adduct,
        )
        if key in existing_keys:
            continue

        fragments = [_hydroxy_fa_fragment(fragment, source.compound_class) for fragment in source.fragments]
        derived = LibraryRecord(
            record_id=next_id,
            compound_class=source.compound_class,
            lipid_name=(
                f"{source.compound_class}(t{base_c + fa_c}:{base_db + fa_db})"
            ),
            lipid_chain_name=chain_name,
            precursor_mz=precursor_mz,
            adduct=source.adduct,
            formula=_add_formula_oxygen(source.formula),
            polarity=source.polarity,
            fragments=fragments,
            metadata={
                **source.metadata,
                "generated_hydroxy_fa_isomer": "true",
                "source_chain_name": source.lipid_chain_name,
            },
        )
        derived = _without_d18_1_300_lcb_marker(derived)
        generated.append(derived)
        existing_keys.add(key)
        next_id += 1

    return normalized + generated


def _library_cache_metadata(path: Path) -> Dict[str, object]:
    stat = path.stat()
    return {
        "version": LIBRARY_CACHE_VERSION,
        "path": _library_cache_identity(path),
        "size": stat.st_size,
        "mtime_ns": 0 if _is_bundled_library(path) else stat.st_mtime_ns,
        "sha256": sha256(path),
        "rules_sha256": code_fingerprint(Path(__file__).parent),
    }


def _is_bundled_library(path: Path) -> bool:
    bundle_root = getattr(sys, "_MEIPASS", None)
    return bool(
        getattr(sys, "frozen", False)
        and bundle_root
        and path.name in {"current_positive.msp.gz", "current_negative.msp.gz"}
        and path.resolve().is_relative_to(Path(bundle_root).resolve())
    )


def _library_cache_identity(path: Path) -> str:
    if _is_bundled_library(path):
        return f"bundled:{path.name}:{sha256(path)}"
    return str(path.resolve())


def _library_cache_path(path: Path) -> Path:
    digest = hashlib.sha1(_library_cache_identity(path).encode("utf-8")).hexdigest()[:16]
    cache_root_override = os.environ.get("LIPIDGATE_CACHE_DIR")
    if cache_root_override:
        cache_root = Path(cache_root_override)
    elif os.environ.get("LOCALAPPDATA"):
        cache_root = Path(os.environ["LOCALAPPDATA"]) / "LipidGate" / "Cache"
    else:
        cache_root = Path.home() / ".cache" / "lipidgate"
    return cache_root / "ms2_libraries" / f"library_{digest}.pkl"


def _library_cache_valid(path: Path) -> bool:
    cache_path = _library_cache_path(path)
    meta_path = cache_path.with_suffix(".json")
    if not cache_path.is_file() or not meta_path.is_file():
        return False
    try:
        return json.loads(meta_path.read_text(encoding="utf-8")) == _library_cache_metadata(path)
    except (OSError, ValueError):
        return False


def _load_cached_standard_msp(path: Path) -> List[LibraryRecord] | None:
    cache_path = _library_cache_path(path)
    if not _library_cache_valid(path):
        return None
    try:
        with cache_path.open("rb") as handle:
            records = safe_pickle.load(handle)
    except Exception as exc:
        logging.getLogger(__name__).warning("Cannot read library cache %s: %s", cache_path, exc)
        return None
    return records if isinstance(records, list) else None


def _write_cached_standard_msp(path: Path, records: List[LibraryRecord]) -> None:
    cache_path = _library_cache_path(path)
    meta_path = cache_path.with_suffix(".json")
    temp_path = None
    temp_meta = None
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="wb", dir=cache_path.parent, prefix=cache_path.stem, suffix=".tmp", delete=False) as handle:
            temp_path = Path(handle.name)
            pickle.dump(records, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=cache_path.parent, prefix=cache_path.stem, suffix=".tmp", delete=False) as handle:
            temp_meta = Path(handle.name)
            json.dump(_library_cache_metadata(path), handle)
        temp_path.replace(cache_path)
        temp_meta.replace(meta_path)
    except Exception as exc:
        logging.getLogger(__name__).warning("Cannot write library cache %s: %s", cache_path, exc)
    finally:
        for temporary in (temp_path, temp_meta):
            if temporary is not None:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError as exc:
                    logging.getLogger(__name__).warning("Cannot clean library cache temporary file %s: %s", temporary, exc)


def _bundled_prebuilt_archive(path: Path) -> Path | None:
    if not _is_bundled_library(path):
        return None
    base = Path(sys._MEIPASS) / "libraries" / "ms2" / "prebuilt"
    stem = path.name.removesuffix(".msp.gz")
    archive = base / f"{stem}.pkl.gz"
    metadata = base / f"{stem}.json"
    if not archive.is_file() or not metadata.is_file():
        return None
    expected = {
        "version": LIBRARY_CACHE_VERSION,
        "source_sha256": sha256(path),
        "rules_sha256": code_fingerprint(Path(__file__).parent),
    }
    try:
        if json.loads(metadata.read_text(encoding="utf-8")) != expected:
            return None
        return archive
    except (OSError, ValueError) as exc:
        logging.getLogger(__name__).warning("Cannot read bundled library snapshot %s: %s", archive, exc)
        return None


def _install_bundled_prebuilt_cache(path: Path) -> bool:
    """Expand a validated snapshot to the persistent cache without loading records."""
    archive = _bundled_prebuilt_archive(path)
    if archive is None:
        return False
    cache_path = _library_cache_path(path)
    meta_path = cache_path.with_suffix(".json")
    temporary = None
    temp_meta = None
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="wb", dir=cache_path.parent, prefix=cache_path.stem, suffix=".tmp", delete=False) as destination:
            temporary = Path(destination.name)
            with gzip.open(archive, "rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=cache_path.parent, prefix=cache_path.stem, suffix=".tmp", delete=False) as destination:
            temp_meta = Path(destination.name)
            json.dump(_library_cache_metadata(path), destination)
        temporary.replace(cache_path)
        temp_meta.replace(meta_path)
        return True
    except (OSError, EOFError, ValueError) as exc:
        logging.getLogger(__name__).warning("Cannot install bundled library snapshot %s: %s", archive, exc)
        return False
    finally:
        for leftover in (temporary, temp_meta):
            if leftover is not None:
                leftover.unlink(missing_ok=True)


def _load_bundled_prebuilt_msp(path: Path) -> List[LibraryRecord] | None:
    archive = _bundled_prebuilt_archive(path)
    if archive is None:
        return None
    try:
        with gzip.open(archive, "rb") as handle:
            records = safe_pickle.load(handle)
        return records if isinstance(records, list) else None
    except (OSError, ValueError, pickle.UnpicklingError) as exc:
        logging.getLogger(__name__).warning("Cannot read bundled library snapshot %s: %s", archive, exc)
        return None


EXCEL_COLUMN_ALIASES = {
    "main_class": "main_class",
    "lipid_name": "lipid_name",
    "lipid_chain_name": "lipid_chain_name",
    "化学式": "formula",
    "加合物类型": "adduct",
    "加合物m/z": "precursor_mz",
    "碎片名": "fragment_name",
    "碎片m/z": "fragment_mz",
    "Fragment_Type": "fragment_type",
}


def _matches_special_hg_fragment(fragment_name: str, fragment_mz: float | None, target_mz: float) -> bool:
    if fragment_mz is not None and abs(float(fragment_mz) - target_mz) <= 0.01:
        return True
    text = str(fragment_name or "")
    if not text:
        return False
    return f"{target_mz:.4f}" in text or f"{target_mz:.1f}" in text


def _is_named_negative_headgroup_fragment(
    fragment_type: str,
    fragment_name: str,
    adduct: str,
) -> bool:
    if str(fragment_type or "").strip() not in {"Common", "Candidate_HG", "Diagnostic_HG"}:
        return False
    if not str(adduct or "").strip().endswith("-"):
        return False
    return str(fragment_name or "").strip() in {
        "[C4H11NO4P]-",
        "M-CH3",
        "[M-CH3]-",
    }


def _is_positive_mg_adduct(adduct: str) -> bool:
    return str(adduct or "").strip().endswith("+")


def _is_positive_mg_hg_fragment(fragment_name: str) -> bool:
    name = str(fragment_name or "").strip()
    if name == "[R1C=O-H2O]+":
        return True
    return bool(name.startswith("(R=O)+(") and name.endswith(")"))


def _is_positive_mg_negative_mode_fragment(fragment_name: str) -> bool:
    return str(fragment_name or "").strip().startswith("[RCOO]-")


def _is_class_specific_hg_fragment(
    compound_class: str,
    fragment_type: str,
    fragment_name: str,
    fragment_mz: float | None,
    adduct: str,
) -> bool:
    cls = str(compound_class or "").strip()
    ftype = str(fragment_type or "").strip()
    if ftype not in {"Common", "Candidate_HG"}:
        return False
    if cls == "MG" and _is_positive_mg_adduct(adduct):
        return str(fragment_name or "").strip() == "[M-H2O+H]+"
    if fragment_mz is None or not str(adduct or "").strip().endswith("-"):
        return False
    class_targets = {
        "PG": (152.9933, 171.0064, 209.0221, 227.0326),
        "PEtOH": (181.0271, 181.0280),
        "PMeOH": (167.0109,),
        "DMPE": (168.0411, 168.0431),
    }
    return any(abs(float(fragment_mz) - target_mz) <= 0.02 for target_mz in class_targets.get(cls, ()))


def _normalize_fragment_type(
    compound_class: str,
    fragment_type: str,
    fragment_name: str = "",
    fragment_mz: float | None = None,
    adduct: str = "",
) -> str | None:
    cls = str(compound_class or "").strip()
    ftype = str(fragment_type or "").strip()
    name = str(fragment_name or "").strip()
    if cls.upper() == "PE" and str(adduct or "").strip() == "[M+Na]+" and name in {"[M+Na-43]+", "[M+Na-141]+"}:
        return "Diagnostic_HG"
    if (
        cls.upper() == "PE" and str(adduct or "").strip() == "[M+Na]+"
        and re.match(r"^(?:\(R=O\)\+|RCO\(|\[RCO\]\+)", re.sub(r"\s+", "", name), flags=re.I)
    ):
        return "FA_Frag"
    if (
        cls == "PS"
        and str(adduct or "").strip() == "[M+H]+"
        and re.sub(r"\s+", "", name).upper()
        in {"[M-C3H8O6NP+H]+", "[M-C3H8NO6P+H]+"}
    ):
        return "Diagnostic_HG"
    if _is_named_negative_headgroup_fragment(ftype, name, adduct):
        return "Diagnostic_HG"
    if cls == "MG" and _is_positive_mg_adduct(adduct):
        if _is_positive_mg_negative_mode_fragment(name) and ftype == "Diagnostic_FA":
            return None
        if _is_positive_mg_hg_fragment(name):
            return "Diagnostic_HG"
    if _is_class_specific_hg_fragment(cls, ftype, name, fragment_mz, adduct):
        return "Diagnostic_HG"
    if ftype == "头基诊断碎片" or ("璇婃柇" in ftype and "纰庣墖" in ftype):
        return "Diagnostic_HG"
    if cls in {"Cer1P", "SM", "LSM", "LacCer", "Hex2Cer"} and ftype == "头基诊断碎片":
        return "Diagnostic_HG"
    if cls in {"PI", "LPI", "PIO", "LPIO", "PI-O", "LPI-O"} and name == "[M-C6H13O9P+H]+":
        return "Diagnostic_HG"
    if cls == "NAGly" and _matches_special_hg_fragment(fragment_name, fragment_mz, 76.0393):
        return "Diagnostic_HG"
    if cls == "NAGlySer" and _matches_special_hg_fragment(fragment_name, fragment_mz, 106.0499):
        return "Diagnostic_HG"
    if cls == "NAOrn" and _matches_special_hg_fragment(fragment_name, fragment_mz, 115.0866):
        return "Diagnostic_HG"
    if cls == "CE" and _matches_special_hg_fragment(fragment_name, fragment_mz, 369.3516):
        return "Diagnostic_HG"
    return ftype


def _standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        if column in EXCEL_COLUMN_ALIASES:
            rename_map[column] = EXCEL_COLUMN_ALIASES[column]
    standardized = df.rename(columns=rename_map).copy()
    if "formula" not in standardized.columns and len(standardized.columns) >= 4:
        standardized = standardized.rename(columns={standardized.columns[3]: "formula"})
    if "adduct" not in standardized.columns and len(standardized.columns) >= 5:
        standardized = standardized.rename(columns={standardized.columns[4]: "adduct"})
    if "precursor_mz" not in standardized.columns and len(standardized.columns) >= 6:
        standardized = standardized.rename(columns={standardized.columns[5]: "precursor_mz"})
    if "fragment_name" not in standardized.columns and len(standardized.columns) >= 7:
        standardized = standardized.rename(columns={standardized.columns[6]: "fragment_name"})
    if "fragment_mz" not in standardized.columns and len(standardized.columns) >= 8:
        standardized = standardized.rename(columns={standardized.columns[7]: "fragment_mz"})
    if "fragment_type" not in standardized.columns and len(standardized.columns) >= 9:
        standardized = standardized.rename(columns={standardized.columns[8]: "fragment_type"})
    return standardized


def _required_group_for_fragment(fragment_type: str) -> str | None:
    if fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}:
        return "fah"
    if fragment_type == "Diagnostic_HG":
        return "hg"
    return None


def _ensure_positive_choline_common_fragments(
    compound_class: str,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> List[FragmentRecord]:
    result = list(fragments)
    class_key = str(compound_class or "").strip().upper().replace("-", "")
    if class_key not in POSITIVE_CHOLINE_CLASS_KEYS:
        return result
    if str(adduct or "").strip() != "[M+H]+":
        return result
    for target_mz, name in POSITIVE_CHOLINE_COMMON_FRAGMENTS:
        matching_fragments = [
            fragment
            for fragment in result
            if abs(float(fragment.mz) - target_mz) <= 0.02
        ]
        result = [
            fragment
            for fragment in result
            if abs(float(fragment.mz) - target_mz) > 0.02
        ]
        source_fragment = matching_fragments[0] if matching_fragments else None
        result.append(
            FragmentRecord(
                mz=float(source_fragment.mz) if source_fragment is not None else target_mz,
                intensity=float(source_fragment.intensity) if source_fragment is not None else 100.0,
                name=name,
                fragment_type="Common",
                weight=float(source_fragment.weight) if source_fragment is not None else 1.0,
                required_group=None,
            )
        )
    return result


def _infer_polarity_from_adduct(adduct: str, default: str = "-") -> str:
    text = str(adduct or "").strip()
    if text.endswith("+"):
        return "+"
    if text.endswith("-"):
        return "-"
    return default


def load_excel_directory(directory: str | Path) -> List[LibraryRecord]:
    directory = Path(directory)
    files = sorted(
        file_path
        for file_path in directory.glob("*.xlsx")
        if not file_path.name.startswith("~$")
    )
    records: List[LibraryRecord] = []
    next_id = 0
    for file_path in files:
        df = _standardize_columns(pd.read_excel(file_path))
        required_columns = {
            "main_class",
            "lipid_name",
            "lipid_chain_name",
            "adduct",
            "precursor_mz",
            "fragment_name",
            "fragment_mz",
            "fragment_type",
        }
        missing = required_columns.difference(df.columns)
        if missing:
            raise ValueError(f"{file_path} 缺少列: {sorted(missing)}")
        grouped = df.groupby(["main_class", "lipid_name", "lipid_chain_name", "adduct", "precursor_mz"], dropna=False)
        for (main_class, lipid_name, lipid_chain_name, adduct, precursor_mz), group_df in grouped:
            adduct = str(adduct or "").strip()
            main_class, lipid_name, lipid_chain_name = _canonicalize_sphingoid_base_identity(
                main_class,
                lipid_name,
                lipid_chain_name,
            )
            main_class, lipid_name, lipid_chain_name = canonicalize_glc_gp_identity(
                main_class, lipid_name, lipid_chain_name,
            )
            lipid_name = canonicalize_n_acyl_glycerophospholipid_name(lipid_name, main_class)
            lipid_chain_name = canonicalize_n_acyl_glycerophospholipid_name(
                lipid_chain_name,
                main_class,
            )
            formula = ""
            if "formula" in group_df.columns and not group_df["formula"].isna().all():
                formula = str(group_df["formula"].dropna().iloc[0])
            fragments = []
            for row in group_df.itertuples(index=False):
                fragment_mz = float(getattr(row, "fragment_mz"))
                fragment_name = str(getattr(row, "fragment_name"))
                fragment_type = _normalize_fragment_type(
                    str(main_class),
                    str(getattr(row, "fragment_type")),
                    fragment_name=fragment_name,
                    fragment_mz=fragment_mz,
                    adduct=str(adduct),
                )
                if fragment_type is None:
                    continue
                fragments.append(
                    FragmentRecord(
                        mz=fragment_mz,
                        name=fragment_name,
                        fragment_type=fragment_type,
                        intensity=100.0,
                        weight=1.0,
                        required_group=_required_group_for_fragment(fragment_type),
                    )
                )
            fragments = _ensure_positive_choline_common_fragments(
                str(main_class),
                str(adduct),
                fragments,
            )
            fragments = _normalize_targeted_positive_sphingolipid_fragments(
                str(main_class),
                str(lipid_chain_name),
                float(precursor_mz),
                str(adduct),
                fragments,
            )
            fragments = _normalize_special_aminophospholipid_fragments(
                str(main_class),
                str(lipid_chain_name),
                float(precursor_mz),
                str(adduct),
                fragments,
            )
            fragments = _normalize_tg_o_fragments(
                str(main_class),
                str(lipid_chain_name),
                str(adduct),
                fragments,
            )
            fragments = _normalize_tg_est_fragments(
                str(main_class),
                str(lipid_chain_name),
                float(precursor_mz),
                str(adduct),
                fragments,
            )
            fragments.sort(key=lambda fragment: fragment.mz)
            records.append(
                LibraryRecord(
                    record_id=next_id,
                    compound_class=str(main_class),
                    lipid_name=str(lipid_name),
                    lipid_chain_name=str(lipid_chain_name),
                    precursor_mz=float(precursor_mz),
                    adduct=str(adduct),
                    formula=formula,
                    polarity=_infer_polarity_from_adduct(str(adduct)),
                    fragments=fragments,
                    metadata={"source_file": file_path.name},
                )
            )
            next_id += 1
    return _expand_hydroxy_fa_sphingolipids([
        normalize_protonated_single_chain_record(normalize_nae_record(record)) for record in records
    ])


def _quote_msp_field(value: object) -> str:
    text = str(value if value is not None else "")
    text = text.replace('"', "'")
    return f'"{text}"'


def _encode_fragment_comment(fragment: FragmentRecord) -> str:
    payload = {
        "name": fragment.name,
        "type": fragment.fragment_type,
    }
    inferred_group = _required_group_for_fragment(fragment.fragment_type)
    actual_group = fragment.required_group or ""
    if actual_group and actual_group != (inferred_group or ""):
        payload["required_group"] = actual_group
    rounded_weight = round(float(fragment.weight), 6)
    if abs(rounded_weight - 1.0) > 1e-9:
        payload["weight"] = rounded_weight
    return repr(payload)


def _encode_fragment_fields(fragment: FragmentRecord) -> str:
    return f"{_quote_msp_field(fragment.name)} {_quote_msp_field(fragment.fragment_type)}"


def _decode_fragment_comment(comment: str) -> Dict[str, str]:
    comment = comment.strip()
    if not comment:
        return {}
    try:
        payload = ast.literal_eval(comment)
    except (ValueError, SyntaxError):
        return {"name": comment, "type": "Common", "required_group": "", "weight": 1.0}
    if not isinstance(payload, dict):
        return {"name": comment, "type": "Common", "required_group": "", "weight": 1.0}
    return payload


def _parse_fragment_payload(line: str) -> tuple[float, float, Dict[str, str]] | None:
    pieces = line.split(None, 2)
    if len(pieces) < 2:
        return None
    mz = float(pieces[0])
    intensity = float(pieces[1])
    remainder = pieces[2].strip() if len(pieces) > 2 else ""
    if remainder.startswith('"') and remainder.endswith('"') and len(remainder) >= 2:
        inner_comment = remainder[1:-1]
        if inner_comment.lstrip().startswith("{") or remainder.count('"') == 2:
            payload = _decode_fragment_comment(inner_comment)
            return mz, intensity, payload
    quoted = re.findall(r'"([^"]*)"', remainder)
    if len(quoted) >= 2:
        payload: Dict[str, str] = {
            "name": quoted[0],
            "type": quoted[1] or "Common",
            "required_group": "",
            "weight": 1.0,
        }
        return mz, intensity, payload
    if len(quoted) == 1:
        payload = _decode_fragment_comment(quoted[0])
        return mz, intensity, payload
    payload = _decode_fragment_comment(remainder)
    return mz, intensity, payload


def _decode_record_comment(comment: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    matches = list(re.finditer(r"(?:^|;)\s*(ms1_name|polarity)=", str(comment), flags=re.IGNORECASE))
    for index, match in enumerate(matches):
        key = match.group(1).strip().lower()
        value_start = match.end()
        value_end = matches[index + 1].start() if index + 1 < len(matches) else len(comment)
        value = str(comment)[value_start:value_end]
        if value.startswith(";"):
            value = value[1:]
        fields[key] = value.strip()
    return fields


def _iter_standard_msp_lines(msp_path: Path) -> Iterator[str]:
    if msp_path.name.lower().endswith(".msp.gz"):
        with gzip.open(msp_path, "rt", encoding="utf-8") as handle:
            yield from handle
        return
    with msp_path.open("rt", encoding="utf-8") as handle:
        yield from handle


def write_standard_msp(records: Iterable[LibraryRecord], output_path: str | Path) -> Path:
    output_path = Path(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f"Name: {record.lipid_chain_name}\n")
            handle.write(f"PrecursorMZ: {record.precursor_mz:.4f}\n")
            handle.write(f"PrecursorType: {record.adduct}\n")
            handle.write(f"CompoundClass: {record.compound_class}\n")
            if record.formula:
                handle.write(f"Formula: {record.formula}\n")
            handle.write(f"Comment: MS1_name={record.lipid_name};polarity={record.polarity}\n")
            handle.write(f"Num Peaks: {len(record.fragments)}\n")
            for fragment in record.fragments:
                handle.write(
                    f"{fragment.mz:.4f} {fragment.intensity:.2f} {_encode_fragment_fields(fragment)}\n"
                )
            handle.write("\n")
    return output_path


def convert_excel_directory_to_msp(input_directory: str | Path, output_path: str | Path) -> Path:
    records = load_excel_directory(input_directory)
    return write_standard_msp(records, output_path)


def normalize_imported_record(record: LibraryRecord) -> LibraryRecord:
    """Explicit legacy import policy, separate from raw MSP decoding."""
    from .oxidized_fragment_policy import normalize_oxidized_fragments
    current = {"compoundclass": record.compound_class, "name": record.lipid_chain_name,
               "ms1_name": record.lipid_name, "precursormz": record.precursor_mz,
               "precursortype": record.adduct}
    fragments = []
    for fragment in record.fragments:
        normalized_type = _normalize_fragment_type(
            record.compound_class, fragment.fragment_type, fragment_name=fragment.name,
            fragment_mz=fragment.mz, adduct=record.adduct,
        )
        if normalized_type is not None:
            fragments.append(replace(fragment, fragment_type=normalized_type,
                                     required_group=_required_group_for_fragment(normalized_type) or fragment.required_group))
    normalized_fragments = list(fragments)
    adduct = str(current.get("precursortype", "") or "").strip()
    if current.get("compoundclass", "").strip() == "SL":
        normalized_fragments = annotate_positive_sl_fragments(
            current.get("name", ""),
            float(current["precursormz"]),
            adduct,
            normalized_fragments,
        )
    compound_class, lipid_name, lipid_chain_name = _canonicalize_sphingoid_base_identity(
        current.get("compoundclass", ""),
        current.get("ms1_name", current.get("name", "")),
        current.get("name", ""),
    )
    compound_class, lipid_name, lipid_chain_name = canonicalize_glc_gp_identity(
        compound_class, lipid_name, lipid_chain_name,
    )
    lipid_name = canonicalize_n_acyl_glycerophospholipid_name(lipid_name, compound_class)
    lipid_chain_name = canonicalize_n_acyl_glycerophospholipid_name(
        lipid_chain_name,
        compound_class,
    )
    normalized_fragments = _ensure_positive_choline_common_fragments(
        compound_class,
        adduct,
        normalized_fragments,
    )
    normalized_fragments = _normalize_targeted_positive_sphingolipid_fragments(
        compound_class,
        lipid_chain_name,
        float(current["precursormz"]),
        adduct,
        normalized_fragments,
    )
    normalized_fragments = _normalize_special_aminophospholipid_fragments(
        compound_class,
        lipid_chain_name,
        float(current["precursormz"]),
        adduct,
        normalized_fragments,
    )
    normalized_fragments = _normalize_complex_ceramide_fragments(compound_class,adduct,lipid_chain_name,normalized_fragments)
    normalized_fragments = _normalize_tg_o_fragments(
        compound_class,
        lipid_chain_name,
        adduct,
        normalized_fragments,
    )
    normalized_fragments = _normalize_tg_est_fragments(
        compound_class,
        lipid_chain_name,
        float(current["precursormz"]),
        adduct,
        normalized_fragments,
    )
    normalized_fragments = normalize_oxidized_fragments(
        compound_class, lipid_chain_name, record.precursor_mz, adduct, normalized_fragments,
    )
    normalized_fragments = _curate_positive_lyso_pe_fragments(
        compound_class, record.precursor_mz, adduct, normalized_fragments,
    )
    normalized_fragments = _curate_positive_pe_p_fragments(
        compound_class, lipid_chain_name, adduct, normalized_fragments,
    )
    normalized_fragments = _without_pi_ammonium_acyl_losses(
        compound_class, adduct, normalized_fragments,
    )
    return normalize_positive_ether_lpc(normalize_negative_lpc_record(normalize_negative_lpe_record(normalize_lpc_sodium(normalize_protonated_single_chain_record(normalize_nae_record(annotate_record_chain_labels(replace(record, compound_class=compound_class, lipid_name=lipid_name,
                   lipid_chain_name=lipid_chain_name,
                   fragments=sorted(normalized_fragments, key=lambda fragment: fragment.mz)))))))))


def load_standard_msp(msp_path: str | Path, *, normalize: bool = True) -> List[LibraryRecord]:
    """Decode MSP; normalize=False exposes disk records for migration/audit."""
    records = list(iter_standard_msp_records(msp_path, normalize=normalize))
    return _expand_hydroxy_fa_sphingolipids(records) if normalize else records


def iter_standard_msp_records(msp_path: str | Path, *, normalize: bool = True) -> Iterator[LibraryRecord]:
    """Decode one MSP record at a time, before cohort-wide hydroxy-FA expansion."""
    msp_path = Path(msp_path)
    current: Dict[str, str] = {}
    fragments: List[FragmentRecord] = []
    next_id = 0
    reading_peaks = False

    def flush_current() -> LibraryRecord | None:
        nonlocal current, fragments, next_id, reading_peaks
        if not current:
            return
        if "precursormz" not in current:
            current = {}
            fragments = []
            reading_peaks = False
            return
        record = LibraryRecord(
            record_id=next_id, compound_class=current.get("compoundclass", ""),
            lipid_name=current.get("ms1_name", current.get("name", "")),
            lipid_chain_name=current.get("name", ""), precursor_mz=float(current["precursormz"]),
            adduct=current.get("precursortype", ""), formula=current.get("formula", ""),
            polarity=current.get("polarity", _infer_polarity_from_adduct(current.get("precursortype", ""))),
            fragments=list(fragments),
        )
        record = normalize_imported_record(record) if normalize else record
        next_id += 1
        current = {}
        fragments = []
        reading_peaks = False
        return record

    for raw_line in _iter_standard_msp_lines(msp_path):
        line = raw_line.strip()
        if not line:
            record = flush_current()
            if record is not None:
                yield record
            continue
        if ":" in line and not reading_peaks:
            key, value = line.split(":", 1)
            normalized_key = key.strip().lower()
            value = value.strip()
            if normalized_key == "comment":
                current.update(_decode_record_comment(value))
            else:
                current[normalized_key] = value
            if normalized_key == "num peaks":
                reading_peaks = True
            continue
        if reading_peaks:
            try:
                parsed_fragment = _parse_fragment_payload(line)
            except ValueError:
                if ":" in line:
                    record = flush_current()
                    if record is not None:
                        yield record
                    key, value = line.split(":", 1)
                    normalized_key = key.strip().lower()
                    value = value.strip()
                    if normalized_key == "comment":
                        current.update(_decode_record_comment(value))
                    else:
                        current[normalized_key] = value
                    if normalized_key == "num peaks":
                        reading_peaks = True
                    continue
                continue
            if parsed_fragment is None:
                continue
            mz, intensity, payload = parsed_fragment
            fragments.append(
                FragmentRecord(
                    mz=mz,
                    intensity=intensity,
                    name=str(payload.get("name", "")),
                    fragment_type=str(payload.get("type", "Common")),
                    required_group=str(payload.get("required_group", "")) or None,
                    weight=float(payload.get("weight", 1.0)),
                )
            )
    record = flush_current()
    if record is not None:
        yield record


def load_library(path: str | Path, use_cache: bool = True) -> List[LibraryRecord]:
    path = Path(path)
    if path.is_dir():
        return load_excel_directory(path)
    if path.name.lower().endswith((".msp", ".msp.gz")):
        if use_cache:
            cached_records = _load_cached_standard_msp(path)
            if cached_records is not None:
                return cached_records
            if _install_bundled_prebuilt_cache(path):
                cached_records = _load_cached_standard_msp(path)
                if cached_records is not None:
                    return cached_records
        records = _load_bundled_prebuilt_msp(path) if use_cache else None
        if records is None:
            records = load_standard_msp(path)
        if use_cache:
            _write_cached_standard_msp(path, records)
        return records
    raise ValueError(f"不支持的库路径（需要 .msp 或 .msp.gz）: {path}")
