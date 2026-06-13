"""Input/output helpers and validation for SphinGOlipID workflows."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

TARGET_REQUIRED_COLUMNS = ("target", "下限", "上限", "RT", "Abund")
TARGET_NUMERIC_COLUMNS = ("target", "下限", "上限", "RT", "Abund")
MS1_LIBRARY_REQUIRED_COLUMNS = ("理论值", "classy", "name", "structure")
MS1_LIBRARY_NUMERIC_COLUMNS = ("理论值",)


class InputValidationError(ValueError):
    """Raised when an input table or text export is not usable."""


@dataclass(frozen=True)
class MS2InputPaths:
    """Resolved file paths for one MS2 batch item."""

    text_file: Path
    target_file: Path
    result_file: Path


def require_file(path: str | Path, label: str) -> Path:
    """Return ``path`` as a Path, or raise a clear missing-file error."""

    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Missing required {label}: {resolved}")
    if not resolved.is_file():
        raise InputValidationError(f"Expected {label} to be a file: {resolved}")
    return resolved


def ensure_output_dir(path: str | Path) -> Path:
    """Create and check that an output directory is writable."""

    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    probe = output_dir / ".sphingolipid_write_test"
    try:
        probe.write_text("ok", encoding="utf-8")
    except OSError as exc:
        raise OSError(f"Output directory is not writable: {output_dir}") from exc
    finally:
        if probe.exists():
            probe.unlink()
    return output_dir


def require_columns(data: pd.DataFrame, required_columns: Iterable[str], context: str) -> None:
    """Raise ``InputValidationError`` when a table is missing columns."""

    missing = [col for col in required_columns if col not in data.columns]
    if missing:
        details = "; ".join(f"Missing required column: {col}" for col in missing)
        raise InputValidationError(f"{context}: {details}")


def coerce_numeric_columns(data: pd.DataFrame, columns: Iterable[str], context: str) -> pd.DataFrame:
    """Convert columns to numeric values and fail on invalid non-empty cells."""

    out = data.copy()
    for col in columns:
        if col not in out.columns:
            raise InputValidationError(f"{context}: Missing required column: {col}")
        original = out[col]
        converted = pd.to_numeric(original, errors="coerce")
        invalid = converted.isna() & original.notna() & (original.astype(str).str.strip() != "")
        if invalid.any():
            bad_value = original[invalid].iloc[0]
            raise InputValidationError(f"{context}: Column {col} contains non-numeric value: {bad_value!r}")
        out[col] = converted
    return out


def read_excel_table(
    path: str | Path,
    required_columns: Sequence[str] = (),
    numeric_columns: Sequence[str] = (),
    context: str | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Read an Excel sheet and validate required/numeric columns."""

    table_path = require_file(path, context or "Excel file")
    label = context or str(table_path)
    data = pd.read_excel(table_path, **kwargs)
    require_columns(data, required_columns, label)
    if numeric_columns:
        data = coerce_numeric_columns(data, numeric_columns, label)
    return data


def resolve_ms2_input_paths(config: object, file_index: int) -> MS2InputPaths:
    """Resolve text, precursor, and result paths from either config class."""

    raw_ms2_dir = Path(_get_required_attr(config, "raw_ms2_dir", "input_dir"))
    precursor_dir = Path(getattr(config, "precursor_dir", raw_ms2_dir) or raw_ms2_dir)
    output_dir = Path(getattr(config, "output_dir"))
    text_pattern = getattr(config, "text_pattern")
    target_pattern = getattr(config, "target_pattern")
    result_pattern = getattr(config, "result_pattern")
    return MS2InputPaths(
        text_file=raw_ms2_dir / text_pattern.format(i=file_index),
        target_file=precursor_dir / target_pattern.format(i=file_index),
        result_file=output_dir / result_pattern.format(i=file_index),
    )


def get_ms1_library_path(config: object) -> Path:
    """Return the MS1 library path from either the unified or legacy config."""

    return Path(_get_required_attr(config, "ms1_library_path", "ms1_db_path"))


def validate_ms2_file_inputs(config: object, file_index: int) -> tuple[MS2InputPaths, pd.DataFrame]:
    """Validate all inputs needed for one MS2 file and return target rows."""

    paths = resolve_ms2_input_paths(config, file_index)
    require_file(paths.text_file, f"MS/MS txt file for index {file_index}")
    target_df = read_excel_table(
        paths.target_file,
        required_columns=TARGET_REQUIRED_COLUMNS,
        numeric_columns=TARGET_NUMERIC_COLUMNS,
        context=f"precursor target Excel for index {file_index}",
    )
    read_excel_table(
        get_ms1_library_path(config),
        required_columns=MS1_LIBRARY_REQUIRED_COLUMNS,
        numeric_columns=MS1_LIBRARY_NUMERIC_COLUMNS,
        context="MS1 theoretical library",
    )
    ensure_output_dir(getattr(config, "output_dir"))
    return paths, target_df


def validate_input_paths(config: object, file_indices: Sequence[int] | None = None) -> list[MS2InputPaths]:
    """Validate configured MS2 inputs for all requested file indices."""

    indices = tuple(file_indices or getattr(config, "file_indices"))
    validated: list[MS2InputPaths] = []
    for file_index in indices:
        paths, _ = validate_ms2_file_inputs(config, int(file_index))
        validated.append(paths)
    return validated


def validate_required_columns(data: pd.DataFrame, required_columns: Iterable[str], context: str) -> None:
    """Compatibility wrapper with the function name used in the v0.3 plan."""

    require_columns(data, required_columns, context)


def read_precursor_table(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read and validate a precursor target Excel table."""

    return read_excel_table(
        path,
        required_columns=TARGET_REQUIRED_COLUMNS,
        numeric_columns=TARGET_NUMERIC_COLUMNS,
        context="precursor target Excel",
        **kwargs,
    )


def read_ms1_library(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read and validate the theoretical MS1 library."""

    return read_excel_table(
        path,
        required_columns=MS1_LIBRARY_REQUIRED_COLUMNS,
        numeric_columns=MS1_LIBRARY_NUMERIC_COLUMNS,
        context="MS1 theoretical library",
        **kwargs,
    )


def read_ms2_txt_standardized(text_file: str | Path, encoding: str = "GBK") -> pd.DataFrame:
    """Compatibility wrapper for standardized MS2 txt parsing."""

    return parse_ms2_text_to_dataframe(text_file, encoding=encoding)


def read_ms2_spectra_text(text_file: str | Path, encoding: str = "GBK") -> list[str]:
    """Read a legacy Agilent-style MS/MS txt export as spectrum blocks."""

    path = require_file(text_file, "MS/MS txt file")
    data = path.read_text(encoding=encoding, errors="ignore")
    spectra = data.split("spectrum:")[1:]
    if not spectra:
        raise InputValidationError(f"MS/MS txt file contains no spectrum blocks: {path}")
    return spectra


def parse_ms2_text_to_dataframe(text_file: str | Path, encoding: str = "GBK") -> pd.DataFrame:
    """Convert a supported MS/MS txt export to a normalized fragment table.

    Returned columns follow the planned backend contract:
    ``file_name``, ``scan_id``, ``precursor_mz``, ``rt``, ``fragment_mz``,
    and ``fragment_intensity``.
    """

    path = Path(text_file)
    rows: list[dict[str, object]] = []
    for block in read_ms2_spectra_text(path, encoding=encoding):
        target_match = re.search(r"target m/z,\s*([^,]+),\s*m/z", block)
        index_match = re.search(r"index:(.*?)\n", block)
        rt_match = re.search(r"scan start time,\s*([^,]+),", block)
        arrays = block.split("binaryDataArray:")[1:]
        if not target_match or not index_match or len(arrays) < 2:
            continue
        mz_values = _parse_binary_values(arrays[0])
        intensity_values = _parse_binary_values(arrays[1])
        precursor_mz = _to_float(target_match.group(1))
        rt = _to_float(rt_match.group(1)) if rt_match else None
        for mz, intensity in zip(mz_values, intensity_values):
            rows.append(
                {
                    "file_name": path.name,
                    "scan_id": index_match.group(1).strip(),
                    "precursor_mz": precursor_mz,
                    "rt": rt,
                    "fragment_mz": mz,
                    "fragment_intensity": intensity,
                }
            )

    return pd.DataFrame(
        rows,
        columns=["file_name", "scan_id", "precursor_mz", "rt", "fragment_mz", "fragment_intensity"],
    )


def _parse_binary_values(block: str) -> list[float]:
    match = re.search(r"binary:\s*(.*)", block)
    if not match:
        return []
    tokens = match.group(1).split()[1:]
    if tokens and tokens[-1] == "":
        tokens = tokens[:-1]
    values: list[float] = []
    for token in tokens:
        value = _to_float(token)
        if value is not None:
            values.append(value)
    return values


def _to_float(value: object) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _get_required_attr(config: object, *names: str) -> object:
    for name in names:
        if hasattr(config, name):
            return getattr(config, name)
    joined = " or ".join(names)
    raise AttributeError(f"Config object must define {joined}")
