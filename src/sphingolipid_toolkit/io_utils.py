"""Input/output helpers and validation for SphinGOlipID workflows."""

from __future__ import annotations

import re
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from tempfile import NamedTemporaryFile
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
    try:
        with NamedTemporaryFile(dir=output_dir, prefix=".sphingolipid_write_test_", mode="w") as probe:
            probe.write("ok")
    except OSError as exc:
        raise OSError(f"Output directory is not writable: {output_dir}") from exc
    return output_dir


def write_excel_atomic(data: pd.DataFrame, path: str | Path, **kwargs) -> Path:
    """Replace an Excel result only after a complete workbook has been written."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=destination.parent, prefix=".sphingolipid_result_", suffix=destination.suffix, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        data.to_excel(temporary, **kwargs)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


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
        invalid = (
            converted.isna() & original.notna() & (original.astype(str).str.strip() != "")
        ) | converted.isin([float("inf"), -float("inf")])
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


def protected_library_paths(config: object) -> set[Path]:
    """Protect selected reference inputs and the bundled scientific assets."""
    from importlib.resources import files
    protected=set()
    for name in ('ms1_library_path','ms1_db_path','ms2_library_path'):
        selected=getattr(config,name,None)
        if selected is not None:protected.add(Path(selected).resolve())
    bundled=Path(str(files('sphingolipid_toolkit').joinpath('data','ms2')))
    protected.update((bundled/f'sphingolipid_{mode}{suffix}').resolve()
        for mode in ('positive','negative') for suffix in ('.sqlite','.msp.gz'))
    return protected


def validate_output_paths(config: object, file_indices: Sequence[int]) -> None:
    """Check every batch output before cleanup can overwrite another input."""

    paths = [resolve_ms2_input_paths(config, i) for i in file_indices]
    protected=protected_library_paths(config)
    protected.update(path.resolve() for item in paths for path in (item.text_file, item.target_file))
    output_dir = Path(config.output_dir)
    intermediate = output_dir / config.intermediate_dir_name
    outputs = [item.result_file for item in paths]
    outputs.extend(output_dir / name for name in (config.final_result_name, config.log_file_name))
    outputs.extend(intermediate / name for name in ("Isomer.xlsx", "out_DB.csv", "out-query.csv", "目标母离子.xlsx"))
    outputs.extend(intermediate / f"{prefix}-{i}.xlsx" for i in file_indices for prefix in ("actual_ms2_data", "DB_sheet"))
    if any(path.resolve() in protected for path in outputs):
        raise InputValidationError("An output file would overwrite a selected input file; choose another output folder or result name")
    resolved_outputs = [path.resolve() for path in outputs]
    if len(set(resolved_outputs)) != len(resolved_outputs):
        raise InputValidationError("Result, log and intermediate files must have distinct output paths")


def validate_ms2_file_inputs(config: object, file_index: int) -> tuple[MS2InputPaths, pd.DataFrame]:
    """Validate all inputs needed for one MS2 file and return target rows."""

    paths = resolve_ms2_input_paths(config, file_index)
    require_file(paths.text_file, f"MS/MS txt file for index {file_index}")
    target_df = read_precursor_table(paths.target_file)
    validate_output_paths(config, (file_index,))
    ensure_output_dir(getattr(config, "output_dir"))
    return paths, target_df


def validate_input_paths(config: object, file_indices: Sequence[int] | None = None) -> list[MS2InputPaths]:
    """Validate configured MS2 inputs for all requested file indices."""

    indices = tuple(getattr(config, "file_indices") if file_indices is None else file_indices)
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

    data = read_excel_table(
        path,
        required_columns=TARGET_REQUIRED_COLUMNS,
        numeric_columns=TARGET_NUMERIC_COLUMNS,
        context="precursor target Excel",
        **kwargs,
    )
    if data[list(TARGET_NUMERIC_COLUMNS)].isna().any().any():
        raise InputValidationError("precursor target Excel: required numeric cells must not be empty")
    if (data["下限"] > data["上限"]).any():
        raise InputValidationError("precursor target Excel: lower m/z limit must not exceed upper limit")
    return data


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
    return parse_ms2_spectrum_blocks(read_ms2_spectra_text(path, encoding=encoding), file_name=path.name)


def parse_ms2_spectrum_blocks(blocks: Iterable[str], file_name: str = "") -> pd.DataFrame:
    """Parse fragment arrays once for both the normalized API and legacy pipeline."""

    rows: list[dict[str, object]] = []
    for block in blocks:
        target_match = re.search(r"target m/z,\s*([^,]+),\s*m/z", block)
        index_match = re.search(r"index:(.*?)\n", block)
        rt_match = re.search(r"scan start time,\s*([^,]+),", block)
        arrays = block.split("binaryDataArray:")[1:]
        if not target_match:
            continue
        if not index_match or len(arrays) < 2:
            raise InputValidationError("MS/MS spectrum is missing a scan index or fragment arrays")
        mz_values = _parse_binary_values(arrays[0])
        intensity_values = _parse_binary_values(arrays[1])
        if len(mz_values) != len(intensity_values):
            raise InputValidationError(f"MS/MS scan {index_match.group(1).strip()}: m/z and intensity arrays have different lengths")
        precursor_mz = _to_float(target_match.group(1))
        if precursor_mz is None:
            raise InputValidationError("MS/MS spectrum contains an invalid precursor m/z")
        rt = _to_float(rt_match.group(1)) if rt_match else None
        if rt_match and rt is None:
            raise InputValidationError("MS/MS spectrum contains an invalid scan start time")
        for mz, intensity in zip(mz_values, intensity_values):
            rows.append(
                {
                    "file_name": file_name,
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
        raise InputValidationError("MS/MS fragment array is missing binary values")
    tokens = match.group(1).split()[1:]
    values: list[float] = []
    for token in tokens:
        value = _to_float(token)
        if value is None:
            raise InputValidationError(f"MS/MS fragment array contains an invalid numeric value: {token!r}")
        values.append(value)
    return values


def _to_float(value: object) -> float | None:
    try:
        converted = float(str(value).strip())
        return converted if isfinite(converted) else None
    except (TypeError, ValueError):
        return None


def _get_required_attr(config: object, *names: str) -> object:
    for name in names:
        if hasattr(config, name):
            return getattr(config, name)
    joined = " or ".join(names)
    raise AttributeError(f"Config object must define {joined}")
