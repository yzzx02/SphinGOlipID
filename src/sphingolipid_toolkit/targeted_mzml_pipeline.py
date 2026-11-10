"""Targeted mzML workflow that annotates MS2 results back to MS1 features."""

from __future__ import annotations

import logging
import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable, Iterator, Sequence

import numpy as np
import pandas as pd
from pyteomics import mzml

from . import ms2_legacy_core as core
from .config import TARGETED_FRAGMENT_PPM, validate_matching_parameters
from .fragment_assignment import select_one_to_one_edges
from .glycan_encoding import parse_glycan_encoding
from .io_utils import (
    InputValidationError, MS1_LIBRARY_REQUIRED_COLUMNS, coerce_numeric_columns,
    ensure_output_dir, read_ms1_library, require_file, write_excel_atomic,
)
from .logging_utils import attach_file_handler, configure_logging, get_logger
from .rt_validation import add_series_columns, fit_ransac_models, predict_rt
from .scoring import add_standard_score_columns


@dataclass(frozen=True)
class TargetedMzMLConfig:
    """Configuration for annotating targeted mzML MS2 data to MS1 features."""

    data_dir: Path
    ms1_library_path: Path
    output_dir: Path
    targetlist_dir: Path | None = None
    mzml_dir: Path | None = None
    fragment_ppm: float = TARGETED_FRAGMENT_PPM
    min_fragment_intensity: float = 20.0
    min_matched_fragments: int = 2
    min_match_score: float = 0.35
    top_n: int = 3
    ms1_candidate_ppm: float = 10.0
    feature_mz_ppm: float = 10.0
    feature_rt_tolerance_min: float = 0.5
    feature_merge_ppm: float = 10.0
    feature_merge_mz_tolerance_da: float = 0.005
    feature_merge_rt_tolerance_min: float = 0.05
    rt_fit_residual_threshold: float = 0.1
    rt_fit_min_points: int = 3
    log_file_name: str = "targeted_mzml_run_log.txt"
    raw_result_name: str = "targeted_mzml_raw_identifications.xlsx"
    fitted_result_name: str = "targeted_mzml_rt_fitted_identifications.xlsx"
    feature_table_name: str = "combined_ms1_feature_table.xlsx"
    raw_targetlist_name: str = "combined_targetlist_raw.xlsx"
    backend: str = "lipidgate"
    ms2_library_path: Path | None = None
    fragment_da: float | None = None
    min_total_score: float = 50.0

    def __post_init__(self) -> None:
        validate_matching_parameters(self)
        if self.backend not in {"lipidgate", "legacy"}:
            raise ValueError("backend must be lipidgate or legacy")
        if self.backend == "lipidgate":
            from .engine import config_from_workflow
            config_from_workflow(self)
        for name in ("ms1_candidate_ppm", "feature_mz_ppm", "feature_rt_tolerance_min", "feature_merge_ppm",
                     "feature_merge_mz_tolerance_da", "feature_merge_rt_tolerance_min", "rt_fit_residual_threshold"):
            value = getattr(self, name)
            if not np.isfinite(value) or value < 0 or (name == "rt_fit_residual_threshold" and value == 0):
                raise ValueError(f"Invalid {name}: {value}")
        if not isinstance(self.rt_fit_min_points, int) or isinstance(self.rt_fit_min_points, bool) or self.rt_fit_min_points < 2:
            raise ValueError("rt_fit_min_points must be an integer of at least 2")
        object.__setattr__(self, "data_dir", Path(self.data_dir))
        if self.ms1_library_path is not None:
            object.__setattr__(self, "ms1_library_path", Path(self.ms1_library_path))
        elif self.backend == 'legacy':
            raise ValueError('The legacy backend requires an MS1 library')
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        if self.targetlist_dir is None:
            object.__setattr__(self, "targetlist_dir", Path(self.data_dir) / "targetlist")
        else:
            object.__setattr__(self, "targetlist_dir", Path(self.targetlist_dir))
        if self.mzml_dir is None:
            object.__setattr__(self, "mzml_dir", Path(self.data_dir))
        else:
            object.__setattr__(self, "mzml_dir", Path(self.mzml_dir))

    @property
    def fragment_tolerance_fraction(self) -> float:
        return float(self.fragment_ppm) / 1_000_000


def _ppm_to_da(mz: float, ppm: float) -> float:
    return abs(float(mz)) * float(ppm) / 1_000_000


def _ppm_error(observed: float | pd.Series, reference: float) -> float | pd.Series:
    denominator = max(abs(float(reference)), 1e-12)
    if isinstance(observed, pd.Series):
        return (observed.astype(float) - float(reference)).abs() / denominator * 1_000_000
    return abs(float(observed) - float(reference)) / denominator * 1_000_000


@dataclass(frozen=True)
class MS2Spectrum:
    mzml_file: Path
    scan_id: str
    precursor_mz: float
    ms2_rt: float
    collision_energy: float | None
    fragment_mz: np.ndarray
    fragment_intensity: np.ndarray

    def __post_init__(self) -> None:
        masses = np.asarray(self.fragment_mz, dtype=float)
        intensities = np.asarray(self.fragment_intensity, dtype=float)
        if masses.ndim != 1 or intensities.ndim != 1 or len(masses) != len(intensities):
            raise InputValidationError(f"MS2 scan {self.scan_id}: fragment arrays must be one-dimensional with equal lengths")
        if not np.isfinite(self.precursor_mz) or self.precursor_mz <= 0 or not np.isfinite(self.ms2_rt):
            raise InputValidationError(f"MS2 scan {self.scan_id}: precursor m/z and RT must be finite, and precursor m/z positive")
        object.__setattr__(self, "mzml_file", Path(self.mzml_file))
        object.__setattr__(self, "fragment_mz", masses)
        object.__setattr__(self, "fragment_intensity", intensities)


def _find_input_files(directory: str | Path, suffix: str) -> list[Path]:
    root = Path(directory)
    if not root.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {root}")
    return sorted(path for path in root.rglob("*") if path.is_file() and path.suffix.lower() == suffix.lower())


def _validate_targeted_paths(config: TargetedMzMLConfig) -> None:
    from .engine import library_path
    library = require_file(config.ms1_library_path if config.backend == 'legacy' else
        config.ms2_library_path or library_path('+'), "MS2 library" if config.backend!='legacy' else "MS1 theoretical library")
    targetlists = _find_input_files(config.targetlist_dir, ".csv")
    spectra = _find_input_files(config.mzml_dir, ".mzml")
    if not targetlists or not spectra:
        raise InputValidationError("Targeted workflow requires targetlist CSV files and mzML files")
    inputs = {path.resolve() for path in [library, *targetlists, *spectra]}
    outputs = [(config.output_dir / getattr(config, name)).resolve() for name in
               ("log_file_name", "raw_result_name", "fitted_result_name", "feature_table_name", "raw_targetlist_name")]
    if any(path in inputs for path in outputs):
        raise InputValidationError("An output file would overwrite a targeted input file; choose another output folder or result name")
    if len(set(outputs)) != len(outputs):
        raise InputValidationError("Targeted result and log files must have distinct output paths")


def run_targeted_mzml_batch(
    config: TargetedMzMLConfig,
    logger: logging.Logger | None = None,
) -> dict[str, Path]:
    """Run the complete targeted mzML annotation workflow."""

    _validate_targeted_paths(config)
    logger = logger or get_logger("sphingolipid_toolkit.targeted_mzml")
    output_dir = ensure_output_dir(config.output_dir)
    handler = attach_file_handler(output_dir / config.log_file_name, logger=logger)
    try:
        logger.info("Effective fragment_ppm=%s (targeted mzML workflow config)", config.fragment_ppm)
        return _run_targeted_mzml_batch(config, logger)
    finally:
        logger.removeHandler(handler)
        handler.close()


def _run_targeted_mzml_batch(config: TargetedMzMLConfig, logger: logging.Logger) -> dict[str, Path]:
    output_dir = config.output_dir

    logger.info("Building combined MS1 feature table from targetlists")
    targetlist_raw = read_targetlist_directory(config.targetlist_dir)
    if targetlist_raw.empty:
        raise InputValidationError("No enabled, usable MS1 targets were found in the targetlist CSV files")
    features = build_ms1_feature_table(
        targetlist_raw,
        mz_tolerance_da=config.feature_merge_mz_tolerance_da,
        mz_tolerance_ppm=config.feature_merge_ppm,
        rt_tolerance_min=config.feature_merge_rt_tolerance_min,
    )
    logger.info("Loading MS1 theoretical library")
    library = (_read_legacy_or_normalized_ms1_library(config.ms1_library_path)
               if config.backend == "legacy" else pd.DataFrame())
    write_excel_atomic(targetlist_raw, output_dir / config.raw_targetlist_name, index=False)
    write_excel_atomic(features, output_dir / config.feature_table_name, index=False)

    logger.info("Matching targeted mzML files")
    raw_results = annotate_mzml_directory(features, library, config, logger=logger)
    from .result_writer import write_identifications
    write_identifications(raw_results, output_dir / config.raw_result_name)

    logger.info("Running RT fitting filter with top 1-%s candidates", config.top_n)
    fitted_results = apply_rt_fit_filter(raw_results, config)
    write_identifications(fitted_results, output_dir / config.fitted_result_name)

    logger.info("Targeted mzML workflow completed")
    return {
        "targetlist_raw": output_dir / config.raw_targetlist_name,
        "feature_table": output_dir / config.feature_table_name,
        "raw_identifications": output_dir / config.raw_result_name,
        "rt_fitted_identifications": output_dir / config.fitted_result_name,
        "log": output_dir / config.log_file_name,
    }


def read_targetlist_directory(targetlist_dir: str | Path) -> pd.DataFrame:
    """Read all Agilent targetlist CSV files under a directory."""

    root = Path(targetlist_dir)
    rows: list[pd.DataFrame] = []
    for csv_file in _find_input_files(root, ".csv"):
        frame = read_targetlist_csv(csv_file)
        if frame.empty:
            continue
        frame["source_targetlist_file"] = str(csv_file)
        frame["targetlist_group"] = csv_file.parent.name
        frame["targetlist_name"] = csv_file.stem
        rows.append(frame)
    if not rows:
        return pd.DataFrame(
            columns=[
                "target_mz",
                "ms1_rt",
                "delta_rt",
                "targetlist_collision_energy",
                "source_targetlist_file",
                "targetlist_group",
                "targetlist_name",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def read_targetlist_csv(path: str | Path) -> pd.DataFrame:
    """Read one Agilent TargetedMSMSTable CSV and normalize columns."""

    csv_path = Path(path)
    first_line = ""
    leading_rows = 0
    with csv_path.open(encoding="utf-8-sig") as handle:
        for leading_rows, line in enumerate(handle):
            if line.strip():
                first_line = line
                break
    if not first_line:
        raise InputValidationError(f"Empty targetlist CSV file: {csv_path}")
    skiprows = leading_rows + int(first_line.strip().startswith("TargetedMSMSTable"))
    data = pd.read_csv(csv_path, skiprows=skiprows, encoding="utf-8-sig")
    data.columns = data.columns.str.strip()
    required = {"Prec. m/z", "Ret. Time (min)"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{csv_path}: missing targetlist column(s): {sorted(missing)}")
    if "On" in data.columns:
        data = data[data["On"].astype(str).str.strip().str.lower().isin({"true", "1", "1.0", "yes"})].copy()
    out = pd.DataFrame(
        {
            "target_mz": pd.to_numeric(data["Prec. m/z"], errors="coerce"),
            "ms1_rt": pd.to_numeric(data["Ret. Time (min)"], errors="coerce"),
            "delta_rt": pd.to_numeric(data.get("Delta Ret. Time (min)"), errors="coerce"),
            "targetlist_collision_energy": pd.to_numeric(data.get("Collision Energy"), errors="coerce"),
        }
    )
    return out.dropna(subset=["target_mz", "ms1_rt"]).reset_index(drop=True)


def build_ms1_feature_table(
    targetlist_raw: pd.DataFrame,
    mz_tolerance_da: float = 0.005,
    mz_tolerance_ppm: float | None = None,
    rt_tolerance_min: float = 0.05,
) -> pd.DataFrame:
    """Merge repeated targetlist entries into a feature-level table."""

    if targetlist_raw.empty:
        return pd.DataFrame(columns=["feature_id", "feature_mz", "ms1_rt"])

    data = targetlist_raw.sort_values(["ms1_rt", "target_mz"]).reset_index(drop=True)
    clusters: list[list[int]] = []
    cluster_centers: list[tuple[float, float]] = []
    for idx, row in data.iterrows():
        mz_value = float(row["target_mz"])
        rt_value = float(row["ms1_rt"])
        assigned = False
        for cluster_idx, (center_mz, center_rt) in enumerate(cluster_centers):
            if mz_tolerance_ppm is None:
                mz_match = abs(mz_value - center_mz) <= mz_tolerance_da
            else:
                mz_match = _ppm_error(mz_value, center_mz) <= mz_tolerance_ppm
            if mz_match and abs(rt_value - center_rt) <= rt_tolerance_min:
                clusters[cluster_idx].append(idx)
                group = data.loc[clusters[cluster_idx]]
                cluster_centers[cluster_idx] = (float(group["target_mz"].mean()), float(group["ms1_rt"].mean()))
                assigned = True
                break
        if not assigned:
            clusters.append([idx])
            cluster_centers.append((mz_value, rt_value))

    rows: list[dict[str, object]] = []
    for number, indices in enumerate(clusters, start=1):
        group = data.loc[indices]
        rows.append(
            {
                "feature_id": f"F{number:06d}",
                "feature_mz": float(group["target_mz"].mean()),
                "ms1_rt": float(group["ms1_rt"].mean()),
                "targetlist_row_count": int(len(group)),
                "targetlist_groups": ";".join(sorted(set(group.get("targetlist_group", pd.Series(dtype=str)).dropna().astype(str)))),
                "targetlist_files": ";".join(sorted(set(group.get("source_targetlist_file", pd.Series(dtype=str)).dropna().astype(str)))),
                "mean_delta_rt": float(group["delta_rt"].dropna().mean()) if group["delta_rt"].notna().any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


def annotate_mzml_directory(
    features: pd.DataFrame,
    library: pd.DataFrame,
    config: TargetedMzMLConfig,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Annotate every MS2 spectrum in all mzML files and return long-form rows."""

    logger = logger or get_logger("sphingolipid_toolkit.targeted_mzml")
    if config.backend == "lipidgate":
        from .reference_pipeline import annotate_targeted_directory
        return annotate_targeted_directory(features, config, logger)
    mzml_files = _find_input_files(config.mzml_dir, ".mzml")
    rows: list[pd.DataFrame] = []
    with TemporaryDirectory(prefix="sphingolipid_targeted_") as temp_dir:
        fragment_cache: dict[str, pd.DataFrame] = {}
        for mzml_file in mzml_files:
            logger.info("Reading mzML file: %s", mzml_file)
            for spectrum in iter_ms2_spectra(mzml_file):
                feature = find_matching_feature(
                    features,
                    precursor_mz=spectrum.precursor_mz,
                    ms2_rt=spectrum.ms2_rt,
                    mz_tolerance_ppm=config.feature_mz_ppm,
                    rt_tolerance_min=config.feature_rt_tolerance_min,
                )
                if feature is None:
                    continue
                annotated = annotate_one_spectrum(
                    spectrum,
                    feature,
                    library,
                    config,
                    temp_dir=Path(temp_dir),
                    fragment_cache=fragment_cache,
                    logger=logger,
                )
                if not annotated.empty:
                    rows.append(annotated)
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not result.empty:
        result.insert(0, "identification_id", [f"ID{i:08d}" for i in range(1, len(result) + 1)])
    return result


def iter_ms2_spectra(mzml_file: str | Path) -> Iterator[MS2Spectrum]:
    """Yield MS2 spectra from one mzML file with precursor and RT metadata."""

    path = Path(mzml_file)
    with mzml.MzML(str(path), decode_binary=True) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != 2:
                continue
            precursor = _extract_precursor(spectrum)
            ms2_rt = _extract_scan_rt(spectrum)
            if precursor is None or ms2_rt is None:
                continue
            mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
            intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
            yield MS2Spectrum(
                mzml_file=path,
                scan_id=str(spectrum.get("id", "")),
                precursor_mz=float(precursor["precursor_mz"]),
                ms2_rt=float(ms2_rt),
                collision_energy=precursor.get("collision_energy"),
                fragment_mz=mz_array,
                fragment_intensity=intensity_array,
            )


def find_matching_feature(
    features: pd.DataFrame,
    precursor_mz: float,
    ms2_rt: float,
    mz_tolerance_ppm: float,
    rt_tolerance_min: float,
) -> pd.Series | None:
    """Find the closest MS1 feature for a targeted MS2 spectrum."""

    if features.empty:
        return None
    feature_mz = features["feature_mz"].astype(float)
    mz_error_ppm = (feature_mz - float(precursor_mz)).abs() / max(float(precursor_mz), 1e-12) * 1_000_000
    rt_error = (features["ms1_rt"].astype(float) - float(ms2_rt)).abs()
    candidates = features[(mz_error_ppm <= mz_tolerance_ppm) & (rt_error <= rt_tolerance_min)].copy()
    if candidates.empty:
        return None
    candidates["feature_mz_error_ppm"] = mz_error_ppm.loc[candidates.index]
    candidates["feature_rt_error_min"] = rt_error.loc[candidates.index]
    candidates["_distance"] = (
        candidates["feature_mz_error_ppm"] / max(mz_tolerance_ppm, 1e-12)
        + (candidates["ms1_rt"].astype(float) - float(ms2_rt)).abs() / max(rt_tolerance_min, 1e-12)
    )
    return candidates.sort_values("_distance").iloc[0].drop(labels=["_distance"])


def annotate_one_spectrum(
    spectrum: MS2Spectrum,
    feature: pd.Series,
    library: pd.DataFrame,
    config: TargetedMzMLConfig,
    temp_dir: Path,
    fragment_cache: dict[str, pd.DataFrame] | None = None,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Annotate one MS2 spectrum against legacy-generated candidate fragments."""

    logger = logger or get_logger("sphingolipid_toolkit.targeted_mzml")
    if config.backend == "lipidgate":
        from .engine import SphingolipidIdentifier, config_from_workflow
        from .reference_pipeline import annotate_targeted_spectrum
        with SphingolipidIdentifier(config_from_workflow(config)) as identifier:
            return annotate_targeted_spectrum(spectrum, feature, config, identifier)
    feature_id = str(feature["feature_id"])
    logger.debug("Effective fragment_ppm=%s for scan %s", config.fragment_ppm, spectrum.scan_id)
    candidates = search_ms1_candidates(
        library,
        observed_mz=float(feature["feature_mz"]),
        tolerance_ppm=config.ms1_candidate_ppm,
        ms1_rt=float(feature["ms1_rt"]),
    )
    if candidates.empty:
        return pd.DataFrame()

    cache = fragment_cache if fragment_cache is not None else {}
    if feature_id not in cache:
        cache[feature_id] = generate_legacy_fragments_for_candidates(
            candidates,
            observed_mz=float(feature["feature_mz"]),
            ms1_rt=float(feature["ms1_rt"]),
            abundance=float(feature.get("targetlist_row_count", 1)),
            temp_dir=temp_dir,
            logger=logger,
        )
    theoretical = cache[feature_id]
    if theoretical.empty:
        return pd.DataFrame()

    matched = match_observed_to_theoretical(
        spectrum,
        theoretical,
        fragment_tolerance_fraction=config.fragment_tolerance_fraction,
        min_fragment_intensity=config.min_fragment_intensity,
    )
    result = finalize_spectrum_matches(
        matched,
        min_matched_fragments=config.min_matched_fragments,
        min_match_score=config.min_match_score,
        top_n=config.top_n,
    )
    if result.empty:
        return result

    result = add_standard_score_columns(result)
    lipid_col = "lipid_name" if "lipid_name" in result.columns else "注释"
    candidate_errors = candidates.groupby("name")["ms1_library_mz_error_ppm"].min()
    errors = result["母离子"].map(candidate_errors) if "母离子" in result else result[lipid_col].map(candidate_errors)
    result["ms1_library_mz_error_ppm"] = errors.fillna(result[lipid_col].map(candidate_errors))
    for col, value in feature.items():
        result[col] = value
    result["ms2_rt"] = spectrum.ms2_rt
    result["ms2_precursor_mz"] = spectrum.precursor_mz
    result["ms2_scan_id"] = spectrum.scan_id
    result["ms2_file"] = str(spectrum.mzml_file)
    result["ms2_file_name"] = spectrum.mzml_file.name
    result["collision_energy"] = spectrum.collision_energy
    result["ms1_rt"] = float(feature["ms1_rt"])
    result["feature_mz"] = float(feature["feature_mz"])
    return result


def search_ms1_candidates(
    library: pd.DataFrame,
    observed_mz: float,
    tolerance_ppm: float,
    ms1_rt: float,
) -> pd.DataFrame:
    """Search theoretical MS1 library candidates for one feature."""

    tolerance_da = _ppm_to_da(observed_mz, tolerance_ppm)
    low = float(observed_mz) - tolerance_da
    high = float(observed_mz) + tolerance_da
    candidates = library[(library["理论值"] >= low) & (library["理论值"] <= high)].copy()
    if candidates.empty:
        return candidates
    candidates["ms1_library_mz_error_ppm"] = _ppm_error(candidates["理论值"].astype(float), observed_mz)
    candidates["target"] = float(observed_mz)
    candidates["RT"] = float(ms1_rt)
    candidates["Abund"] = 1.0
    return candidates


@core.synchronized_legacy_runtime
def generate_legacy_fragments_for_candidates(
    candidates: pd.DataFrame,
    observed_mz: float,
    ms1_rt: float,
    abundance: float,
    temp_dir: Path,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Generate theoretical fragments by calling the legacy fragment rules."""

    logger = logger or get_logger("sphingolipid_toolkit.targeted_mzml")
    temp_dir.mkdir(parents=True, exist_ok=True)
    out_db = temp_dir / "out_DB.csv"
    if out_db.exists():
        out_db.unlink()
    old_folder = core.folder
    core.folder = str(temp_dir) + os.sep
    try:
        for _, row in candidates.iterrows():
            core.frag_id1 = ["[C2H5NO+H]+"]
            core.frag_mz1 = [60.044]
            structure = row.get("structure", "")
            classy = row["classy"]
            if isinstance(structure, str) and structure and not re.match(r"^[a-zA-Z]", str(classy)):
                # A malformed encoding is an input error, not an absent match.
                parse_glycan_encoding(structure, classy)
            try:
                core.learn_fuc(
                    row["classy"],
                    row["name"],
                    structure,
                    float(observed_mz),
                    float(ms1_rt),
                    float(abundance),
                )
            except Exception as exc:
                logger.debug("Skipping candidate %s during legacy fragment generation: %s", row.get("name"), exc)
    finally:
        core.folder = old_folder
    if not out_db.exists() or out_db.stat().st_size == 0:
        return pd.DataFrame()
    fragments = pd.read_csv(out_db, header=None)
    fragments.columns = ["idf", "mz", "anno", "total", "target", "lenDB", "RT", "Abund"]
    fragments["mz"] = pd.to_numeric(fragments["mz"], errors="coerce")
    fragments = fragments.dropna(subset=["mz"])
    deduplicated = []
    for _, group in fragments.groupby("anno", sort=False):
        group = group.copy()
        labels = group.groupby("mz")["idf"].agg(lambda xs: " | ".join(sorted(set(xs.astype(str)))))
        group["idf"] = group["mz"].map(labels)
        deduplicated.append(group.drop_duplicates(subset="mz"))
    return pd.concat(deduplicated, ignore_index=True) if deduplicated else pd.DataFrame(columns=fragments.columns)


def match_observed_to_theoretical(
    spectrum: MS2Spectrum,
    theoretical: pd.DataFrame,
    fragment_tolerance_fraction: float,
    min_fragment_intensity: float,
) -> pd.DataFrame:
    """Match one observed MS2 spectrum to legacy theoretical fragments."""

    observed_mz = spectrum.fragment_mz
    observed_intensity = spectrum.fragment_intensity
    keep = np.isfinite(observed_mz) & np.isfinite(observed_intensity) & (observed_intensity >= min_fragment_intensity)
    observed_mz = observed_mz[keep]
    observed_intensity = observed_intensity[keep]
    if len(observed_mz) == 0 or theoretical.empty:
        return pd.DataFrame()

    theory = theoretical.copy()
    theory["low"] = theory["mz"] * (1 - fragment_tolerance_fraction)
    theory["up"] = theory["mz"] * (1 + fragment_tolerance_fraction)
    rows: list[pd.DataFrame] = []
    lows = theory["low"].to_numpy(dtype=float)
    ups = theory["up"].to_numpy(dtype=float)
    for observed_id, (mz_value, intensity_value) in enumerate(zip(observed_mz, observed_intensity)):
        mask = (lows <= mz_value) & (ups >= mz_value)
        if not np.any(mask):
            continue
        matched = theory.loc[mask].copy()
        matched["idx"] = spectrum.scan_id
        matched["_observed_id"] = observed_id
        matched["_observed_mz"] = float(mz_value)
        matched["_observed_intensity"] = float(intensity_value)
        matched["relmz"] = round(float(mz_value), 4)
        matched["intense"] = round(float(intensity_value), 4)
        matched["fragment_error_ppm"] = (matched["mz"].astype(float) - float(mz_value)).abs() / matched["mz"].astype(float).abs().clip(lower=1e-12) * 1_000_000
        rows.append(matched)
    if not rows:
        return pd.DataFrame()
    return select_one_to_one_edges(
        pd.concat(rows, ignore_index=True), observed_cols=["idx", "_observed_id"],
        observed_mz_col="_observed_mz", intensity_col="_observed_intensity",
        theoretical_mz_col="mz", candidate_cols=["anno"], label_col="idf",
    )


def finalize_spectrum_matches(
    matched: pd.DataFrame,
    min_matched_fragments: int,
    min_match_score: float,
    top_n: int,
) -> pd.DataFrame:
    """Apply the same coverage/intensity ranking used by the legacy workflow."""

    if matched.empty:
        return pd.DataFrame()
    matched = select_one_to_one_edges(
        matched, observed_cols=["idx", "_observed_id"] if "_observed_id" in matched else ["idx", "relmz", "intense"],
        observed_mz_col="_observed_mz" if "_observed_mz" in matched else "relmz",
        intensity_col="_observed_intensity" if "_observed_intensity" in matched else "intense",
        theoretical_mz_col="mz", candidate_cols=["anno"], label_col="idf",
    )
    cleaned_frames = [
        group.sort_values("intense", ascending=False).drop_duplicates("mz", keep="first")
        for _, group in matched.groupby("anno", sort=False)
    ]
    cleaned_data = pd.concat(cleaned_frames, ignore_index=True) if cleaned_frames else pd.DataFrame(columns=matched.columns)
    fragment_error_summary = pd.DataFrame()
    if "fragment_error_ppm" in cleaned_data.columns and not cleaned_data.empty:
        fragment_error_summary = cleaned_data.groupby("anno")["fragment_error_ppm"].agg(
            fragment_error_ppm_mean="mean",
            fragment_error_ppm_max="max",
        )
    middle_frames = [core.process_group(group) for _, group in cleaned_data.groupby("anno", sort=False)]
    middle = pd.concat(middle_frames, ignore_index=True) if middle_frames else pd.DataFrame()
    if middle.empty:
        return middle
    if not fragment_error_summary.empty:
        middle = middle.merge(fragment_error_summary, left_on="注释", right_index=True, how="left")
    middle = middle[middle["实际mz"].apply(lambda values: len(values) >= min_matched_fragments)]
    middle = middle[middle["匹配度分数"] >= min_match_score]
    if middle.empty:
        return middle
    middle["强度总和"] = middle["强度"].apply(sum)
    result = pd.DataFrame()
    for _, group in middle.groupby("target"):
        group = group.sort_values(by="强度总和", ascending=False)
        rank = group["强度总和"].rank(method="dense", ascending=False)
        group["总分数"] = (100 - (rank - 1) * 10) * group["匹配度分数"]
        max_intensity = group["强度总和"].max()
        group["相对强度"] = (group["强度总和"] / max_intensity).round(2)
        selected = group.sort_values(
            ["总分数", "匹配度分数", "强度总和", "注释"],
            ascending=[False, False, False, True], kind="stable",
        ).head(top_n)
        result = pd.concat([result, selected], ignore_index=True)
    return result.reset_index(drop=True)


def apply_rt_fit_filter(results: pd.DataFrame, config: TargetedMzMLConfig) -> pd.DataFrame:
    """Remove RANSAC-excluded RT outliers while keeping not-evaluated rows."""

    if results.empty:
        return results.copy()
    data = results.copy()
    annotation_col = "lipid_name" if "lipid_name" in data.columns else "注释"
    data = add_series_columns(data, annotation_col=annotation_col, drop_unparsed=False)
    data["rt_validation_status"] = "not_evaluated"
    data["predicted_rt"] = np.nan
    data["rt_error"] = np.nan

    modelable = data.dropna(subset=["series", "x", "ms1_rt"]).copy()
    if modelable.empty:
        return data
    counts = modelable.groupby("series")["identification_id"].transform("count")
    modelable = modelable[counts >= config.rt_fit_min_points].copy()
    if modelable.empty:
        return data

    inliers, models = fit_ransac_models(
        modelable,
        series_col="series",
        x_col="x",
        y_col="ms1_rt",
        residual_threshold=config.rt_fit_residual_threshold,
        min_points_small_series=config.rt_fit_min_points,
        min_points_large_series=config.rt_fit_min_points,
    )
    if models.empty:
        return data

    predicted = predict_rt(data, models, series_col="series", x_col="x")
    data["predicted_rt"] = predicted["pred_y"]
    data["rt_error"] = data["ms1_rt"] - data["predicted_rt"]
    modeled_ids = set(modelable.loc[modelable["series"].isin(models["series"]), "identification_id"])
    inlier_ids = set(inliers["identification_id"]) if "identification_id" in inliers.columns else set()
    data.loc[data["identification_id"].isin(modeled_ids), "rt_validation_status"] = "outlier"
    data.loc[data["identification_id"].isin(inlier_ids), "rt_validation_status"] = "pass"
    return data[data["rt_validation_status"] != "outlier"].reset_index(drop=True)


def _extract_scan_rt(spectrum: dict) -> float | None:
    scans = spectrum.get("scanList", {}).get("scan", [])
    if not scans:
        return None
    value = scans[0].get("scan start time")
    if value is None:
        return None
    unit = str(getattr(value, "unit_info", "minute")).lower()
    return float(value) / 60 if unit in {"second", "seconds", "s", "uo:0000010"} else float(value)


def _extract_precursor(spectrum: dict) -> dict[str, float | None] | None:
    precursors = spectrum.get("precursorList", {}).get("precursor", [])
    if not precursors:
        return None
    precursor = precursors[0]
    selected_list = precursor.get("selectedIonList", {}).get("selectedIon", [])
    selected_mz = None
    if selected_list:
        selected_mz = selected_list[0].get("selected ion m/z")
    isolation_mz = precursor.get("isolationWindow", {}).get("isolation window target m/z")
    precursor_mz = selected_mz if selected_mz is not None else isolation_mz
    if precursor_mz is None:
        return None
    activation = precursor.get("activation", {})
    return {
        "precursor_mz": float(precursor_mz),
        "collision_energy": float(activation["collision energy"]) if "collision energy" in activation else None,
    }


def _read_legacy_or_normalized_ms1_library(path: str | Path) -> pd.DataFrame:
    path = require_file(path, "MS1 theoretical library")

    if path.suffix.lower() == ".csv":
        return _read_masshunter_compound_csv(path)

    with pd.ExcelFile(path) as xls:
        frames = [xls.parse(sheet) for sheet in xls.sheet_names]
    if set(MS1_LIBRARY_REQUIRED_COLUMNS).issubset(frames[0].columns):
        return read_ms1_library(path)
    data = pd.concat(frames, ignore_index=True)
    mz_col = _first_existing(data, ["理论值", "calculated_m/z", "理论 m/z", "theoretical_mz"])
    name_col = _first_existing(data, ["name", "lipid_name", "Name"])
    class_col = _first_existing(data, ["classy", "main_class", "lipid_class"])
    if mz_col is None or name_col is None:
        raise ValueError(f"MS1 library cannot be normalized; missing m/z or lipid name columns: {path}")
    data = coerce_numeric_columns(data, (mz_col,), "MS1 theoretical library")
    normalized = pd.DataFrame(
        {
            "理论值": pd.to_numeric(data[mz_col], errors="coerce"),
            "classy": data[class_col] if class_col else data[name_col].astype(str).str.extract(r"^([^(]+)", expand=False),
            "name": data[name_col],
            "structure": data["structure"] if "structure" in data.columns else "",
        }
    )
    return normalized.dropna(subset=["理论值", "name"]).reset_index(drop=True)


def _read_masshunter_compound_csv(path: Path) -> pd.DataFrame:
    """Read Agilent MassHunter compound database CSV as a legacy MS1 library."""

    header_line = None
    with path.open("r", encoding="utf-8-sig", errors="ignore") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("# Formula"):
                header_line = line_number
                header_columns = [c.strip() for c in next(csv.reader([line.lstrip("# ")]))]
                break
    if header_line is None:
        raise ValueError(f"Cannot find '# Formula, RT, Mass, Cpd, Comments' header in {path}")

    last_error: Exception | None = None
    for encoding in ("utf-8-sig", "gb18030", "cp1252"):
        try:
            data = pd.read_csv(
                path,
                encoding=encoding,
                encoding_errors="ignore",
                skiprows=header_line+1,
                header=None,
                names=header_columns,
                dtype=str,
            )
            break
        except UnicodeError as exc:
            last_error = exc
    else:
        raise ValueError(f"Cannot decode MassHunter compound CSV: {path}") from last_error
    data = data.dropna(subset=["Formula", "Cpd"]).copy()
    mass = pd.to_numeric(data["Mass"], errors="coerce")
    computed_mass = data["Formula"].apply(_monoisotopic_mass_from_formula)
    neutral_mass = mass.fillna(computed_mass)
    normalized = pd.DataFrame(
        {
            "理论值": neutral_mass + PROTON_MASS,
            "classy": data["Cpd"].apply(_derive_legacy_classy),
            "name": data["Cpd"],
            "structure": data["Cpd"].apply(_derive_legacy_structure),
            "formula": data["Formula"],
            "neutral_mass": neutral_mass,
        }
    )
    structure_col = _first_existing(data,["structure","glycan_encoding"])
    position_col = _first_existing(data,["branch_positions","classy"])
    if structure_col:
        supplied = data[structure_col].fillna("").astype(str)
        mask = supplied.str.strip().ne("")
        normalized.loc[mask,"structure"] = supplied[mask]
        normalized.loc[mask,"classy"] = ""
        if position_col:
            normalized.loc[mask,"classy"] = data.loc[mask,position_col].fillna("").astype(str)
        for _, row in normalized.loc[mask].iterrows():
            parse_glycan_encoding(row["structure"],row["classy"])
    return normalized.dropna(subset=["理论值", "name"]).reset_index(drop=True)


ELEMENT_MASSES = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "Na": 22.9897692820,
}
PROTON_MASS = 1.007276466621


def _monoisotopic_mass_from_formula(formula: object) -> float:
    text = str(formula).replace(" ", "")
    if not text:
        return np.nan
    total = 0.0
    for element, count_text in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        if element not in ELEMENT_MASSES:
            return np.nan
        count = int(count_text) if count_text else 1
        total += ELEMENT_MASSES[element] * count
    return total


def _derive_legacy_classy(name: object) -> str:
    head = _lipid_head(name)
    if head in {"SM", "CerP", "S1P", "Lyso_SM", "Lyso-SM", "Lyso_sulfo", "CAEP", "N_CAEP"}:
        return head.replace("-", "_")
    if head in {"PE_Cer", "PECer"}:
        return "PE_cer"
    if head in {"PI_Cer", "PICer"}:
        return "PI_cer"
    if head in {"PG_Cer", "PGCer"}:
        return "PG_cer"
    if head in {"GlcSo", "Glu_So"}:
        return "Glu_So"
    if head in {"Gb3So", "Gb3_So"}:
        return "Gb3_So"
    if head in {"So", "Sa"}:
        return "So"
    if _derive_legacy_structure(name):
        return ""
    return head


def _derive_legacy_structure(name: object) -> str:
    head = _lipid_head(name)
    # Explicit, historically evidenced templates; do not infer branch topology
    # for GM1a/GM1b/GD1 isomers from a shared name prefix.
    if head == "GM1":
        return "Gal-GalNAc-Gal(-NeuAc)-Glc"
    if head.lower() == "type i b antigen":
        return "Gal-Gal(-Fuc)-GlcNAc-Gal-Glc"
    glycan_map = {
        "GlcCer": "-Glc",
        "GalCer": "-Gal",
        "LacCer": "-Gal -Glc",
        "Gb3Cer": "-Gal -Gal -Glc",
        "Gb4Cer": "-GalNAc -Gal -Gal -Glc",
        "GM3": "-NeuAc -Gal -Glc",
        "GM2": "-GalNAc -NeuAc -Gal -Glc",
        "GM1": "-Gal -GalNAc -NeuAc -Gal -Glc",
        "GD3": "-NeuAc -NeuAc -Gal -Glc",
        "GD2": "-GalNAc -NeuAc -NeuAc -Gal -Glc",
        "GD1": "-Gal -GalNAc -NeuAc -NeuAc -Gal -Glc",
    }
    for prefix, structure in glycan_map.items():
        if head.startswith(prefix):
            return structure
    return ""


def _lipid_head(name: object) -> str:
    text = str(name)
    return text.split("(", 1)[0]


def _first_existing(data: pd.DataFrame, names: Sequence[str]) -> str | None:
    for name in names:
        if name in data.columns:
            return name
    return None


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for targeted mzML annotation."""

    import argparse

    parser = argparse.ArgumentParser(description="Annotate targeted mzML MS2 spectra back to combined MS1 targetlist features.")
    parser.add_argument("--data-dir", required=True, type=Path, help="Root folder containing targetlist/ and mzML files.")
    parser.add_argument("--targetlist-dir", type=Path, default=None, help="Folder containing Agilent targetlist CSV files.")
    parser.add_argument("--mzml-dir", type=Path, default=None, help="Folder containing mzML files. Defaults to --data-dir.")
    parser.add_argument("--ms1-db", type=Path, help="MS1 Excel/CSV library; required only for --backend legacy.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output folder.")
    parser.add_argument("--fragment-ppm", type=float, default=TARGETED_FRAGMENT_PPM, help=f"MS2 fragment matching tolerance in ppm. Default: {TARGETED_FRAGMENT_PPM:g}.")
    parser.add_argument("--ms1-ppm", type=float, default=10.0, help="MS1 library candidate tolerance in ppm. Default: 10.")
    parser.add_argument("--feature-ppm", type=float, default=10.0, help="mzML precursor to targetlist feature tolerance in ppm. Default: 10.")
    parser.add_argument("--min-intensity", type=float, default=20.0)
    parser.add_argument("--min-fragments", type=int, default=2)
    parser.add_argument("--min-score", type=float, default=0.35)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--backend",choices=("lipidgate","legacy"),default="lipidgate")
    parser.add_argument("--ms2-library",type=Path)
    parser.add_argument("--fragment-da",type=float)
    parser.add_argument("--min-total-score",type=float,default=50.)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args(argv)

    logger = configure_logging(args.log_level.upper(), logger_name="sphingolipid_toolkit.targeted_mzml")
    config = TargetedMzMLConfig(
        data_dir=args.data_dir,
        targetlist_dir=args.targetlist_dir,
        mzml_dir=args.mzml_dir,
        ms1_library_path=args.ms1_db,
        output_dir=args.output_dir,
        fragment_ppm=args.fragment_ppm,
        ms1_candidate_ppm=args.ms1_ppm,
        feature_mz_ppm=args.feature_ppm,
        min_fragment_intensity=args.min_intensity,
        min_matched_fragments=args.min_fragments,
        min_match_score=args.min_score,
        top_n=args.top_n,
        backend=args.backend,
        ms2_library_path=args.ms2_library,
        fragment_da=args.fragment_da,
        min_total_score=args.min_total_score,
    )
    outputs = run_targeted_mzml_batch(config, logger=logger)
    print("Generated output files:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
