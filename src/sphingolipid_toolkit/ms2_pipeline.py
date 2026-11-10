"""Batch pipeline for the SphinGOlipID MS/MS matching module."""

from __future__ import annotations

import os
import re
import logging
from dataclasses import dataclass
from pathlib import Path
from time import time
from typing import List, Sequence, Tuple

import pandas as pd

from . import ms2_legacy_core as core
from .config import LEGACY_FRAGMENT_PPM, SphinGOlipIDConfig, normalize_file_indices, validate_batch_settings
from .io_utils import (
    read_ms2_spectra_text,
    validate_ms2_file_inputs,
    validate_output_paths,
    require_file,
    resolve_ms2_input_paths,
    write_excel_atomic,
)
from .logging_utils import get_logger
from .logging_utils import attach_file_handler
from .scoring import add_standard_score_columns


@dataclass
class MS2PipelineConfig:
    """Configuration for batch MS/MS matching.

    The defaults mirror the original notebook as closely as possible while
    making hard-coded paths configurable.
    """

    input_dir: Path
    output_dir: Path
    ms1_db_path: Path
    precursor_dir: Path | None = None
    file_indices: Sequence[int] = tuple(range(1, 7))
    text_pattern: str = "hilic-msms-{i}.txt"
    target_pattern: str = "hilic-msms-{i}.xlsx"
    result_pattern: str = "result_msms-{i}.xlsx"
    final_result_name: str = "ms2_annotation_results.xlsx"
    log_file_name: str = "run_log.txt"
    intermediate_dir_name: str = "intermediate"
    encoding: str = "GBK"
    fragment_ppm: float = LEGACY_FRAGMENT_PPM
    min_fragment_intensity: float = 20.0
    min_matched_fragments: int = 2
    min_match_score: float = 0.35
    top_n: int = 3
    keep_intermediate: bool = False
    backend: str = "legacy"
    ms2_library_path: Path | None = None
    precursor_ppm: float = 10.0
    fragment_da: float | None = None
    min_total_score: float = 50.0
    feature_rt_tolerance_min: float = 0.5

    def __post_init__(self) -> None:
        self.input_dir = Path(self.input_dir)
        self.output_dir = Path(self.output_dir)
        self.ms1_db_path = Path(self.ms1_db_path)
        if self.precursor_dir is not None:
            self.precursor_dir = Path(self.precursor_dir)
        self.file_indices = normalize_file_indices(self.file_indices)
        validate_batch_settings(self)

    @property
    def raw_ms2_dir(self) -> Path:
        return self.input_dir

    @property
    def ms1_library_path(self) -> Path:
        return self.ms1_db_path

    @property
    def save_intermediate(self) -> bool:
        return self.keep_intermediate

    @property
    def fragment_tolerance_fraction(self) -> float:
        return self.fragment_ppm / 1_000_000


ConfigLike = MS2PipelineConfig | SphinGOlipIDConfig


def _coerce_config(config: ConfigLike) -> MS2PipelineConfig:
    if isinstance(config, MS2PipelineConfig):
        return config
    from .engine import library_path
    return MS2PipelineConfig(
        backend=config.backend,
        ms2_library_path=config.ms2_library_path,
        precursor_ppm=config.precursor_ppm,
        fragment_da=config.fragment_da,
        min_total_score=config.min_total_score,
        feature_rt_tolerance_min=config.feature_rt_tolerance_min,
        input_dir=config.raw_ms2_dir,
        precursor_dir=config.precursor_dir,
        output_dir=config.output_dir,
        ms1_db_path=config.ms1_library_path or config.ms2_library_path or library_path('+'),
        file_indices=config.file_indices,
        text_pattern=config.text_pattern,
        target_pattern=config.target_pattern,
        result_pattern=config.result_pattern,
        final_result_name=config.final_result_name,
        log_file_name=config.log_file_name,
        intermediate_dir_name=config.intermediate_dir_name,
        encoding=config.encoding,
        fragment_ppm=config.fragment_ppm,
        min_fragment_intensity=config.min_fragment_intensity,
        min_matched_fragments=config.min_matched_fragments,
        min_match_score=config.min_match_score,
        top_n=config.top_n,
        keep_intermediate=config.save_intermediate,
    )


def _ensure_clean_file(path: Path) -> None:
    if path.exists():
        path.unlink()


def _read_spectrum_text(text_file: Path, encoding: str) -> List[str]:
    return read_ms2_spectra_text(text_file, encoding=encoding)


def _write_empty_result(result_xlsx: Path) -> pd.DataFrame:
    columns = [
        "注释", "target", "匹配度分数", "母离子", "RT", "Abund",
        "实际mz", "强度", "碎片", "强度总和", "总分数", "相对强度",
    ]
    result = pd.DataFrame(columns=columns)
    result = add_standard_score_columns(result)
    write_excel_atomic(result, result_xlsx, index=False)
    return result


def _write_db_sheets(out_db_csv: Path, db_data_xlsx: Path, parameter: float) -> bool:
    """Convert generated theoretical fragments to target-grouped sheets.

    Returns ``False`` when no theoretical fragments were generated. This avoids
    a hard crash for files with no MS1-library candidates.
    """

    if not out_db_csv.exists() or out_db_csv.stat().st_size == 0:
        return False

    # Original notebook writes out_DB.csv without a header. Reading with
    # header=None avoids silently dropping the first theoretical fragment row.
    dfs = pd.read_csv(out_db_csv, header=None)
    if dfs.empty:
        return False
    dfs.columns = ["idf", "mz", "anno", "total", "target", "lenDB", "RT", "Abund"]
    dfs["mz"] = pd.to_numeric(dfs["mz"], errors="coerce").round(4)
    dfs = dfs.dropna(subset=["mz"])
    dfs = dfs.assign(
        下限=lambda x: x["mz"] * (1 - parameter),
        上限=lambda x: x["mz"] * (1 + parameter),
    ).fillna("")

    # Keep the original scoring idea: duplicated theoretical fragment m/z values
    # within the same annotation are counted once.
    def merge_labels(group: pd.DataFrame) -> pd.DataFrame:
        group = group.copy()
        labels = group.groupby("mz")["idf"].agg(lambda xs: " | ".join(sorted(set(xs.astype(str)))))
        group["idf"] = group["mz"].map(labels)
        return group.drop_duplicates(subset="mz")
    resultp = pd.concat(
        [merge_labels(group) for _, group in dfs.groupby(["target", "RT", "Abund", "anno"], sort=False)],
        ignore_index=True,
    ) if not dfs.empty else dfs
    if resultp.empty:
        return False
    grouped_s = resultp.groupby("target")

    with pd.ExcelWriter(db_data_xlsx, engine="openpyxl") as writer:
        for target, group in grouped_s:
            group.to_excel(writer, sheet_name=str(float(target)), index=False)
    return True


def _finalize_results(out_query_csv: Path, result_xlsx: Path, min_matched_fragments: int, min_match_score: float, top_n: int) -> pd.DataFrame:
    if not out_query_csv.exists() or out_query_csv.stat().st_size == 0:
        return _write_empty_result(result_xlsx)

    # Original notebook writes out-query.csv without a header. header=None keeps
    # the first matched fragment row.
    df_max = pd.read_csv(out_query_csv, header=None)
    if df_max.empty:
        return _write_empty_result(result_xlsx)
    df_max.columns = [
        "idf", "mz", "anno", "total", "target", "lenDB", "RT", "Abund",
        "low", "up", "idx", "relmz", "intense",
    ]
    df_max["relmz"] = df_max["relmz"].round(4)
    df_max["intense"] = df_max["intense"].round(4)

    # For duplicated theoretical fragments in one annotation, keep the measured
    # fragment with the highest intensity. This follows the notebook logic.
    candidate_keys = ["target", "RT", "Abund", "anno"]
    cleaned_data = df_max.sort_values("intense", ascending=False).drop_duplicates(candidate_keys + ["mz"])
    middle_df = pd.concat(
        [core.process_group(group) for _, group in cleaned_data.groupby(candidate_keys, sort=False)],
        ignore_index=True,
    )
    if middle_df.empty:
        return _write_empty_result(result_xlsx)
    middle_df = middle_df[middle_df["实际mz"].apply(lambda x: len(x) >= min_matched_fragments)]
    middle_df = middle_df[middle_df["匹配度分数"] >= min_match_score]
    if middle_df.empty:
        return _write_empty_result(result_xlsx)
    middle_df["强度总和"] = middle_df["强度"].apply(lambda x: sum(x))

    result = pd.DataFrame()
    for _, group in middle_df.groupby(["target", "RT", "Abund"], sort=False):
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

    result.reset_index(drop=True, inplace=True)
    result = add_standard_score_columns(result)
    write_excel_atomic(result, result_xlsx, index=False)
    return result


@core.synchronized_legacy_runtime
def run_one_file(config: ConfigLike, file_index: int, logger: logging.Logger | None = None) -> Path:
    """Run the MS/MS matching workflow for one HILIC-MS/MS export file."""

    config = _coerce_config(config)
    logger = logger or get_logger("sphingolipid_toolkit.ms2")
    if config.backend == "lipidgate":
        from .reference_pipeline import run_text_file
        return run_text_file(config, file_index, logger)
    logger.info("Effective fragment_ppm=%s (legacy text workflow config)", config.fragment_ppm)
    input_dir = Path(config.input_dir)
    precursor_dir = Path(config.precursor_dir or config.input_dir)
    output_dir = Path(config.output_dir)
    intermediate_dir = output_dir / config.intermediate_dir_name
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    text_file = input_dir / config.text_pattern.format(i=file_index)
    target_file = precursor_dir / config.target_pattern.format(i=file_index)
    ms2_data = intermediate_dir / f"actual_ms2_data-{file_index}.xlsx"
    db_data = intermediate_dir / f"DB_sheet-{file_index}.xlsx"
    result_file = output_dir / config.result_pattern.format(i=file_index)
    isomer_file = intermediate_dir / "Isomer.xlsx"
    out_db_csv = intermediate_dir / "out_DB.csv"
    out_query_csv = intermediate_dir / "out-query.csv"
    parent_ion_file = intermediate_dir / "目标母离子.xlsx"

    try:
        logger.info("Loading MS1 library")
        logger.info("Validating inputs for file index %s", file_index)
        _, target_df = validate_ms2_file_inputs(config, file_index)
    except Exception as exc:
        logger.error("%s", exc)
        raise

    for p in [isomer_file, out_db_csv, out_query_csv, ms2_data, db_data, parent_ion_file]:
        _ensure_clean_file(p)

    # Set legacy-module runtime globals. This keeps the original function bodies
    # largely unchanged but removes absolute D:/ paths from normal use.
    core.folder = str(intermediate_dir) + os.sep
    core.ms2_data = str(ms2_data)
    core.DB_data = str(db_data)
    core.MS1_DB_PATH = str(config.ms1_db_path)
    core.MIN_FRAGMENT_INTENSITY = config.min_fragment_intensity
    core.OH_dict = {'d': 2, '5': 2, 'ω': 2, 'm': 1, '6': 3, 't': 3, 'q': 4}

    logger.info("Reading MS2 file: %s", text_file)
    data_list = _read_spectrum_text(text_file, config.encoding)
    core.target_df = target_df

    logger.info("Extracting precursor spectra from %s", target_file)
    core.real_mz_func(data_list)
    logger.info("Searching MS1 candidates from %s", config.ms1_db_path)
    core.isomer_search(target_df)

    df_isomer = pd.read_excel(isomer_file).fillna("")
    logger.info("Matching precursor candidates for file index %s", file_index)
    if df_isomer.empty:
        logger.warning("No candidate found for file index %s", file_index)
    logger.info("Generating theoretical fragments for %s candidate rows", len(df_isomer))
    for _, row in df_isomer.iterrows():
        core.frag_id1 = ["[C2H5NO+H]+"]
        core.frag_mz1 = [60.044]
        core.learn_fuc(
            row["classy"],
            row["name"],
            row["structure"],
            row["target"],
            row["RT"],
            row["Abund"],
        )

    has_theoretical_fragments = _write_db_sheets(out_db_csv, db_data, config.fragment_tolerance_fraction)
    if not has_theoretical_fragments:
        logger.warning("No theoretical fragments generated for file index %s", file_index)
        _write_empty_result(result_file)
        if not config.keep_intermediate:
            for p in [isomer_file, out_db_csv, out_query_csv, ms2_data, db_data, parent_ion_file]:
                if p.exists():
                    p.unlink()
        return result_file

    with pd.ExcelFile(ms2_data) as df_rel, pd.ExcelFile(db_data, engine="openpyxl") as attr_df:
        start_time = time()
        logger.info("Matching precursor fragment sheets for file index %s", file_index)
        for sheet_name in attr_df.sheet_names:
            if sheet_name not in df_rel.sheet_names:
                logger.warning("No measured spectrum for precursor %s", sheet_name)
                continue
            core.process_sheet(sheet_name, df_rel, attr_df)
        logger.info("MS2 matching for file %s finished in %.2f s", file_index, time() - start_time)

    logger.info("Exporting results: %s", result_file)
    _finalize_results(
        out_query_csv,
        result_file,
        min_matched_fragments=config.min_matched_fragments,
        min_match_score=config.min_match_score,
        top_n=config.top_n,
    )

    if not config.keep_intermediate:
        for p in [isomer_file, out_db_csv, out_query_csv, ms2_data, db_data, parent_ion_file]:
            if p.exists():
                p.unlink()

    logger.info("MS2 matching completed for file index %s", file_index)
    return result_file


@core.synchronized_legacy_runtime
def run_batch(config: ConfigLike, logger: logging.Logger | None = None) -> List[Path]:
    """Run matching for every index in ``config.file_indices``."""

    config = _coerce_config(config)
    validate_batch_settings(config)
    validate_output_paths(config, config.file_indices)
    require_file(config.ms1_db_path, "MS1 theoretical library")
    for i in config.file_indices:
        paths = resolve_ms2_input_paths(config, i)
        require_file(paths.text_file, f"MS/MS txt file for index {i}")
        require_file(paths.target_file, f"precursor target Excel for index {i}")
    logger = logger or get_logger("sphingolipid_toolkit.ms2")
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / config.intermediate_dir_name).mkdir(parents=True, exist_ok=True)
    file_handler = attach_file_handler(output_dir / config.log_file_name, logger=logger)
    results = []
    try:
        logger.info("Starting MS2 batch")
        for i in config.file_indices:
            results.append(run_one_file(config, i, logger=logger))
        _write_batch_result(results, output_dir / config.final_result_name, logger)
        logger.info("MS2 batch completed")
    finally:
        logger.removeHandler(file_handler)
        file_handler.close()
    return results


def _write_batch_result(result_files: Sequence[Path], final_result_file: Path, logger: logging.Logger) -> pd.DataFrame:
    """Write the v0.3 final aggregate result file."""

    frames: list[pd.DataFrame] = []
    for result_file in result_files:
        if not Path(result_file).exists():
            continue
        frame = pd.read_excel(result_file)
        frame["source_result_file"] = Path(result_file).name
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True) if frames else add_standard_score_columns(pd.DataFrame())
    result = add_standard_score_columns(result)
    logger.info("Exporting results: %s", final_result_file)
    write_excel_atomic(result, final_result_file, index=False)
    return result


def parse_file_indices(value: str) -> Tuple[int, ...]:
    """Parse values like ``1-6`` or ``1,3,5`` for the CLI."""

    value = value.strip()
    if re.fullmatch(r"\d+-\d+", value):
        a, b = map(int, value.split("-"))
        if b < a:
            raise ValueError("file index range end must be >= start")
        return tuple(range(a, b + 1))
    indices = tuple(int(x.strip()) for x in value.split(",") if x.strip())
    if not indices or any(i < 0 for i in indices):
        raise ValueError("file indices must contain at least one non-negative integer")
    return tuple(dict.fromkeys(indices))
