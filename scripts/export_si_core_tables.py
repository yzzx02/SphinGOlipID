"""Export compact supplementary-information tables for the six-fraction study.

The source RT result directories contain wide audit tables intended for internal
quality control.  This script keeps only the columns needed to reproduce and
interpret the reported identification, ECN/IUP filtering, and fitted-line
results.  It writes four flat UTF-8 CSV files; workbook packaging is deliberately
kept separate from the scientific data transformation.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from regenerate_top3_rt_ecn_iup_preview import PLOT_STYLE, prepare_groupings  # noqa: E402
from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_ecn, fit_rt_iup  # noqa: E402


MODE_FOLDERS = {
    "detailed_chain": "详细链分类",
    "lcb_series": "长链碱系列分类",
    "total_carbon_unsaturation": "总碳数总不饱和度分类",
}

RAW_COLUMNS = {
    "合并后行ID": "result_id",
    "时间窗口": "hilic_fraction_min",
    "细类": "lipid_subclass",
    "注释": "lipid_annotation",
    "母离子": "precursor_annotation",
    "候选排名": "candidate_rank",
    "匹配度分数": "match_score",
    "总分数": "total_score",
    "丰度": "ms2_abundance",
    "x碳数": "variable_carbon_number",
    "曲线不饱和度": "curve_unsaturation",
    "归一化保留时间": "normalized_rt_min",
    "二级谱图RT均值": "mean_ms2_rt_min",
    "合并谱图数": "merged_ms2_spectrum_count",
    "特征ID": "ms1_feature_id",
    "MS2扫描ID": "ms2_scan_id",
}

FILTER_COLUMNS = {
    "过滤口径": "filter_definition",
    "初始总候选": "initial_candidate_pairs",
    "过滤后剩余": "retained_candidate_pairs",
    "MS1 feature": "unique_ms1_features",
    "MS1 RT捞回": "ms1_rt_rescued_pairs",
}

FILTERED_COLUMNS = {
    "合并后行ID": "result_id",
    "时间窗口": "hilic_fraction_min",
    "细类": "lipid_subclass",
    "注释": "lipid_annotation",
    "特征ID": "ms1_feature_id",
    "候选排名": "candidate_rank",
    "匹配度分数": "match_score",
    "归一化保留时间": "normalized_rt_min",
    "x碳数": "carbon_number",
    "曲线不饱和度": "unsaturation",
    "拟合类型": "fit_type",
    "R²": "fit_r_squared",
    "保留时间绝对偏差": "absolute_rt_residual_min",
    "IUP_RT漂移校正(min)": "iup_vertical_shift_min",
    "判定状态": "fit_status",
    "保留原因": "retention_reason",
}

LINE_COLUMNS = {
    "细类": "lipid_subclass",
    "不饱和度": "unsaturation",
    "拟合类型": "fit_type",
    "参数": "fit_parameters",
    "R²": "fit_r_squared",
    "点数": "raw_point_count",
    "不同碳数": "distinct_carbon_count",
    "代表点数": "representative_point_count",
    "内点数": "inlier_count",
    "x最小": "carbon_min",
    "x最大": "carbon_max",
    "拟合方程": "fit_equation",
    "IUP_RT漂移校正(min)": "iup_vertical_shift_min",
    "平行性是否通过": "parallelism_passed",
    "平行性参考曲线": "parallelism_reference_unsaturation",
    "斜率绝对差": "median_absolute_slope_difference",
    "斜率相对差": "median_relative_slope_difference",
    "平行性说明": "parallelism_note",
}


def read_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    frame.columns = [str(column).lstrip("\ufeff") for column in frame.columns]
    return frame


def compact_columns(frame: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Select mapped columns, adding blank optional fields when absent."""

    output = pd.DataFrame(index=frame.index)
    for source, destination in mapping.items():
        output[destination] = frame[source] if source in frame.columns else pd.NA
    return output


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def build_raw_identifications(input_path: Path) -> pd.DataFrame:
    raw = read_csv(input_path)
    missing = [column for column in RAW_COLUMNS if column not in raw.columns]
    if missing:
        raise ValueError(f"Input table is missing required columns: {missing}")
    if raw["合并后行ID"].nunique() != len(raw):
        raise ValueError("The six-fraction input must have one row per 合并后行ID")
    output = compact_columns(raw, RAW_COLUMNS)
    return output.sort_values(
        ["hilic_fraction_min", "ms1_feature_id", "candidate_rank", "lipid_annotation"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_filter_summary(results_root: Path) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for mode, folder in MODE_FOLDERS.items():
        frame = read_csv(results_root / folder / "过滤口径统计.csv")
        selected = compact_columns(frame, FILTER_COLUMNS)
        selected.insert(0, "classification_mode", mode)
        parts.append(selected)
    return pd.concat(parts, ignore_index=True)


def build_filtered_results(
    input_path: Path,
    results_root: Path,
    score_threshold: float,
    include_low_score: bool,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for mode, folder in MODE_FOLDERS.items():
        for rule in ("ECN", "IUP"):
            frame = read_csv(results_root / folder / f"{rule}_表" / "最终保留明细.csv")
            score = pd.to_numeric(frame["匹配度分数"], errors="coerce")
            if score.lt(score_threshold).any():
                raise ValueError(f"{mode}/{rule} contains rows below score {score_threshold}")
            selected = compact_columns(frame, FILTERED_COLUMNS)
            if rule == "ECN":
                selected["iup_vertical_shift_min"] = 0.0
            selected.insert(0, "analysis_scope", "formal_score_ge_0.50_fit")
            selected.insert(0, "score_tier", "high_score_ge_0.50")
            selected.insert(0, "rt_rule", rule)
            selected.insert(0, "classification_mode", mode)
            parts.append(selected)

    if include_low_score:
        raw = read_csv(input_path)
        groupings = prepare_groupings(raw, -np.inf)
        config = replace(RTIUPConfig(), **PLOT_STYLE)
        folder_to_mode = {folder: mode for mode, folder in MODE_FOLDERS.items()}
        for folder, universe in groupings.items():
            mode = folder_to_mode[folder]
            for rule, result in (
                ("ECN", fit_ecn(universe, config)),
                ("IUP", fit_rt_iup(universe, config)),
            ):
                score = pd.to_numeric(result.rows["匹配度分数"], errors="coerce")
                low_rows = result.rows.loc[score.lt(score_threshold)].copy()
                selected = compact_columns(low_rows, FILTERED_COLUMNS)
                if rule == "ECN":
                    selected["iup_vertical_shift_min"] = 0.0
                selected.insert(0, "analysis_scope", "all_score_fit_low_score_rows")
                selected.insert(0, "score_tier", "low_score_lt_0.50")
                selected.insert(0, "rt_rule", rule)
                selected.insert(0, "classification_mode", mode)
                parts.append(selected)
    return pd.concat(parts, ignore_index=True)


def build_fit_lines(results_root: Path) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for mode, folder in MODE_FOLDERS.items():
        for rule in ("ECN", "IUP"):
            frame = read_csv(results_root / folder / f"{rule}_表" / "有效拟合曲线.csv")
            selected = compact_columns(frame, LINE_COLUMNS)
            if rule == "ECN":
                selected["iup_vertical_shift_min"] = 0.0
            selected.insert(0, "rt_rule", rule)
            selected.insert(0, "classification_mode", mode)
            parts.append(selected)
    return pd.concat(parts, ignore_index=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Canonical six-fraction TOP3 RT input CSV.")
    parser.add_argument("--results-root", required=True, type=Path, help="Final ECN/IUP result directory.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--score-threshold", type=float, default=0.50)
    parser.add_argument(
        "--include-low-score",
        action="store_true",
        help="Refit the full input and append RT-retained rows below the score threshold.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "SI01_raw_identifications.csv": build_raw_identifications(args.input),
        "SI02_filter_summary.csv": build_filter_summary(args.results_root),
        "SI03_rt_filtered_all_scores.csv": build_filtered_results(
            args.input,
            args.results_root,
            args.score_threshold,
            args.include_low_score,
        ),
        "SI04_fit_lines.csv": build_fit_lines(args.results_root),
    }
    for filename, frame in tables.items():
        write_csv(frame, args.output_dir / filename)
        print(f"{filename}: {len(frame):,} rows x {len(frame.columns)} columns")


if __name__ == "__main__":
    main()
