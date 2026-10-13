"""Duplicate-removal and result-cleaning utilities.

This module is the cleaned Python version of ``去除重复值.ipynb``. The notebook
contained one incomplete cell and several top-level file paths; those have been
converted into reusable functions that can be called from CLI/GUI code.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from .rt_validation import add_series_columns


def reset_series_index(
    data: pd.DataFrame,
    series_col: str = "series",
    min_count: int = 5,
    index_col: str = "index",
) -> pd.DataFrame:
    """Keep series with at least ``min_count`` rows and assign compact indices."""

    if series_col not in data.columns:
        raise KeyError(f"Missing required column: {series_col}")
    out = data.copy()
    counts = out[series_col].value_counts(dropna=True)
    valid = counts[counts >= min_count].index
    out = out[out[series_col].isin(valid)].sort_values(by=series_col).reset_index(drop=True)
    out[index_col] = out.groupby(series_col, sort=False).ngroup() + 1
    return out


def deduplicate_by_rt_score_intensity(
    data: pd.DataFrame,
    group_cols: Sequence[str] = ("series", "x"),
    rt_col: str = "y",
    score_col: str = "匹配度分数",
    intensity_col: str = "强度总和",
    rt_tolerance: float = 0.1,
) -> pd.DataFrame:
    """Remove near-duplicate annotations within a homologous-series group.

    For rows within ``rt_tolerance`` minutes, keep the row with the higher
    matching score; if scores tie, keep the higher total fragment intensity.
    Select greedily in descending score/intensity order, then ascending RT.
    Compare every candidate with *all retained* RTs in its group. Exact ties
    use serialized row content (never the input index) for reproducibility.
    Missing RTs cannot establish proximity and are retained independently.
    """

    required = set(group_cols) | {rt_col, score_col, intensity_col}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")

    if rt_tolerance < 0:
        raise ValueError("rt_tolerance must be non-negative")
    work = data.reset_index(drop=True)
    kept_groups: list[pd.DataFrame] = []

    for _, group in work.groupby(list(group_cols), group_keys=False, sort=True):
        # A separate ordering table avoids changing input columns or dtypes.
        ordering = pd.DataFrame({
            "score": pd.to_numeric(group[score_col]),
            "intensity": pd.to_numeric(group[intensity_col]),
            "rt": pd.to_numeric(group[rt_col]),
            "tie": group.apply(
                lambda row: repr(tuple(
                    value.tolist() if isinstance(value, np.ndarray) else value
                    for value in row
                )), axis=1,
            ),
        }).sort_values(
            ["score", "intensity", "rt", "tie"],
            ascending=[False, False, True, True], kind="stable", na_position="last",
        )
        keep_indices: list[int] = []
        kept_rts: list[float] = []
        for i, candidate in ordering.iterrows():
            rt = candidate["rt"]
            if pd.isna(rt) or not any(abs(rt - kept) <= rt_tolerance for kept in kept_rts):
                keep_indices.append(i)
                if pd.notna(rt):
                    kept_rts.append(float(rt))
        kept_groups.append(group.loc[keep_indices])

    return pd.concat(kept_groups, ignore_index=True) if kept_groups else pd.DataFrame(columns=data.columns)


def filter_target_abundance_duplicates(
    data: pd.DataFrame,
    group_cols: Sequence[str] = ("target", "Abund"),
    index_col: str = "index",
    total_score_col: str = "总分数",
) -> pd.DataFrame:
    """For each target/abundance group, keep the most consistent series index.

    The notebook selected the most frequent ``index`` within each
    ``target``/``Abund`` group; ties were resolved by the highest ``总分数``.
    """

    required = set(group_cols) | {index_col, total_score_col}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")

    def _filter(group: pd.DataFrame) -> pd.Series:
        counts = group[index_col].value_counts()
        most_common = counts.index[0]
        candidates = group[group[index_col] == most_common]
        if len(candidates) > 1:
            return candidates.loc[candidates[total_score_col].idxmax()]
        return candidates.iloc[0]

    return data.groupby(list(group_cols), group_keys=False).apply(_filter).reset_index(drop=True)


def remove_duplicate_annotations(
    data: pd.DataFrame,
    annotation_col: str = "注释",
    add_series: bool = True,
    rt_col: str = "y",
    min_series_count: int = 5,
    rt_tolerance: float = 0.1,
) -> pd.DataFrame:
    """Convenience wrapper for the common cleaning workflow.

    It can parse ``series``/``x`` from annotation names, remove close RT
    duplicates by score/intensity, then reset homologous-series indices.
    """

    out = data.copy()
    if add_series and ("series" not in out.columns or "x" not in out.columns):
        out = add_series_columns(out, annotation_col=annotation_col, drop_unparsed=True)
    out = deduplicate_by_rt_score_intensity(out, rt_col=rt_col, rt_tolerance=rt_tolerance)
    out = reset_series_index(out, min_count=min_series_count)
    return out


def clean_excel_file(
    input_file: str | Path,
    output_file: str | Path,
    sheet_name: str | int = 0,
    **kwargs,
) -> pd.DataFrame:
    """Read an Excel result table, clean duplicate annotations, and save it."""

    data = pd.read_excel(input_file, sheet_name=sheet_name)
    cleaned = remove_duplicate_annotations(data, **kwargs)
    cleaned.to_excel(output_file, index=False)
    return cleaned
