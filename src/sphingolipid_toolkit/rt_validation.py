"""Retention-time series extraction and RT-based validation utilities.

This module is the cleaned Python version of the RT notebooks originally named
``rentetion time.ipynb`` and ``RT_liner.ipynb``. The notebook cells were
converted into reusable functions and the hard-coded Windows paths were removed.

The core logic is intentionally close to the notebook workflow:
1. parse lipid annotations into homologous ``series`` and variable carbon number
   ``x``;
2. fit one RANSAC linear RT model per homologous series;
3. return inliers, equations, R2, predicted RT and deviations for downstream
   validation or plotting.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import r2_score

LYSO_OR_BASE_CLASSES = {
    "So",
    "So1P",
    "N-methylSo",
    "N,N-dimethylSo",
    "GlcSo",
    "Lyso-SM",
    "Lyso-sulfo",
    "Gb3So",
}


@dataclass(frozen=True)
class RTModelSummary:
    """Summary of one homologous-series RT model."""

    series: str
    slope: float
    intercept: float
    r2: float
    n_points: int
    n_inliers: int


def merge_excel_results(
    input_files: Sequence[str | Path],
    output_file: str | Path | None = None,
    sheet_name: str | int = 0,
    sort_by: str = "RT",
) -> pd.DataFrame:
    """Merge multiple result Excel files and optionally write the merged table.

    This replaces the notebook pattern that looped over ``D:/0-5MIN/result-i``
    and used the deprecated ``DataFrame.append`` method.
    """

    frames = [pd.read_excel(path, sheet_name=sheet_name) for path in input_files]
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if sort_by in merged.columns:
        merged = merged.sort_values(by=sort_by).reset_index(drop=True)
    if output_file is not None:
        merged.to_excel(output_file, index=False)
    return merged


def parse_lipid_series(annotation: str) -> tuple[str | None, int | None]:
    """Return ``(series, x)`` from a lipid annotation string.

    Examples
    --------
    ``Cer(d18:1/24:0)`` -> ``Cer(d18:1/x:0)``, ``24``

    For lyso/base sphingolipids, the varying carbon number is read from the
    first ``C:U`` pattern. For N-acyl sphingolipids, it is read after ``/``.
    """

    if not isinstance(annotation, str) or "(" not in annotation:
        return None, None

    head_match = re.search(r"^(.*?)(?=\()", annotation)
    if not head_match:
        return None, None
    head = head_match.group(1)

    if head in LYSO_OR_BASE_CLASSES:
        matches = re.findall(r"(\d+):", annotation)
        if not matches:
            return None, None
        carbon = int(matches[0])
        return "x".join(annotation.split(matches[0], 1)), carbon

    fa_match = re.search(r"/(\d+):", annotation)
    if not fa_match or "/" not in annotation:
        return None, None
    carbon = int(fa_match.group(1))
    left, right = annotation.split("/", 1)
    return left + "/" + "x".join(right.split(fa_match.group(1), 1)), carbon


def add_series_columns(
    data: pd.DataFrame,
    annotation_col: str = "注释",
    series_col: str = "series",
    carbon_col: str = "x",
    index_col: str = "index",
    drop_unparsed: bool = False,
) -> pd.DataFrame:
    """Add homologous-series labels and variable carbon number to a table."""

    out = data.copy()
    parsed = out[annotation_col].apply(parse_lipid_series)
    out[series_col] = parsed.apply(lambda item: item[0])
    out[carbon_col] = parsed.apply(lambda item: item[1])
    if drop_unparsed:
        out = out.dropna(subset=[series_col, carbon_col]).reset_index(drop=True)
    out[index_col] = pd.factorize(out[series_col], sort=False)[0] + 1
    out.loc[out[series_col].isna(), index_col] = np.nan
    return out


def fit_ransac_models(
    data: pd.DataFrame,
    series_col: str = "series",
    x_col: str = "x",
    y_col: str = "y",
    residual_threshold: float = 0.1,
    min_points_small_series: int = 3,
    min_points_large_series: int = 4,
    large_series_cutoff: int = 8,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit a RANSAC RT model for each homologous series.

    Parameters mirror the notebooks: residual threshold defaults to 0.1 min;
    if too few inliers are selected by RANSAC, the closest 3 or 4 points are
    forced to remain in the training set.

    Returns
    -------
    inlier_table:
        Rows kept as model inliers, with ``pred_y`` and ``rt_deviation``.
    model_table:
        One row per series with slope, intercept, R2 and point counts.
    """

    required = {series_col, x_col, y_col}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")

    inlier_frames: list[pd.DataFrame] = []
    summaries: list[RTModelSummary] = []

    clean = data.dropna(subset=[series_col, x_col, y_col]).copy()
    for series_name, group in clean.groupby(series_col, sort=False):
        group = group.sort_values(by=x_col).reset_index(drop=True)
        if len(group) < 2:
            continue

        x = group[x_col].astype(float).to_numpy().reshape(-1, 1)
        y = group[y_col].astype(float).to_numpy()
        model = RANSACRegressor(
            min_samples=2,
            residual_threshold=residual_threshold,
            random_state=random_state,
        )
        try:
            model.fit(x, y)
        except ValueError:
            continue
        y_pred = model.predict(x)
        inliers_mask = np.abs(y - y_pred) < residual_threshold

        min_required = min_points_small_series if len(group) <= large_series_cutoff else min_points_large_series
        min_required = min(min_required, len(group))
        if int(np.sum(inliers_mask)) < min_required:
            closest = np.argsort(np.abs(y - y_pred))[:min_required]
            inliers_mask[closest] = True

        filtered = group.loc[inliers_mask].copy()
        filtered["pred_y"] = model.predict(filtered[x_col].astype(float).to_numpy().reshape(-1, 1))
        filtered["rt_deviation"] = filtered[y_col].astype(float) - filtered["pred_y"]
        filtered["rt_abs_deviation"] = filtered["rt_deviation"].abs()
        inlier_frames.append(filtered)

        inlier_x = filtered[x_col].astype(float).to_numpy().reshape(-1, 1)
        inlier_y = filtered[y_col].astype(float).to_numpy()
        inlier_pred = model.predict(inlier_x)
        r2 = r2_score(inlier_y, inlier_pred) if len(filtered) >= 2 else np.nan
        estimator = model.estimator_
        summaries.append(
            RTModelSummary(
                series=str(series_name),
                slope=float(estimator.coef_[0]),
                intercept=float(estimator.intercept_),
                r2=float(r2),
                n_points=int(len(group)),
                n_inliers=int(len(filtered)),
            )
        )

    inlier_table = pd.concat(inlier_frames, ignore_index=True) if inlier_frames else pd.DataFrame()
    model_table = pd.DataFrame([summary.__dict__ for summary in summaries])
    return inlier_table, model_table


def predict_rt(data: pd.DataFrame, model_table: pd.DataFrame, series_col: str = "series", x_col: str = "x") -> pd.DataFrame:
    """Add predicted RT values to ``data`` using a model table from RANSAC."""

    out = data.copy()
    models = model_table.set_index("series")[["slope", "intercept"]].to_dict("index")

    def _predict(row: pd.Series) -> float:
        model = models.get(row.get(series_col))
        if model is None or pd.isna(row.get(x_col)):
            return np.nan
        return float(model["slope"] * row[x_col] + model["intercept"])

    out["pred_y"] = out.apply(_predict, axis=1)
    return out


def find_missing_series_members(
    data: pd.DataFrame,
    expected_range: Iterable[int] = range(14, 31),
    series_col: str = "series",
    x_col: str = "x",
) -> pd.DataFrame:
    """Report missing carbon numbers within each homologous series."""

    rows: list[dict[str, object]] = []
    expected = set(int(v) for v in expected_range)
    for series_name, group in data.dropna(subset=[series_col, x_col]).groupby(series_col):
        present = set(group[x_col].astype(int).tolist())
        for missing_x in sorted(expected - present):
            rows.append({series_col: series_name, "lost_x": missing_x, "annotation_template": str(series_name).replace("x", str(missing_x))})
    return pd.DataFrame(rows)


def match_predicted_to_targets(
    predicted: pd.DataFrame,
    targets: pd.DataFrame,
    rt_col_pred: str = "pred_y",
    rt_col_target: str = "RT",
    mz_cols_pred: Sequence[str] = ("MH_T", "MH+H", "MH+2H"),
    mz_col_target: str = "Precur",
    rt_tolerance: float = 0.1,
    mz_tolerance: float = 0.1,
) -> pd.DataFrame:
    """Match predicted RT/mass candidates to a target table.

    This is a vector-safe version of the notebook's binary-search matching cell.
    """

    if predicted.empty or targets.empty:
        return pd.DataFrame()

    rt_values = targets[rt_col_target].to_numpy(dtype=float)
    sorted_idx = np.argsort(rt_values)
    sorted_rt = rt_values[sorted_idx]
    target_mz = targets[mz_col_target].to_numpy(dtype=float)
    rows: list[pd.Series] = []

    for _, pred_row in predicted.iterrows():
        pred_rt = float(pred_row[rt_col_pred])
        lower = np.searchsorted(sorted_rt, pred_rt - rt_tolerance)
        upper = np.searchsorted(sorted_rt, pred_rt + rt_tolerance, side="right")
        candidate_masses = [float(pred_row[col]) for col in mz_cols_pred if col in pred_row and pd.notna(pred_row[col])]
        for target_pos in sorted_idx[lower:upper]:
            if not candidate_masses:
                continue
            if np.any(np.isclose(candidate_masses, target_mz[target_pos], atol=mz_tolerance)):
                rows.append(pd.concat([pred_row, targets.iloc[target_pos]], axis=0))

    return pd.DataFrame(rows)
