"""Retention-time ECN/IUP validation and plotting utilities.

This module contains the retention-time workflow used by the local analysis
scripts: per-unsaturation RANSAC linear/quadratic fitting, IUP ordering checks,
strict rescue of missing curves from raw candidates, and publication-style RT
plots with 95% confidence bands.
"""

from __future__ import annotations

import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import optimize, stats
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
import warnings


warnings.filterwarnings("ignore", category=UndefinedMetricWarning)


DEFAULT_COLORS = [
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#8E44AD",
    "#7E6148",
    "#B09C85",
]


@dataclass(frozen=True)
class RTIUPConfig:
    """Configuration for RT/IUP fitting and plotting."""

    rt_window_min: float = 0.2
    r2_threshold: float = 0.99
    iup_tolerance_min: float = 0.05
    order_tolerance_min: float = 0.08
    flat_y_range_min: float = 0.12
    flat_slope_min: float = 0.03
    rescue_min_distinct_x: int = 3
    short_series_max_crossed_unsats: int = 1
    max_refit_iterations: int = 5
    quadratic_min_distinct_x: int = 4
    quadratic_min_r2_gain: float = 0.002
    parallel_min_overlap_x: float = 2.0
    parallel_abs_slope_tolerance: float = 0.10
    parallel_relative_slope_tolerance: float = 0.25
    enable_iup_rescue: bool = True
    iup_multistart_trials: int = 16
    iup_max_candidates_per_unsaturation: int = 12
    iup_beam_width: int = 192
    random_state: int = 42
    colors: tuple[str, ...] = tuple(DEFAULT_COLORS)
    figure_size: tuple[float, float] = (7.2, 7.2)
    point_size: float = 54.0
    two_point_line_width: float = 1.9
    fit_line_width: float = 2.2
    legend_font_size: int = 14
    axis_label_size: int = 18
    tick_label_size: int = 16


@dataclass(frozen=True)
class RTIUPResult:
    """Tables returned by :func:`fit_rt_iup`."""

    rows: pd.DataFrame
    lines: pd.DataFrame
    plot_rows: pd.DataFrame
    stats: dict[str, object]


def finite(value: object) -> float | None:
    """Return a finite float or ``None``."""

    try:
        parsed = float(value)
    except Exception:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def relative_rt_tolerance_min(config: RTIUPConfig = RTIUPConfig()) -> float:
    """Return the worst-case relative RT allowance for two independent features.

    Each feature keeps the configured ``rt_window_min`` drift allowance. When
    two features or fitted series are compared for IUP order, their shifts can
    occur in opposite directions, and the configured IUP tolerance is added on
    top. The default is therefore ``0.20 + 0.20 + 0.05 = 0.45 min`` without
    changing either locked threshold.
    """

    return 2.0 * float(config.rt_window_min) + float(config.iup_tolerance_min)


def bool_value(value: object) -> bool:
    """Normalize bool-like values used in CSV/Excel payloads."""

    if isinstance(value, bool):
        return value
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except Exception:
        pass
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "是"}


def safe_name(value: object) -> str:
    """Make a Windows-safe plot filename stem."""

    text = re.sub(r'[\\/:*?"<>|]', "_", str(value))
    text = re.sub(r"\s+", "_", text).strip("._ ")
    return text[:150] or "plot"


def normalize_json_value(value: object) -> object:
    """Convert numpy/pandas scalars into JSON-safe values."""

    if value is None:
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, float):
        return None if not math.isfinite(value) else value
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    return value


def parse_series_unsaturation(series: object) -> float | None:
    """Parse the unsaturation value from a homologous series like ``x:2``."""

    if pd.isna(series):
        return None
    match = re.search(r"x:(\d+(?:\.\d+)?)", str(series))
    if not match:
        return None
    return finite(match.group(1))


def series_to_plot_group(series: object) -> str:
    """Replace ``x:<unsaturation>`` with ``x:y`` for grouped plots."""

    if pd.isna(series):
        return ""
    return re.sub(r"x:\d+(?:\.\d+)?", "x:y", str(series))


def color_for_unsaturation(value: object, config: RTIUPConfig = RTIUPConfig()) -> str:
    parsed = finite(value)
    if parsed is None:
        return "#7F7F7F"
    return config.colors[int(parsed) % len(config.colors)]


def _model_params(model: object, fit_type: str) -> list[float]:
    if fit_type == "Linear":
        estimator = model.estimator_
        return [float(estimator.coef_[0]), float(estimator.intercept_)]
    estimator = model.estimator_
    lin = estimator.named_steps["linearregression"]
    return [float(lin.coef_[1]), float(lin.coef_[0]), float(lin.intercept_)]


def predict_values(fit_type: object, params: object, x: np.ndarray | float) -> np.ndarray | float | None:
    """Predict RT values from serialized linear/quadratic parameters."""

    if isinstance(params, str):
        try:
            params = json.loads(params)
        except Exception:
            return None
    if not isinstance(params, (list, tuple)):
        return None
    arr = np.asarray(x, dtype=float)
    if fit_type == "Linear" and len(params) == 2:
        result = float(params[0]) * arr + float(params[1])
    elif fit_type == "Quadratic" and len(params) == 3:
        result = float(params[0]) * arr * arr + float(params[1]) * arr + float(params[2])
    else:
        return None
    if np.isscalar(x):
        return float(np.asarray(result))
    return result


def fit_prediction(line: pd.Series, x: np.ndarray | float) -> np.ndarray | float | None:
    return predict_values(line.get("拟合类型"), line.get("参数"), x)


def _is_positive_monotonic(fit_type: str, params: list[float], x_min: float, x_max: float, config: RTIUPConfig) -> bool:
    if x_max <= x_min:
        return False
    x_grid = np.linspace(x_min, x_max, 100)
    y_grid = predict_values(fit_type, params, x_grid)
    if y_grid is None:
        return False
    y_grid = np.asarray(y_grid, dtype=float)
    if not np.all(np.isfinite(y_grid)):
        return False
    return bool((y_grid[-1] - y_grid[0]) > 0 and np.mean(np.diff(y_grid) < -0.01) <= 0.05)


def _median_points_by_x(points: pd.DataFrame) -> pd.DataFrame:
    clean = points[["x碳数", "归一化保留时间"]].dropna().copy()
    if clean.empty:
        return clean
    return (
        clean.astype({"x碳数": float, "归一化保留时间": float})
        .groupby("x碳数", as_index=False)["归一化保留时间"]
        .median()
        .sort_values("x碳数")
        .reset_index(drop=True)
    )


def _has_clear_decreasing_tail(points: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> bool:
    """Detect obvious non-IUP tail drops in sorted scatter points."""

    ordered = _median_points_by_x(points)
    if len(ordered) < 3:
        return False
    y = ordered["归一化保留时间"].to_numpy(dtype=float)
    diffs = np.diff(y)
    clear_drop = diffs < -config.order_tolerance_min
    if len(diffs) == 2:
        return bool(clear_drop[-1])
    return bool(clear_drop[-1] or np.sum(clear_drop) >= 2)


def _has_clear_negative_short_slope(points: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> bool:
    ordered = _median_points_by_x(points)
    if len(ordered) != 2:
        return False
    x = ordered["x碳数"].to_numpy(dtype=float)
    y = ordered["归一化保留时间"].to_numpy(dtype=float)
    if x[1] <= x[0]:
        return False
    return bool(y[1] + config.iup_tolerance_min < y[0])


def _base_estimator(fit_type: str) -> object:
    if fit_type == "Linear":
        return LinearRegression()
    return make_pipeline(
        PolynomialFeatures(2, include_bias=False),
        LinearRegression(),
    )


def _direct_model_params(model: object, fit_type: str) -> list[float]:
    if fit_type == "Linear":
        return [float(model.coef_[0]), float(model.intercept_)]
    linear = model.named_steps["linearregression"]
    return [float(linear.coef_[1]), float(linear.coef_[0]), float(linear.intercept_)]


def iterative_refit(
    points: pd.DataFrame,
    fit_type: str,
    config: RTIUPConfig = RTIUPConfig(),
) -> tuple[object, pd.DataFrame] | None:
    """RANSAC-seed a fit, remove residual outliers, and refit to stability.

    Input points are first represented by one median RT per carbon number.  A
    final ordinary least-squares fit is always performed on the stable inlier
    set, so serialized parameters never come from an obsolete RANSAC sample.
    """

    active = _median_points_by_x(points)
    min_distinct_x = 3 if fit_type == "Linear" else config.quadratic_min_distinct_x
    min_samples = 2 if fit_type == "Linear" else 3
    if len(active) < min_distinct_x:
        return None

    x = active["x碳数"].to_numpy(dtype=float).reshape(-1, 1)
    y = active["归一化保留时间"].to_numpy(dtype=float)
    ransac = RANSACRegressor(
        estimator=_base_estimator(fit_type),
        min_samples=min_samples,
        residual_threshold=config.rt_window_min,
        random_state=config.random_state,
        max_trials=300,
    )
    try:
        ransac.fit(x, y)
    except Exception:
        return None
    seed_mask = np.asarray(ransac.inlier_mask_, dtype=bool)
    if int(seed_mask.sum()) < min_distinct_x:
        return None
    active = active.loc[seed_mask].copy()

    for _ in range(config.max_refit_iterations):
        model = _base_estimator(fit_type)
        x_active = active["x碳数"].to_numpy(dtype=float).reshape(-1, 1)
        y_active = active["归一化保留时间"].to_numpy(dtype=float)
        model.fit(x_active, y_active)
        residual = np.abs(y_active - model.predict(x_active))
        new_active = active.loc[residual <= config.rt_window_min].copy()
        if new_active.index.equals(active.index):
            break
        if len(new_active) < min_distinct_x:
            return None
        active = new_active

    final_model = _base_estimator(fit_type)
    final_model.fit(
        active["x碳数"].to_numpy(dtype=float).reshape(-1, 1),
        active["归一化保留时间"].to_numpy(dtype=float),
    )
    return final_model, active


def _fit_candidate_model(points: pd.DataFrame, fit_type: str, config: RTIUPConfig) -> dict[str, object] | None:
    result = iterative_refit(points, fit_type, config)
    if result is None:
        return None
    model, inlier_points = result
    x = inlier_points["x碳数"].to_numpy(dtype=float).reshape(-1, 1)
    y = inlier_points["归一化保留时间"].to_numpy(dtype=float)
    pred = model.predict(x)
    r2 = float(r2_score(y, pred))
    params = _direct_model_params(model, fit_type)
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    if not math.isfinite(r2) or r2 < config.r2_threshold:
        return None
    if not _is_positive_monotonic(fit_type, params, x_min, x_max, config):
        return None
    if fit_type == "Quadratic" and _has_clear_decreasing_tail(inlier_points, config):
        return None

    return {
        "拟合类型": fit_type,
        "参数": params,
        "R²": r2,
        "代表点数": int(len(_median_points_by_x(points))),
        "内点数": int(len(inlier_points)),
        "x最小": x_min,
        "x最大": x_max,
    }


def fit_line(group: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> dict[str, object]:
    """Fit one unsaturation curve using linear/quadratic RANSAC candidates."""

    plot_group = str(group["细类"].iloc[0])
    line_unsat = float(group["曲线不饱和度"].iloc[0])
    raw_points = (
        group[["x碳数", "归一化保留时间"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["x碳数", "归一化保留时间"])
        .reset_index(drop=True)
    )
    fit_points = _median_points_by_x(raw_points)
    base = {
        "细类": plot_group,
        "不饱和度": line_unsat,
        "拟合成功": False,
        "拟合类型": None,
        "参数": None,
        "R²": None,
        "点数": int(len(raw_points)),
        "不同碳数": int(len(fit_points)),
        "代表点数": int(len(fit_points)),
        "内点数": 0,
        "x最小": None,
        "x最大": None,
        "失败原因": "",
        "是否有效曲线": False,
        "去除原因": "",
    }
    if len(fit_points) < 3:
        base["失败原因"] = "少于3个点"
        return base

    linear_fit = _fit_candidate_model(fit_points, "Linear", config)
    quadratic_fit = None
    if len(fit_points) >= config.quadratic_min_distinct_x:
        quadratic_fit = _fit_candidate_model(fit_points, "Quadratic", config)
    if linear_fit is None and quadratic_fit is None:
        base["失败原因"] = "未达到R²或单调ECN要求"
        return base
    if linear_fit is None:
        best = quadratic_fit
    elif quadratic_fit is None:
        best = linear_fit
    elif float(quadratic_fit["R²"]) - float(linear_fit["R²"]) < config.quadratic_min_r2_gain:
        best = linear_fit
    else:
        best = quadratic_fit
    assert best is not None
    best["点数"] = int(len(raw_points))
    best["不同碳数"] = int(len(fit_points))
    base.update(best)
    base["拟合成功"] = True
    return base


def derivative_values(line: pd.Series, x: np.ndarray | float) -> np.ndarray:
    """Return the first derivative of a serialized linear/quadratic fit."""

    params = line.get("参数")
    if isinstance(params, str):
        params = json.loads(params)
    arr = np.asarray(x, dtype=float)
    if line.get("拟合类型") == "Linear" and isinstance(params, (list, tuple)) and len(params) == 2:
        return np.full_like(arr, float(params[0]), dtype=float)
    if line.get("拟合类型") == "Quadratic" and isinstance(params, (list, tuple)) and len(params) == 3:
        return 2.0 * float(params[0]) * arr + float(params[1])
    return np.full_like(arr, np.nan, dtype=float)


def pair_parallel_status(
    left: pd.Series,
    right: pd.Series,
    config: RTIUPConfig = RTIUPConfig(),
) -> tuple[bool | None, dict[str, object]]:
    """Evaluate slope parallelism over the shared carbon-number interval."""

    left_min, left_max = finite(left.get("x最小")), finite(left.get("x最大"))
    right_min, right_max = finite(right.get("x最小")), finite(right.get("x最大"))
    if None in (left_min, left_max, right_min, right_max):
        return None, {"平行性说明": "缺少拟合范围，无法判断"}
    start = max(float(left_min), float(right_min))
    end = min(float(left_max), float(right_max))
    if end - start < config.parallel_min_overlap_x:
        return None, {"平行性说明": "重叠碳数范围不足，无法判断"}

    x_grid = np.linspace(start, end, 80)
    left_slope = derivative_values(left, x_grid)
    right_slope = derivative_values(right, x_grid)
    valid = np.isfinite(left_slope) & np.isfinite(right_slope)
    if not valid.any():
        return None, {"平行性说明": "拟合导数不可用，无法判断"}
    left_slope = left_slope[valid]
    right_slope = right_slope[valid]
    abs_diff = float(np.median(np.abs(left_slope - right_slope)))
    reference = float(np.median(np.maximum(np.abs(left_slope), np.abs(right_slope))))
    relative_diff = abs_diff / max(reference, 1e-8)
    passed = bool(
        abs_diff <= config.parallel_abs_slope_tolerance
        or relative_diff <= config.parallel_relative_slope_tolerance
    )
    return passed, {
        "斜率绝对差": abs_diff,
        "斜率相对差": relative_diff,
        "平行性说明": "通过" if passed else "相邻不饱和度曲线斜率差超过阈值",
    }


def apply_parallel_checks(
    lines: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    *,
    invalidate: bool = False,
) -> pd.DataFrame:
    """Annotate adjacent-unsaturation pairs, optionally invalidating failures.

    The higher-unsaturation curve is compared with the immediately lower
    available curve in the same class. ECN uses this as audit metadata only;
    IUP may use it while searching/ranking alternatives. Insufficient overlap
    is recorded as unjudgeable and never removes a curve.
    """

    checked = lines.copy()
    defaults: dict[str, object] = {
        "平行性是否通过": None,
        "平行性参考曲线": "",
        "斜率绝对差": np.nan,
        "斜率相对差": np.nan,
        "平行性说明": "",
    }
    for column, default in defaults.items():
        checked[column] = default

    for _, indices in checked.groupby("细类", sort=False).groups.items():
        ordered = checked.loc[indices].sort_values("不饱和度", kind="mergesort")
        successful = ordered[ordered["拟合成功"].eq(True)]
        previous_index: int | None = None
        for row_index, row in successful.iterrows():
            if previous_index is None:
                checked.at[row_index, "平行性说明"] = "最低不饱和度曲线，无下邻参考曲线"
                previous_index = row_index
                continue
            reference = checked.loc[previous_index]
            passed, detail = pair_parallel_status(reference, row, config)
            checked.at[row_index, "平行性是否通过"] = passed
            checked.at[row_index, "平行性参考曲线"] = str(reference.get("不饱和度"))
            for column, value in detail.items():
                checked.at[row_index, column] = value
            if invalidate and passed is False:
                checked.at[row_index, "是否有效曲线"] = False
                checked.at[row_index, "去除原因"] = "相邻不饱和度曲线未通过平行性检查"
            previous_index = row_index
    return checked


def line_y_range(line: pd.Series) -> float | None:
    x_min = finite(line.get("x最小"))
    x_max = finite(line.get("x最大"))
    if x_min is None or x_max is None or x_max <= x_min:
        return None
    x_grid = np.linspace(x_min, x_max, 80)
    y_grid = fit_prediction(line, x_grid)
    if y_grid is None:
        return None
    y_grid = np.asarray(y_grid, dtype=float)
    if not np.all(np.isfinite(y_grid)):
        return None
    return float(np.max(y_grid) - np.min(y_grid))


def is_near_horizontal(line: pd.Series, config: RTIUPConfig = RTIUPConfig()) -> bool:
    x_min = finite(line.get("x最小"))
    x_max = finite(line.get("x最大"))
    if x_min is None or x_max is None or x_max <= x_min:
        return True
    x_span = x_max - x_min
    y_span = line_y_range(line)
    if y_span is None:
        return True
    if y_span < config.flat_y_range_min:
        return True
    return bool(x_span >= 3 and (y_span / x_span) < config.flat_slope_min)


def pair_has_obvious_order_violation(lower: pd.Series, higher: pd.Series, config: RTIUPConfig = RTIUPConfig()) -> bool:
    low_unsat = finite(lower.get("不饱和度"))
    high_unsat = finite(higher.get("不饱和度"))
    if low_unsat is None or high_unsat is None or low_unsat >= high_unsat:
        return False

    low_min, low_max = finite(lower.get("x最小")), finite(lower.get("x最大"))
    high_min, high_max = finite(higher.get("x最小")), finite(higher.get("x最大"))
    if None in (low_min, low_max, high_min, high_max):
        return False
    start = max(float(low_min), float(high_min))
    end = min(float(low_max), float(high_max))
    if end - start < 0.5:
        return False

    x_grid = np.linspace(start, end, 80)
    low_y = fit_prediction(lower, x_grid)
    high_y = fit_prediction(higher, x_grid)
    if low_y is None or high_y is None:
        return False
    diff = np.asarray(low_y, dtype=float) - np.asarray(high_y, dtype=float)
    diff = diff[np.isfinite(diff)]
    if len(diff) == 0:
        return False
    # This is a pre-correction feasibility test. Each independently fitted
    # series may still move by ±rt_window_min, and the locked IUP tolerance is
    # added only when comparing the two series. A curve is removed here only
    # when even those bounded opposite shifts cannot repair the order.
    tolerance = relative_rt_tolerance_min(config)
    bad_fraction = float(np.mean(diff < -tolerance))
    return bool(
        bad_fraction >= 0.75
        and float(np.median(diff)) < -tolerance
        and float(np.mean(diff)) < -tolerance
    )


def subset_is_ordered(lines: Iterable[pd.Series], config: RTIUPConfig = RTIUPConfig()) -> bool:
    ordered = sorted(lines, key=lambda row: float(row["不饱和度"]))
    for left, right in itertools.combinations(ordered, 2):
        if pair_has_obvious_order_violation(left, right, config):
            return False
    return True


def subset_score(lines: Iterable[pd.Series]) -> float:
    score = 0.0
    for line in lines:
        r2 = finite(line.get("R²")) or 0.0
        n_inliers = finite(line.get("内点数")) or 0.0
        span = line_y_range(line) or 0.0
        score += max(0.0, min(1.0, r2)) * 3.0 + math.log1p(max(0.0, n_inliers)) + span * 0.05
    return score


def ensure_rescue_columns(lines: pd.DataFrame) -> pd.DataFrame:
    lines = lines.copy()
    defaults = {
        "二次捞点": False,
        "救回点数": 0,
        "救回不同碳数": 0,
        "救回下界不饱和度": None,
        "救回上界不饱和度": None,
        "救回说明": "",
    }
    for column, default in defaults.items():
        if column not in lines.columns:
            lines[column] = default
    return lines


def choose_iup_lines(lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> pd.DataFrame:
    """Keep the largest high-scoring subset that satisfies IUP ordering."""

    lines = ensure_rescue_columns(lines)
    if lines.empty:
        return lines

    valid_mask = lines["拟合成功"].eq(True)
    lines.loc[valid_mask, "是否有效曲线"] = True
    rows = [row for _, row in lines[valid_mask].iterrows()]
    flat_keys: set[tuple[str, float]] = set()
    for row in rows:
        key = (str(row["细类"]), float(row["不饱和度"]))
        if is_near_horizontal(row, config):
            flat_keys.add(key)
            mask = lines["细类"].eq(key[0]) & lines["不饱和度"].astype(float).eq(key[1])
            lines.loc[mask, "是否有效曲线"] = False
            lines.loc[mask, "去除原因"] = "近似横线，去除"

    candidates = [row for row in rows if (str(row["细类"]), float(row["不饱和度"])) not in flat_keys]
    if len(candidates) <= 1:
        keep_keys = {(str(row["细类"]), float(row["不饱和度"])) for row in candidates}
    elif len(candidates) <= 12:
        best: tuple[pd.Series, ...] | None = None
        best_score = -1.0
        for size in range(len(candidates), 0, -1):
            for combo in itertools.combinations(candidates, size):
                combo_list = list(combo)
                if not subset_is_ordered(combo_list, config):
                    continue
                current_score = subset_score(combo_list)
                if current_score > best_score:
                    best = combo
                    best_score = current_score
            if best is not None:
                break
        keep_keys = {(str(row["细类"]), float(row["不饱和度"])) for row in (best or [])}
    else:
        remaining = candidates[:]
        changed = True
        while changed:
            changed = False
            for left, right in itertools.combinations(remaining, 2):
                if pair_has_obvious_order_violation(left, right, config):
                    remove = left if subset_score([left]) < subset_score([right]) else right
                    remaining = [row for row in remaining if row["不饱和度"] != remove["不饱和度"]]
                    changed = True
                    break
        keep_keys = {(str(row["细类"]), float(row["不饱和度"])) for row in remaining}

    for row in candidates:
        key = (str(row["细类"]), float(row["不饱和度"]))
        if key not in keep_keys:
            mask = lines["细类"].eq(key[0]) & lines["不饱和度"].astype(float).eq(key[1])
            lines.loc[mask, "是否有效曲线"] = False
            lines.loc[mask, "去除原因"] = "明显违反IUP上下顺序，去除"
    return lines


def _iup_candidate_key(line: pd.Series) -> tuple[object, ...]:
    params = line.get("参数")
    if isinstance(params, str):
        params = json.loads(params)
    rounded = tuple(round(float(value), 5) for value in (params or []))
    return (line.get("拟合类型"), rounded, finite(line.get("x最小")), finite(line.get("x最大")))


def generate_iup_line_candidates(
    group: pd.DataFrame,
    initial_line: pd.Series,
    config: RTIUPConfig = RTIUPConfig(),
    *,
    multistart: bool = True,
) -> list[pd.Series]:
    """Generate alternative real-data fits for one IUP unsaturation series.

    ECN deliberately uses one median per carbon number. IUP additionally runs
    deterministic multi-start RANSAC on the original RT candidates so a
    coherent sub-series is not hidden by the per-carbon median of several
    structural candidates. Every returned alternative still needs at least
    three carbon numbers, R² >= the locked threshold, positive monotonicity,
    and point-to-line residuals inside ``rt_window_min``.
    """

    plot_group = str(group["细类"].iloc[0])
    unsat = float(group["曲线不饱和度"].iloc[0])
    raw_points = (
        group[["x碳数", "归一化保留时间"]]
        .dropna()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    records: list[pd.Series] = []
    seen: set[tuple[object, ...]] = set()

    def append_candidate(payload: dict[str, object] | pd.Series, source: str, rescued: bool) -> None:
        record = pd.Series(payload).copy()
        if "拟合成功" in record.index and not bool_value(record.get("拟合成功")):
            return
        record["细类"] = plot_group
        record["不饱和度"] = unsat
        record["拟合成功"] = True
        record["是否有效曲线"] = True
        record["失败原因"] = ""
        record["去除原因"] = ""
        record["点数"] = int(len(raw_points))
        record["不同碳数"] = int(raw_points["x碳数"].nunique())
        record["IUP候选来源"] = source
        record["IUP重新寻优"] = rescued
        # This is a full raw-candidate re-fit, not the older bracket-only
        # 二次捞点 path; keeping the flags separate avoids applying bracket
        # metadata guards to a globally optimized curve.
        record["二次捞点"] = False
        if rescued:
            record["救回点数"] = int(record.get("内点数") or 0)
            record["救回不同碳数"] = int(record.get("内点数") or 0)
            record["救回说明"] = "IUP多起点原始候选寻优：阈值内重新选择同碳数代表点"
        key = _iup_candidate_key(record)
        if key not in seen:
            seen.add(key)
            records.append(record)

    append_candidate(initial_line, "同碳数中位数初始拟合", False)
    if multistart and len(raw_points) >= 3 and raw_points["x碳数"].nunique() >= 3:
        carbon_groups = [part.reset_index(drop=True) for _, part in raw_points.groupby("x碳数", sort=True)]
        branch_count = math.prod(len(part) for part in carbon_groups)
        if branch_count <= 256:
            for choices in itertools.product(*(range(len(part)) for part in carbon_groups)):
                branch = pd.DataFrame(
                    [part.iloc[choice] for part, choice in zip(carbon_groups, choices)]
                ).reset_index(drop=True)
                branch_candidate = _fit_candidate_model(branch, "Linear", config)
                if branch_candidate is not None:
                    append_candidate(branch_candidate, "逐碳数分支穷举RANSAC", True)

        x_all = raw_points["x碳数"].to_numpy(dtype=float).reshape(-1, 1)
        y_all = raw_points["归一化保留时间"].to_numpy(dtype=float)
        for seed in range(max(0, int(config.iup_multistart_trials))):
            # Also draw one observed candidate per carbon. This explicitly
            # explores alternative structural branches when two or more RTs
            # share the same total carbon number; raw-point RANSAC alone can
            # repeatedly sample duplicate x values and miss that branch.
            rng = np.random.default_rng(config.random_state + seed)
            sampled_rows: list[pd.Series] = []
            for _, carbon_group in raw_points.groupby("x碳数", sort=True):
                sampled_rows.append(carbon_group.iloc[int(rng.integers(0, len(carbon_group)))])
            sampled = pd.DataFrame(sampled_rows).reset_index(drop=True)
            sampled_candidate = _fit_candidate_model(sampled, "Linear", config)
            if sampled_candidate is not None:
                append_candidate(sampled_candidate, "逐碳数分支抽样RANSAC", True)

            ransac = RANSACRegressor(
                estimator=LinearRegression(),
                min_samples=3,
                residual_threshold=config.rt_window_min,
                random_state=config.random_state + seed,
                max_trials=300,
            )
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UndefinedMetricWarning)
                    ransac.fit(x_all, y_all)
            except Exception:
                continue
            mask = np.asarray(ransac.inlier_mask_, dtype=bool)
            if int(mask.sum()) < 3:
                continue
            representatives = _median_points_by_x(raw_points.loc[mask])
            candidate = _fit_candidate_model(representatives, "Linear", config)
            if candidate is not None:
                append_candidate(candidate, "多起点原始候选RANSAC", True)

    def rank(record: pd.Series) -> tuple[int, float, int]:
        linear_bonus = 1 if record.get("拟合类型") == "Linear" else 0
        return (
            int(finite(record.get("内点数")) or 0),
            float(finite(record.get("R²")) or 0.0),
            linear_bonus,
        )

    records.sort(key=rank, reverse=True)
    return records[: max(1, int(config.iup_max_candidates_per_unsaturation))]


def solve_iup_vertical_shifts(
    lines: list[pd.Series],
    config: RTIUPConfig = RTIUPConfig(),
) -> list[float] | None:
    """Find the smallest whole-series shifts satisfying IUP order.

    Each fitted series gets one constant offset bounded by ±0.20 min. The
    slopes and within-series residuals are untouched. After correction, every
    overlapping lower-unsaturation curve must remain above the corresponding
    higher-unsaturation curve. Consecutive unsaturations target the locked
    0.05 min IUP separation so they do not plot on top of one another; across
    a missing-unsaturation gap the locked 0.08 min order tolerance is used.
    ``None`` means no honest bounded correction exists.
    """

    n = len(lines)
    if n == 0:
        return []
    # Variables are offsets followed by their absolute-value auxiliaries.
    objective = np.r_[np.zeros(n), np.ones(n)]
    a_ub: list[list[float]] = []
    b_ub: list[float] = []
    for left_index, left in enumerate(lines):
        for right_index in range(left_index + 1, n):
            right = lines[right_index]
            start = max(float(left["x最小"]), float(right["x最小"]))
            end = min(float(left["x最大"]), float(right["x最大"]))
            if end - start < 0.5:
                continue
            x_grid = np.linspace(start, end, 80)
            left_y = fit_prediction(left, x_grid)
            right_y = fit_prediction(right, x_grid)
            if left_y is None or right_y is None:
                continue
            left_unsat = float(left["不饱和度"])
            right_unsat = float(right["不饱和度"])
            # Consecutive unsaturations should be visibly separated by the
            # locked 0.05 min IUP amount. Across a larger missing-unsaturation
            # gap, enforce non-inversion with the locked 0.08 min tolerance but
            # do not manufacture equal spacing through an unobserved series.
            required_gap = (
                config.iup_tolerance_min
                if abs(right_unsat - left_unsat) <= 1.01
                else -config.order_tolerance_min
            )
            rhs = float(
                np.min(
                    np.asarray(left_y, dtype=float)
                    - np.asarray(right_y, dtype=float)
                    - required_gap
                )
            )
            constraint = [0.0] * (2 * n)
            constraint[right_index] = 1.0
            constraint[left_index] = -1.0
            a_ub.append(constraint)
            b_ub.append(rhs)

    for index in range(n):
        positive = [0.0] * (2 * n)
        positive[index] = 1.0
        positive[n + index] = -1.0
        a_ub.append(positive)
        b_ub.append(0.0)
        negative = [0.0] * (2 * n)
        negative[index] = -1.0
        negative[n + index] = -1.0
        a_ub.append(negative)
        b_ub.append(0.0)

    result = optimize.linprog(
        objective,
        A_ub=np.asarray(a_ub, dtype=float),
        b_ub=np.asarray(b_ub, dtype=float),
        bounds=[(-config.rt_window_min, config.rt_window_min)] * n
        + [(0.0, config.rt_window_min)] * n,
        method="highs",
    )
    if not result.success:
        return None
    return [float(value) for value in result.x[:n]]


def _iup_candidate_set_score(
    lines: list[pd.Series],
    shifts: list[float],
    config: RTIUPConfig,
) -> tuple[float, ...]:
    parallel_passes = 0
    parallel_penalty = 0.0
    for left, right in zip(lines, lines[1:]):
        passed, detail = pair_parallel_status(left, right, config)
        parallel_passes += int(passed is True)
        relative = finite(detail.get("斜率相对差"))
        parallel_penalty += relative if relative is not None else 2.0
    return (
        float(len(lines)),
        float(sum(int(finite(line.get("内点数")) or 0) for line in lines)),
        float(parallel_passes),
        -parallel_penalty,
        float(sum(finite(line.get("R²")) or 0.0 for line in lines)),
        float(sum(line.get("拟合类型") == "Linear" for line in lines)),
        -float(sum(abs(value) for value in shifts)),
    )


def optimize_iup_lines_for_group(
    source: pd.DataFrame,
    initial_lines: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    *,
    multistart: bool = True,
) -> pd.DataFrame:
    """Choose the largest bounded, ordered, preferentially parallel IUP set."""

    candidate_groups: list[tuple[float, list[pd.Series]]] = []
    for unsat, group in source.groupby("曲线不饱和度", sort=True):
        match = initial_lines[initial_lines["不饱和度"].astype(float).eq(float(unsat))]
        if match.empty:
            continue
        candidates = generate_iup_line_candidates(
            group,
            match.iloc[0],
            config,
            multistart=multistart,
        )
        candidate_groups.append((float(unsat), candidates))

    beam: list[tuple[list[pd.Series], list[float]]] = [([], [])]
    for _, candidates in candidate_groups:
        skipped_states = list(beam)
        next_beam = list(skipped_states)  # Skipping a genuinely infeasible unsaturation is allowed.
        for selected, _ in beam:
            for candidate in candidates:
                proposed = selected + [candidate]
                shifts = solve_iup_vertical_shifts(proposed, config)
                if shifts is not None:
                    next_beam.append((proposed, shifts))
        next_beam.sort(
            key=lambda item: _iup_candidate_set_score(item[0], item[1], config),
            reverse=True,
        )
        width = max(1, int(config.iup_beam_width))
        reserve = min(len(skipped_states), max(1, width // 8))
        skipped_states.sort(
            key=lambda item: _iup_candidate_set_score(item[0], item[1], config),
            reverse=True,
        )
        combined = next_beam[: max(0, width - reserve)] + skipped_states[:reserve]
        deduplicated: list[tuple[list[pd.Series], list[float]]] = []
        seen_states: set[tuple[object, ...]] = set()
        for state in combined:
            signature = tuple(
                (float(line["不饱和度"]), _iup_candidate_key(line))
                for line in state[0]
            )
            if signature in seen_states:
                continue
            seen_states.add(signature)
            deduplicated.append(state)
        beam = deduplicated[:width]

    selected, shifts = max(
        beam,
        key=lambda item: _iup_candidate_set_score(item[0], item[1], config),
    )
    optimized = ensure_rescue_columns(initial_lines)
    optimized["是否有效曲线"] = False
    successful = optimized["拟合成功"].eq(True)
    optimized.loc[successful, "去除原因"] = "IUP全局寻优未选中：固定漂移范围内无法与最终曲线组同时满足顺序"
    for line, shift in zip(selected, shifts):
        record = line.to_dict()
        raw_params = record.get("参数")
        if isinstance(raw_params, str):
            raw_params = json.loads(raw_params)
        shifted_params = list(raw_params or [])
        if shifted_params:
            shifted_params[-1] = float(shifted_params[-1]) + float(shift)
        record["IUP原始参数"] = list(raw_params or [])
        record["IUP_RT漂移校正(min)"] = float(shift)
        record["参数"] = shifted_params
        record["是否有效曲线"] = True
        record["去除原因"] = ""
        mask = optimized["细类"].eq(record["细类"]) & optimized["不饱和度"].astype(float).eq(float(record["不饱和度"]))
        if not mask.any():
            optimized = pd.concat([optimized, pd.DataFrame([record])], ignore_index=True)
            continue
        row_index = optimized.index[mask][0]
        for column, value in record.items():
            if column not in optimized.columns:
                optimized[column] = None
            optimized.at[row_index, column] = value
    if "IUP_RT漂移校正(min)" not in optimized.columns:
        optimized["IUP_RT漂移校正(min)"] = 0.0
    optimized["IUP_RT漂移校正(min)"] = pd.to_numeric(
        optimized["IUP_RT漂移校正(min)"], errors="coerce"
    ).fillna(0.0)
    return optimized


def point_iup_status(point: pd.Series, valid_lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> tuple[bool, str]:
    """Return whether a 1-2 point series can be shown by IUP order alone."""

    x = finite(point.get("x碳数"))
    y = finite(point.get("归一化保留时间"))
    unsat = finite(point.get("曲线不饱和度"))
    if x is None or y is None or unsat is None:
        return False, "缺少x/保留时间/不饱和度"

    tolerance = relative_rt_tolerance_min(config)
    supports: list[str] = []
    violations: list[str] = []
    for _, line in valid_lines.iterrows():
        other_unsat = finite(line.get("不饱和度"))
        if other_unsat is None or other_unsat == unsat:
            continue
        x_min, x_max = finite(line.get("x最小")), finite(line.get("x最大"))
        if x_min is not None and x_max is not None and not (x_min <= x <= x_max):
            continue
        pred = finite(fit_prediction(line, x))
        if pred is None:
            continue
        if other_unsat < unsat:
            if pred + tolerance < y:
                violations.append(f"低不饱和度{other_unsat:g}曲线未在上方")
            else:
                supports.append(f"低不饱和度{other_unsat:g}曲线在上方")
        else:
            if y + tolerance < pred:
                violations.append(f"高不饱和度{other_unsat:g}曲线未在下方")
            else:
                supports.append(f"高不饱和度{other_unsat:g}曲线在下方")
    if violations:
        return False, "IUP不满足：" + "；".join(violations[:3])
    if supports:
        return True, "IUP满足：" + "；".join(supports[:3])
    return False, "无相邻有效曲线可判断IUP"


def apply_short_series_guard(rows: pd.DataFrame, lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> pd.DataFrame:
    """Remove short 1-2 point series with clear negative slope or IUP violations."""

    if "点数不足保留" not in rows.columns:
        return rows
    guarded = rows.copy()
    if "点数不足过滤原因" not in guarded.columns:
        guarded["点数不足过滤原因"] = ""
    else:
        guarded["点数不足过滤原因"] = guarded["点数不足过滤原因"].fillna("")
    short_rows = guarded[guarded["点数不足保留"].eq(True)]
    if short_rows.empty:
        return guarded

    for (plot_group, unsat), group in short_rows.groupby(["细类", "曲线不饱和度"], sort=False):
        idx = group.index
        notes = group.get("IUP说明", pd.Series("", index=idx)).astype(str)
        negative_slope = _has_clear_negative_short_slope(group, config)
        explicit_violation = notes.str.startswith("IUP不满足").any()
        if negative_slope:
            reason = "点数不足曲线负斜率，去除"
        elif explicit_violation:
            reason = "点数不足曲线明显违反IUP，去除"
        else:
            continue
        guarded.loc[idx, "点数不足保留"] = False
        guarded.loc[idx, "是否作图"] = False
        guarded.loc[idx, "IUP说明"] = reason
        guarded.loc[idx, "点数不足过滤原因"] = reason
    return guarded


def find_bracket_lines(valid_lines: pd.DataFrame, unsat: float) -> tuple[pd.Series | None, pd.Series | None]:
    valid = valid_lines.copy()
    valid["_unsat"] = pd.to_numeric(valid["不饱和度"], errors="coerce")
    lower = valid[valid["_unsat"].lt(unsat)].sort_values("_unsat")
    higher = valid[valid["_unsat"].gt(unsat)].sort_values("_unsat")
    if lower.empty or higher.empty:
        return None, None
    return lower.iloc[-1].drop(labels=["_unsat"], errors="ignore"), higher.iloc[0].drop(labels=["_unsat"], errors="ignore")


def point_between_bracket(
    point: pd.Series,
    lower_line: pd.Series,
    higher_line: pd.Series,
    config: RTIUPConfig = RTIUPConfig(),
) -> tuple[bool, str, float | None]:
    x = finite(point.get("x碳数"))
    y = finite(point.get("归一化保留时间"))
    if x is None or y is None:
        return False, "缺少x或保留时间", None

    low_x_min, low_x_max = finite(lower_line.get("x最小")), finite(lower_line.get("x最大"))
    high_x_min, high_x_max = finite(higher_line.get("x最小")), finite(higher_line.get("x最大"))
    if None in (low_x_min, low_x_max, high_x_min, high_x_max):
        return False, "相邻曲线缺少x范围", None
    if not (float(low_x_min) <= x <= float(low_x_max) and float(high_x_min) <= x <= float(high_x_max)):
        return False, "x不在相邻有效曲线重叠范围", None

    upper_y = finite(fit_prediction(lower_line, x))
    lower_y = finite(fit_prediction(higher_line, x))
    if upper_y is None or lower_y is None:
        return False, "相邻曲线无法预测", None
    tolerance = relative_rt_tolerance_min(config)
    if float(upper_y) + tolerance < float(lower_y):
        return False, "相邻曲线在该x处交叉", None
    if y > float(upper_y) + tolerance:
        return False, "高于低不饱和度相邻曲线", None
    if y < float(lower_y) - tolerance:
        return False, "低于高不饱和度相邻曲线", None

    low_unsat = finite(lower_line.get("不饱和度"))
    high_unsat = finite(higher_line.get("不饱和度"))
    target_unsat = finite(point.get("曲线不饱和度"))
    if low_unsat is None or high_unsat is None or target_unsat is None or high_unsat <= low_unsat:
        expected_y = (float(upper_y) + float(lower_y)) / 2.0
    else:
        fraction = (target_unsat - low_unsat) / (high_unsat - low_unsat)
        expected_y = float(upper_y) - fraction * (float(upper_y) - float(lower_y))
    return True, f"位于{low_unsat:g}-{high_unsat:g}不饱和度曲线之间", abs(y - expected_y)


def crossed_unsaturation_count(lower_unsat: object, target_unsat: object, higher_unsat: object) -> int | None:
    """Count skipped integer unsaturation values inside a bracket, excluding the target."""

    low = finite(lower_unsat)
    target = finite(target_unsat)
    high = finite(higher_unsat)
    if low is None or target is None or high is None or not (low < target < high):
        return None
    low_i = int(round(low))
    target_i = int(round(target))
    high_i = int(round(high))
    return max(0, target_i - low_i - 1) + max(0, high_i - target_i - 1)


def two_point_series_iup_status(
    points: pd.DataFrame,
    valid_lines: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
) -> tuple[bool, str]:
    """Validate a short two-point connector against a nearby IUP bracket."""

    if points.empty:
        return False, "短两点线缺少候选点"
    target_unsat = finite(points["曲线不饱和度"].iloc[0])
    if target_unsat is None:
        return False, "短两点线缺少不饱和度"
    lower_line, higher_line = find_bracket_lines(valid_lines, target_unsat)
    if lower_line is None or higher_line is None:
        return False, "短两点线过滤：缺少上下相邻有效曲线"

    lower_unsat = finite(lower_line.get("不饱和度"))
    higher_unsat = finite(higher_line.get("不饱和度"))
    crossed = crossed_unsaturation_count(lower_unsat, target_unsat, higher_unsat)
    if crossed is None:
        return False, "短两点线过滤：无法计算上下相邻不饱和度跨度"
    if crossed > config.short_series_max_crossed_unsats:
        return (
            False,
            f"短两点线过滤：横跨{crossed}个不饱和度，超过最多{config.short_series_max_crossed_unsats}个",
        )

    unique_points = (
        points[["x碳数", "归一化保留时间", "曲线不饱和度"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["x碳数", "归一化保留时间"])
    )
    if len(unique_points) != 2:
        return False, f"短两点线过滤：唯一散点数为{len(unique_points)}，不是2"

    notes: list[str] = []
    for _, point in unique_points.iterrows():
        ok, note, _ = point_between_bracket(point, lower_line, higher_line, config)
        notes.append(note)
        if not ok:
            return False, "短两点线IUP不满足：" + note
    return True, f"短两点线IUP满足：位于{lower_unsat:g}-{higher_unsat:g}不饱和度曲线之间，横跨{crossed}个不饱和度"


def _is_judgeable_iup_point(note: str) -> bool:
    skipped = (
        "缺少x或保留时间",
        "相邻曲线缺少x范围",
        "x不在相邻有效曲线重叠范围",
        "相邻曲线无法预测",
        "相邻曲线在该x处交叉",
        "无相邻有效曲线可判断IUP",
    )
    return not any(text in note for text in skipped)


def point_level_iup_status(point: pd.Series, valid_lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> tuple[bool, str]:
    """Check one fitted inlier point against adjacent or one-sided IUP lines."""

    unsat = finite(point.get("曲线不饱和度"))
    if unsat is None:
        return False, "缺少不饱和度"
    lower_line, higher_line = find_bracket_lines(valid_lines, unsat)
    if lower_line is not None and higher_line is not None:
        ok, note, _ = point_between_bracket(point, lower_line, higher_line, config)
        if ok or _is_judgeable_iup_point(note):
            return ok, note
    return point_iup_status(point, valid_lines, config)


def apply_point_level_iup_guard(
    rows: pd.DataFrame,
    lines: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
) -> pd.DataFrame:
    """Remove only fitted inlier points that violate IUP while retaining good points."""

    if rows.empty or lines.empty or "是否有效曲线" not in lines.columns:
        return rows
    guarded = rows.copy()
    guarded["拟合点IUP支持数"] = np.nan
    guarded["拟合点IUP判断数"] = np.nan
    guarded["点级IUP是否通过"] = ""
    guarded["点级IUP说明"] = ""
    guarded["点级IUP过滤原因"] = ""

    valid_lines = lines[lines["是否有效曲线"].eq(True)].copy()
    for plot_group, group_lines in valid_lines.groupby("细类", sort=False):
        group_lines = group_lines.sort_values("不饱和度")
        for _, line in group_lines.iterrows():
            unsat = finite(line.get("不饱和度"))
            if unsat is None:
                continue
            point_mask = (
                guarded["细类"].eq(plot_group)
                & guarded["曲线不饱和度"].astype(float).eq(float(unsat))
                & guarded["同曲线0.2min内点"].eq(True)
            )
            points = guarded[point_mask].copy()
            if points.empty:
                continue

            decisions: dict[int, tuple[str, str]] = {}
            supported = 0
            judged = 0
            for row_idx, point in points.iterrows():
                ok, note = point_level_iup_status(point, group_lines, config)
                if not _is_judgeable_iup_point(note):
                    decisions[row_idx] = ("unjudged", note)
                    continue
                judged += 1
                if ok:
                    supported += 1
                    decisions[row_idx] = ("pass", note)
                else:
                    decisions[row_idx] = ("fail", note)

            guarded.loc[point_mask, "拟合点IUP支持数"] = supported
            guarded.loc[point_mask, "拟合点IUP判断数"] = judged
            for row_idx, (status, note) in decisions.items():
                guarded.at[row_idx, "点级IUP说明"] = note
                if status == "pass":
                    guarded.at[row_idx, "点级IUP是否通过"] = True
                    existing_note = str(guarded.at[row_idx, "IUP说明"] or "")
                    if "点级IUP满足" not in existing_note:
                        guarded.at[row_idx, "IUP说明"] = f"{existing_note}；点级IUP满足：{note}" if existing_note else f"点级IUP满足：{note}"
                elif status == "fail":
                    reason = "点级IUP不满足：" + note
                    guarded.at[row_idx, "点级IUP是否通过"] = False
                    guarded.at[row_idx, "点级IUP过滤原因"] = reason
                    guarded.at[row_idx, "同曲线0.2min内点"] = False
                    guarded.at[row_idx, "是否作图"] = False
                    guarded.at[row_idx, "IUP说明"] = reason
                    if "去除原因" in guarded.columns:
                        guarded.at[row_idx, "去除原因"] = reason
                else:
                    guarded.at[row_idx, "点级IUP是否通过"] = ""
    return guarded


def apply_fitted_point_iup_guard(
    rows: pd.DataFrame,
    lines: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    min_supported_points: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop fitted curves unless at least two inlier points satisfy bracket IUP."""

    if rows.empty or lines.empty or "是否有效曲线" not in lines.columns:
        return rows, lines

    guarded_rows = rows.copy()
    guarded_lines = lines.copy()
    guarded_rows["拟合点IUP支持数"] = np.nan
    guarded_rows["拟合点IUP判断数"] = np.nan
    guarded_lines["拟合点IUP支持数"] = np.nan
    guarded_lines["拟合点IUP判断数"] = np.nan

    changed = True
    while changed:
        changed = False
        valid_lines = guarded_lines[guarded_lines["是否有效曲线"].eq(True)].copy()
        for plot_group, group_lines in valid_lines.groupby("细类", sort=False):
            group_lines = group_lines.sort_values("不饱和度")
            for line_idx, line in group_lines.iterrows():
                unsat = finite(line.get("不饱和度"))
                if unsat is None:
                    continue
                lower_line, higher_line = find_bracket_lines(group_lines, unsat)
                if lower_line is None or higher_line is None:
                    continue

                point_mask = (
                    guarded_rows["细类"].eq(plot_group)
                    & guarded_rows["曲线不饱和度"].astype(float).eq(float(unsat))
                    & guarded_rows["同曲线0.2min内点"].eq(True)
                )
                points = guarded_rows[point_mask].copy()
                if points.empty:
                    continue

                supported = 0
                judged = 0
                for _, point in points.iterrows():
                    ok, note, _ = point_between_bracket(point, lower_line, higher_line, config)
                    if not _is_judgeable_iup_point(note):
                        continue
                    judged += 1
                    if ok:
                        supported += 1

                guarded_rows.loc[point_mask, "拟合点IUP支持数"] = supported
                guarded_rows.loc[point_mask, "拟合点IUP判断数"] = judged
                guarded_lines.loc[line_idx, "拟合点IUP支持数"] = supported
                guarded_lines.loc[line_idx, "拟合点IUP判断数"] = judged

                if judged >= min_supported_points and supported < min_supported_points:
                    reason = f"拟合曲线IUP支持点不足：{supported}/{judged}个点位于相邻曲线之间"
                    guarded_lines.loc[line_idx, "是否有效曲线"] = False
                    guarded_lines.loc[line_idx, "去除原因"] = reason
                    guarded_rows.loc[point_mask, "是否有效曲线"] = False
                    guarded_rows.loc[point_mask, "去除原因"] = reason
                    guarded_rows.loc[point_mask, "同曲线0.2min内点"] = False
                    guarded_rows.loc[point_mask, "是否作图"] = False
                    guarded_rows.loc[point_mask, "IUP说明"] = reason
                    changed = True
                    break
            if changed:
                break
    return guarded_rows, guarded_lines


def line_between_bracket(
    line: pd.Series,
    lower_line: pd.Series,
    higher_line: pd.Series,
    config: RTIUPConfig = RTIUPConfig(),
) -> bool:
    x_min = finite(line.get("x最小"))
    x_max = finite(line.get("x最大"))
    low_x_min, low_x_max = finite(lower_line.get("x最小")), finite(lower_line.get("x最大"))
    high_x_min, high_x_max = finite(higher_line.get("x最小")), finite(higher_line.get("x最大"))
    if None in (x_min, x_max, low_x_min, low_x_max, high_x_min, high_x_max):
        return False
    start = max(float(x_min), float(low_x_min), float(high_x_min))
    end = min(float(x_max), float(low_x_max), float(high_x_max))
    if end - start < 0.5:
        return False

    x_grid = np.linspace(start, end, 100)
    y = fit_prediction(line, x_grid)
    upper_y = fit_prediction(lower_line, x_grid)
    lower_y = fit_prediction(higher_line, x_grid)
    if y is None or upper_y is None or lower_y is None:
        return False
    y = np.asarray(y, dtype=float)
    upper_y = np.asarray(upper_y, dtype=float)
    lower_y = np.asarray(lower_y, dtype=float)
    ok = np.isfinite(y) & np.isfinite(upper_y) & np.isfinite(lower_y)
    if not np.any(ok):
        return False
    tolerance = relative_rt_tolerance_min(config)
    if np.mean((upper_y[ok] + tolerance) < lower_y[ok]) > 0.10:
        return False
    inside = (y[ok] <= upper_y[ok] + tolerance) & (y[ok] >= lower_y[ok] - tolerance)
    return bool(np.mean(inside) >= 0.90)


def pick_rescue_points(points: pd.DataFrame, lower_line: pd.Series, higher_line: pd.Series, config: RTIUPConfig) -> pd.DataFrame:
    candidates = points.dropna(subset=["x碳数", "归一化保留时间", "曲线不饱和度"]).copy()
    if candidates.empty:
        return candidates
    keep_records: list[pd.Series] = []
    distances: list[float] = []
    notes: list[str] = []
    for _, point in candidates.iterrows():
        ok, note, distance = point_between_bracket(point, lower_line, higher_line, config)
        if ok:
            keep_records.append(point)
            distances.append(float(distance or 0.0))
            notes.append(note)
    if not keep_records:
        return pd.DataFrame(columns=points.columns)
    filtered = pd.DataFrame(keep_records).copy()
    filtered["_救回距离"] = distances
    filtered["_救回IUP说明"] = notes
    filtered["_总分数"] = _numeric_series(filtered, "总分数", -np.inf)
    filtered["_匹配度分数"] = _numeric_series(filtered, "匹配度分数", -np.inf)
    filtered["_丰度"] = _numeric_series(filtered, "丰度", -np.inf)
    filtered = filtered.sort_values(
        ["x碳数", "_救回距离", "_总分数", "_匹配度分数", "_丰度"],
        ascending=[True, True, False, False, False],
    )
    return filtered.groupby("x碳数", as_index=False, sort=False).head(1).reset_index(drop=True)


def _numeric_series(frame: pd.DataFrame, column: str, default: float) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def _fit_rescue_candidate(points: pd.DataFrame, fit_type: str, config: RTIUPConfig) -> dict[str, object] | None:
    # Rescue candidates receive the same RANSAC -> iterative inlier removal ->
    # final refit treatment as ordinary curves. Requiring every bracket-picked
    # point to lie on one initial OLS curve was too brittle and could discard an
    # otherwise coherent alternative because of a single local candidate.
    candidate = _fit_candidate_model(points, fit_type, config)
    if candidate is None:
        return None
    candidate["点数"] = int(len(points))
    return candidate


def build_rescue_line(
    plot_group: str,
    unsat: float,
    points: pd.DataFrame,
    lower_line: pd.Series,
    higher_line: pd.Series,
    config: RTIUPConfig = RTIUPConfig(),
) -> dict[str, object] | None:
    rescued_points = pick_rescue_points(points, lower_line, higher_line, config)
    if len(rescued_points) < config.rescue_min_distinct_x or rescued_points["x碳数"].nunique() < config.rescue_min_distinct_x:
        return None
    candidates = []
    for fit_type in ["Linear", "Quadratic"]:
        candidate = _fit_rescue_candidate(rescued_points, fit_type, config)
        if candidate is None:
            continue
        line = pd.Series(
            {
                "细类": plot_group,
                "不饱和度": unsat,
                "拟合成功": True,
                "是否有效曲线": True,
                "失败原因": "",
                "去除原因": "",
                **candidate,
            }
        )
        if is_near_horizontal(line, config):
            continue
        if not line_between_bracket(line, lower_line, higher_line, config):
            continue
        candidates.append(line.to_dict())
    if not candidates:
        return None

    def score(item: dict[str, object]) -> tuple[int, float, float, int, int]:
        candidate_line = pd.Series(item)
        parallel_results = [
            pair_parallel_status(lower_line, candidate_line, config)[0],
            pair_parallel_status(candidate_line, higher_line, config)[0],
        ]
        parallel_passes = sum(result is True for result in parallel_results)
        slope_penalty = 0.0
        for left, right in ((lower_line, candidate_line), (candidate_line, higher_line)):
            _, detail = pair_parallel_status(left, right, config)
            relative = finite(detail.get("斜率相对差"))
            if relative is not None:
                slope_penalty += relative
        complexity_bonus = 0 if item["拟合类型"] == "Linear" else -1
        return (
            parallel_passes,
            -slope_penalty,
            float(item["R²"]),
            int(item["内点数"]),
            complexity_bonus,
        )

    best = sorted(candidates, key=score, reverse=True)[0]
    lower_unsat = finite(lower_line.get("不饱和度"))
    higher_unsat = finite(higher_line.get("不饱和度"))
    best.update(
        {
            "细类": plot_group,
            "不饱和度": unsat,
            "拟合成功": True,
            "是否有效曲线": True,
            "失败原因": "",
            "去除原因": "",
            "二次捞点": True,
            "救回点数": int(len(rescued_points)),
            "救回不同碳数": int(rescued_points["x碳数"].nunique()),
            "救回下界不饱和度": lower_unsat,
            "救回上界不饱和度": higher_unsat,
            "救回说明": f"二次IUP捞点：位于{lower_unsat:g}-{higher_unsat:g}不饱和度有效曲线之间，且至少3个不同碳数",
        }
    )
    return best


def rescue_missing_iup_lines(source_points: pd.DataFrame, lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> pd.DataFrame:
    """Rescue missing unsaturation curves that fit between adjacent IUP lines."""

    lines = ensure_rescue_columns(lines)
    if source_points.empty or lines.empty:
        return lines

    selected_records: list[dict[str, object]] = []
    for plot_group, group_lines in lines.groupby("细类", sort=False):
        group_points = source_points[source_points["细类"].eq(plot_group)].copy()
        if group_points.empty:
            continue
        valid_lines = group_lines[group_lines["是否有效曲线"].eq(True)].copy()
        if len(valid_lines) < 2:
            continue
        valid_unsats = {float(value) for value in pd.to_numeric(valid_lines["不饱和度"], errors="coerce").dropna()}
        candidate_records: list[dict[str, object]] = []
        for unsat in sorted(pd.to_numeric(group_points["曲线不饱和度"], errors="coerce").dropna().unique()):
            unsat = float(unsat)
            if unsat in valid_unsats:
                continue
            lower_line, higher_line = find_bracket_lines(valid_lines, unsat)
            if lower_line is None or higher_line is None:
                continue
            points = group_points[group_points["曲线不饱和度"].astype(float).eq(unsat)].copy()
            rescue_line = build_rescue_line(str(plot_group), unsat, points, lower_line, higher_line, config)
            if rescue_line is not None:
                candidate_records.append(rescue_line)
        if not candidate_records:
            continue

        existing = [row for _, row in valid_lines.iterrows()]
        best_combo: tuple[dict[str, object], ...] | None = None
        best_score: tuple[int, float, int] = (-1, -1.0, -1)
        for size in range(len(candidate_records), 0, -1):
            for combo in itertools.combinations(candidate_records, size):
                combo_series = [pd.Series(item) for item in combo]
                if not subset_is_ordered(existing + combo_series, config):
                    continue
                rescue_score = subset_score(combo_series)
                point_count = sum(int(item.get("救回点数") or 0) for item in combo)
                current = (len(combo), rescue_score, point_count)
                if current > best_score:
                    best_combo = combo
                    best_score = current
            if best_combo is not None:
                break
        selected_records.extend(list(best_combo or []))

    for record in selected_records:
        mask = lines["细类"].eq(record["细类"]) & lines["不饱和度"].astype(float).eq(float(record["不饱和度"]))
        if not mask.any():
            lines = pd.concat([lines, pd.DataFrame([record])], ignore_index=True)
            continue
        row_index = lines.index[mask][0]
        for key, value in record.items():
            if key not in lines.columns:
                lines[key] = None
            lines.at[row_index, key] = value
    return lines


def apply_rescue_point_guard(
    rows: pd.DataFrame,
    lines: pd.DataFrame,
    inlier_column: str,
    note_column: str,
    config: RTIUPConfig = RTIUPConfig(),
) -> pd.DataFrame:
    """Re-check rescued points against their adjacent valid IUP lines."""

    if "二次捞点" not in rows.columns:
        return rows
    rows = rows.copy()
    valid_lines = lines[lines["是否有效曲线"].eq(True)].copy()
    rescue_mask = rows["二次捞点"].map(bool_value) & rows[inlier_column].eq(True)
    for row_idx, point in rows[rescue_mask].iterrows():
        plot_group = point.get("细类")
        unsat = finite(point.get("曲线不饱和度"))
        if unsat is None:
            rows.at[row_idx, inlier_column] = False
            continue
        line_match = valid_lines[
            valid_lines["细类"].eq(plot_group)
            & valid_lines["不饱和度"].astype(float).eq(float(unsat))
            & valid_lines["二次捞点"].map(bool_value)
        ]
        if line_match.empty:
            continue
        line = line_match.iloc[0]
        lower_unsat = finite(line.get("救回下界不饱和度"))
        higher_unsat = finite(line.get("救回上界不饱和度"))
        lower_line = valid_lines[
            valid_lines["细类"].eq(plot_group)
            & valid_lines["不饱和度"].astype(float).eq(float(lower_unsat))
        ] if lower_unsat is not None else pd.DataFrame()
        higher_line = valid_lines[
            valid_lines["细类"].eq(plot_group)
            & valid_lines["不饱和度"].astype(float).eq(float(higher_unsat))
        ] if higher_unsat is not None else pd.DataFrame()
        if lower_line.empty or higher_line.empty:
            rows.at[row_idx, inlier_column] = False
            rows.at[row_idx, note_column] = "二次IUP捞点失败：缺少相邻有效曲线"
            continue
        ok, note, _ = point_between_bracket(point, lower_line.iloc[0], higher_line.iloc[0], config)
        if not ok:
            rows.at[row_idx, inlier_column] = False
            rows.at[row_idx, note_column] = "二次IUP捞点失败：" + note
        else:
            rows.at[row_idx, note_column] = "二次IUP捞点：" + note
    return rows


def pearson_stats(x: np.ndarray, y: np.ndarray) -> tuple[float | None, str]:
    """Return Pearson r and formatted p-value for plot legends."""

    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return None, ""
    try:
        corr = stats.pearsonr(x, y)
    except Exception:
        return None, ""
    r_value = finite(corr.statistic)
    p_value = finite(corr.pvalue)
    if r_value is None or p_value is None:
        return None, ""
    if p_value < 0.0001:
        return r_value, "<0.0001"
    return r_value, f"{p_value:.4f}"


def equation_text(line: pd.Series) -> str:
    """Return a compact ``y``/``x`` equation string."""

    params = line.get("参数")
    if isinstance(params, str):
        params = json.loads(params)
    if line.get("拟合类型") == "Linear" and isinstance(params, list) and len(params) == 2:
        return f"y = {params[0]:.6g} * x + {params[1]:.6g}"
    if line.get("拟合类型") == "Quadratic" and isinstance(params, list) and len(params) == 3:
        return f"y = {params[0]:.6g} * x^2 + {params[1]:.6g} * x + {params[2]:.6g}"
    return ""


def prepare_rank_table(
    data: pd.DataFrame,
    dataset_name: str = "",
    series_col: str = "series",
    rt_col: str | None = None,
) -> pd.DataFrame:
    """Prepare SphinGOlipID rank tables for :func:`fit_rt_iup`."""

    out = data.copy()
    rt_source = rt_col or ("rt_2d_min" if "rt_2d_min" in out.columns else "RT_2d_min")
    if dataset_name:
        out["数据集"] = dataset_name
    out["细类"] = out[series_col].apply(series_to_plot_group)
    out["曲线不饱和度"] = out[series_col].apply(parse_series_unsaturation)
    out["x碳数"] = pd.to_numeric(out["x"], errors="coerce")
    out["归一化保留时间"] = pd.to_numeric(out[rt_source], errors="coerce")
    if "rank" in out.columns:
        out["候选排名"] = pd.to_numeric(out.get("rank"), errors="coerce")
    if "Abund" in out.columns:
        out["丰度"] = pd.to_numeric(out.get("Abund"), errors="coerce")
    out = out.dropna(subset=["细类", "曲线不饱和度", "x碳数", "归一化保留时间"]).copy()
    return out[out["细类"].astype(str).str.len() > 0].reset_index(drop=True)


def fit_ecn(
    data: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    keep_short_series: bool = True,
) -> RTIUPResult:
    """Fit ECN curves first, using carbon-median representatives.

    Curve parameters are estimated from median RT values at distinct carbon
    numbers. The final ±RT-window decision is then applied to every original
    row. Every curve is judged independently: cross-unsaturation parallelism
    is recorded for audit only and never removes an ECN curve.
    """

    required = {"细类", "曲线不饱和度", "x碳数", "归一化保留时间"}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")
    source = data.dropna(subset=list(required)).copy()
    line_records = [
        fit_line(group, config)
        for _, group in source.groupby(["细类", "曲线不饱和度"], sort=False)
    ]
    lines = pd.DataFrame(line_records)
    if lines.empty:
        return RTIUPResult(rows=source.iloc[0:0], lines=lines, plot_rows=source.iloc[0:0], stats={})
    lines["是否有效曲线"] = lines["拟合成功"].eq(True)
    lines = apply_parallel_checks(lines, config, invalidate=False)

    line_key_columns = [
        "细类", "不饱和度", "拟合成功", "拟合类型", "参数", "R²", "点数",
        "不同碳数", "代表点数", "内点数", "x最小", "x最大", "是否有效曲线",
        "去除原因", "平行性是否通过", "平行性参考曲线", "斜率绝对差",
        "斜率相对差", "平行性说明",
    ]
    line_key = lines[[column for column in line_key_columns if column in lines.columns]].copy()
    line_key = line_key.rename(columns={"不饱和度": "曲线不饱和度"})
    rows = source.merge(line_key, on=["细类", "曲线不饱和度"], how="left")

    predictions: list[float | None] = []
    for _, row in rows.iterrows():
        predictions.append(
            finite(predict_values(row.get("拟合类型"), row.get("参数"), row.get("x碳数")))
        )
    rows["预测保留时间"] = predictions
    rows["保留时间偏差"] = rows["归一化保留时间"] - rows["预测保留时间"]
    rows["保留时间绝对偏差"] = rows["保留时间偏差"].abs()
    rows["同曲线0.2min内点"] = (
        rows["拟合成功"].eq(True)
        & rows["是否有效曲线"].eq(True)
        & rows["保留时间绝对偏差"].le(config.rt_window_min)
    )
    rows["点数不足保留"] = False
    if keep_short_series:
        short_mask = rows["拟合成功"].eq(False) & rows["不同碳数"].between(1, 2, inclusive="both")
        for _, indices in rows[short_mask].groupby(["细类", "曲线不饱和度"], sort=False).groups.items():
            points = _median_points_by_x(rows.loc[indices, ["x碳数", "归一化保留时间"]])
            keep = not _has_clear_negative_short_slope(points, config)
            rows.loc[indices, "点数不足保留"] = keep

    rows["最终保留"] = rows["同曲线0.2min内点"] | rows["点数不足保留"]
    rows["是否作图"] = rows["最终保留"]
    rows["保留原因"] = np.where(
        rows["点数不足保留"],
        "点数不足但ECN未见明确下降，保留",
        np.where(rows["同曲线0.2min内点"], "ECN拟合通过且在0.2min内", ""),
    )
    final_rows = rows[rows["最终保留"]].copy()

    line_stats: list[dict[str, object]] = []
    for _, line in lines[lines["是否有效曲线"].eq(True)].sort_values(["细类", "不饱和度"]).iterrows():
        points = final_rows[
            final_rows["细类"].eq(line["细类"])
            & final_rows["曲线不饱和度"].astype(float).eq(float(line["不饱和度"]))
            & final_rows["同曲线0.2min内点"]
        ]
        x = points["x碳数"].astype(float).to_numpy()
        y = points["归一化保留时间"].astype(float).to_numpy()
        if len(points) < 3 or len(np.unique(x)) < 3:
            continue
        r_value, p_value = pearson_stats(x, y)
        record = line.to_dict()
        record["最终IUP点数"] = int(len(points))
        record["最终IUP不同碳数"] = int(len(np.unique(x)))
        record["Pearson r"] = r_value
        record["Pearson p"] = p_value
        record["拟合方程"] = equation_text(line)
        line_stats.append(record)
    valid_lines = pd.DataFrame(line_stats)
    stats = {
        "输入行数": int(len(source)),
        "最终保留行数": int(len(final_rows)),
        "拟合内点保留行数": int(final_rows["同曲线0.2min内点"].sum()),
        "点数不足保留行数": int(final_rows["点数不足保留"].sum()),
        "作图点数": int(len(final_rows)),
        "有效拟合曲线数": int(len(valid_lines)),
        "平行性未通过曲线数": int(lines["平行性是否通过"].eq(False).sum()),
    }
    return RTIUPResult(rows=final_rows, lines=valid_lines, plot_rows=final_rows, stats=stats)


def fit_rt_iup(
    data: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    keep_short_series: bool = True,
) -> RTIUPResult:
    """Fit IUP independently from raw points, then coordinate nearby curves.

    This is intentionally not derived from :func:`fit_ecn`. Initial curves are
    fitted from the full source data, missing/order-conflicting unsaturations
    are re-searched inside the locked RT window, and rescue candidates prefer
    slopes that are closer to their neighbouring IUP curves. Parallelism is a
    search preference and audit field, not an immediate reason to discard an
    entire unsaturation.
    """

    required = {"细类", "曲线不饱和度", "x碳数", "归一化保留时间"}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")

    data = data.dropna(subset=list(required)).copy()
    line_records = [fit_line(group, config) for _, group in data.groupby(["细类", "曲线不饱和度"], sort=False)]
    lines = pd.DataFrame(line_records)
    optimized_groups: list[pd.DataFrame] = []
    for plot_group, group_lines in lines.groupby("细类", sort=False):
        group_source = data[data["细类"].eq(plot_group)].copy()
        initial_choice = choose_iup_lines(group_lines, config)
        audit_choice = apply_parallel_checks(initial_choice, config, invalidate=False)
        successful_but_removed = (
            audit_choice["拟合成功"].eq(True)
            & audit_choice["是否有效曲线"].ne(True)
        ).any()
        parallel_issue = audit_choice["平行性是否通过"].eq(False).any()
        failed_but_searchable = (
            audit_choice["拟合成功"].ne(True)
            & pd.to_numeric(audit_choice.get("不同碳数"), errors="coerce").ge(3)
        ).any()
        needs_multistart = bool(
            config.enable_iup_rescue
            and (successful_but_removed or parallel_issue or failed_but_searchable)
        )
        optimized_groups.append(
            optimize_iup_lines_for_group(
                group_source,
                group_lines,
                config,
                multistart=needs_multistart,
            )
        )
    lines = pd.concat(optimized_groups, ignore_index=True) if optimized_groups else ensure_rescue_columns(lines)
    lines = apply_parallel_checks(lines, config, invalidate=False)

    line_key_columns = [
        "细类",
        "不饱和度",
        "拟合成功",
        "拟合类型",
        "参数",
        "R²",
        "点数",
        "内点数",
        "x最小",
        "x最大",
        "是否有效曲线",
        "去除原因",
        "二次捞点",
        "救回点数",
        "救回不同碳数",
        "救回下界不饱和度",
        "救回上界不饱和度",
        "救回说明",
        "平行性是否通过",
        "平行性参考曲线",
        "斜率绝对差",
        "斜率相对差",
        "平行性说明",
        "IUP候选来源",
        "IUP重新寻优",
        "IUP原始参数",
        "IUP_RT漂移校正(min)",
    ]
    line_key = lines[[col for col in line_key_columns if col in lines.columns]].copy()
    line_key = line_key.rename(columns={"不饱和度": "曲线不饱和度"})
    rows = data.merge(line_key, on=["细类", "曲线不饱和度"], how="left")
    rows["IUP_RT漂移校正(min)"] = pd.to_numeric(
        rows.get("IUP_RT漂移校正(min)", 0.0), errors="coerce"
    ).fillna(0.0)
    rows["原始归一化保留时间"] = rows["归一化保留时间"]
    rows["归一化保留时间"] = (
        pd.to_numeric(rows["归一化保留时间"], errors="coerce")
        + rows["IUP_RT漂移校正(min)"]
    )

    pred_values: list[float | None] = []
    for _, row in rows.iterrows():
        pred = predict_values(row.get("拟合类型"), row.get("参数"), row.get("x碳数"))
        pred_values.append(finite(pred))
    rows["预测保留时间"] = pred_values
    rows["保留时间偏差"] = rows["归一化保留时间"] - rows["预测保留时间"]
    rows["保留时间绝对偏差"] = rows["保留时间偏差"].abs()
    rows["同曲线0.2min内点"] = (
        rows["拟合成功"].eq(True)
        & rows["是否有效曲线"].eq(True)
        & rows["保留时间绝对偏差"].le(config.rt_window_min)
    )
    rows["点数不足保留"] = False
    if keep_short_series:
        rows["点数不足保留"] = rows["拟合成功"].eq(False) & rows["点数"].between(1, 2, inclusive="both")

    rows["IUP说明"] = np.where(rows["同曲线0.2min内点"], "有效拟合曲线，0.2min内点", "")
    rows = apply_rescue_point_guard(rows, lines, "同曲线0.2min内点", "IUP说明", config)
    rows["是否作图"] = rows["同曲线0.2min内点"]
    normal_mask = rows["是否作图"].eq(True) & ~rows.get("二次捞点", pd.Series(False, index=rows.index)).map(bool_value)
    rows.loc[normal_mask, "IUP说明"] = "有效拟合曲线，0.2min内点"
    rows = apply_point_level_iup_guard(rows, lines, config)

    if keep_short_series:
        if "点数不足过滤原因" not in rows.columns:
            rows["点数不足过滤原因"] = ""
        short_rows = rows[rows["点数不足保留"]]
        for (plot_group, unsat), idx in short_rows.groupby(["细类", "曲线不饱和度"], sort=False).groups.items():
            valid = lines[lines["细类"].eq(plot_group) & lines["是否有效曲线"].eq(True)].copy()
            group = rows.loc[idx].copy()
            unique_points = group[["x碳数", "归一化保留时间"]].dropna().drop_duplicates()
            if len(unique_points) == 2:
                ok, note = two_point_series_iup_status(group, valid, config)
                rows.loc[idx, "是否作图"] = ok
                rows.loc[idx, "IUP说明"] = note
                if not ok:
                    rows.loc[idx, "点数不足保留"] = False
                    rows.loc[idx, "点数不足过滤原因"] = note
                continue
            for row_idx in idx:
                ok, note = point_iup_status(rows.loc[row_idx], valid, config)
                rows.at[row_idx, "是否作图"] = ok
                rows.at[row_idx, "IUP说明"] = note
        rows = apply_short_series_guard(rows, lines, config)

    rows["最终保留"] = rows["同曲线0.2min内点"] | rows["点数不足保留"]
    rows["保留原因"] = np.where(
        rows["点数不足保留"],
        "点数不足，保留",
        np.where(
            rows["同曲线0.2min内点"] & rows.get("二次捞点", pd.Series(False, index=rows.index)).map(bool_value),
            "二次IUP捞点恢复",
            np.where(rows["同曲线0.2min内点"], "拟合成功且符合IUP，0.2min内点", ""),
        ),
    )
    final_rows = rows[rows["最终保留"]].copy()

    line_stats: list[dict[str, object]] = []
    for _, line in lines[lines["是否有效曲线"].eq(True)].sort_values(["细类", "不饱和度"]).iterrows():
        pts = rows[
            rows["细类"].eq(line["细类"])
            & rows["曲线不饱和度"].astype(float).eq(float(line["不饱和度"]))
            & rows["同曲线0.2min内点"]
        ].copy()
        x = pts["x碳数"].astype(float).to_numpy()
        y = pts["归一化保留时间"].astype(float).to_numpy()
        final_point_count = int(len(pts))
        final_distinct_x = int(pd.to_numeric(pts["x碳数"], errors="coerce").dropna().nunique())
        if final_point_count < 3 or final_distinct_x < 2:
            continue
        r_value, p_value = pearson_stats(x, y)
        record = line.to_dict()
        record["最终IUP点数"] = final_point_count
        record["最终IUP不同碳数"] = final_distinct_x
        record["Pearson r"] = r_value
        record["Pearson p"] = p_value
        record["拟合方程"] = equation_text(line)
        line_stats.append(record)
    valid_lines = pd.DataFrame(line_stats)
    plot_rows = final_rows[final_rows["是否作图"].eq(True)].copy()
    stats_payload = {
        "输入行数": int(len(data)),
        "最终保留行数": int(len(final_rows)),
        "拟合内点保留行数": int(final_rows["同曲线0.2min内点"].sum()),
        "点级IUP过滤行数": int(rows.get("点级IUP过滤原因", pd.Series("", index=rows.index)).fillna("").astype(str).ne("").sum()),
        "点数不足保留行数": int(final_rows["点数不足保留"].sum()),
        "点数不足过滤行数": int(rows.get("点数不足过滤原因", pd.Series("", index=rows.index)).astype(str).ne("").sum()),
        "拟合点IUP过滤曲线数": int(
            lines.get("去除原因", pd.Series("", index=lines.index)).astype(str).str.contains("拟合曲线IUP支持点不足", regex=False).sum()
        ),
        "作图点数": int(len(plot_rows)),
        "有效拟合曲线数": int(len(valid_lines)),
        "去除曲线数": int(lines["拟合成功"].eq(True).sum() - lines["是否有效曲线"].eq(True).sum()),
        "二次捞点曲线数": int(lines.get("二次捞点", pd.Series(False, index=lines.index)).map(bool_value).sum()),
        "二次捞点保留行数": int(final_rows.get("二次捞点", pd.Series(False, index=final_rows.index)).map(bool_value).sum()),
        "IUP重新寻优曲线数": int(lines.get("IUP重新寻优", pd.Series(False, index=lines.index)).map(bool_value).sum()),
        "IUP重新寻优保留行数": int(final_rows.get("IUP重新寻优", pd.Series(False, index=final_rows.index)).map(bool_value).sum()),
    }
    return RTIUPResult(rows=final_rows, lines=valid_lines, plot_rows=plot_rows, stats=stats_payload)


def filter_iup_from_ecn(
    ecn: RTIUPResult,
    config: RTIUPConfig = RTIUPConfig(),
) -> RTIUPResult:
    """Apply IUP ordering only to an already completed ECN result.

    This function cannot fit or rescue a curve.  Both curve keys and row IDs
    are asserted to remain strict subsets of their ECN inputs.
    """

    ecn_lines = ecn.lines.copy()
    ecn_rows = ecn.rows.copy()
    if ecn_lines.empty or ecn_rows.empty:
        return RTIUPResult(
            rows=ecn_rows.iloc[0:0].copy(),
            lines=ecn_lines.iloc[0:0].copy(),
            plot_rows=ecn_rows.iloc[0:0].copy(),
            stats={"输入行数": int(len(ecn_rows)), "最终保留行数": 0, "有效拟合曲线数": 0},
        )

    selected_groups: list[pd.DataFrame] = []
    for _, group in ecn_lines.groupby("细类", sort=False):
        candidate = group.copy()
        candidate["拟合成功"] = True
        candidate["是否有效曲线"] = True
        selected_groups.append(choose_iup_lines(candidate, config))
    selected = pd.concat(selected_groups, ignore_index=True) if selected_groups else ecn_lines.iloc[0:0].copy()
    valid_lines = selected[selected["是否有效曲线"].eq(True)].copy()
    allowed_line_keys = set(
        zip(valid_lines["细类"].astype(str), valid_lines["不饱和度"].astype(float))
    )

    rows = ecn_rows.copy()
    row_keys = list(zip(rows["细类"].astype(str), rows["曲线不饱和度"].astype(float)))
    fitted = rows.get("拟合成功", pd.Series(False, index=rows.index)).eq(True)
    short = rows.get("点数不足保留", pd.Series(False, index=rows.index)).eq(True)
    keep = pd.Series(
        [((not is_fitted) or key in allowed_line_keys) and (is_fitted or is_short)
         for key, is_fitted, is_short in zip(row_keys, fitted, short)],
        index=rows.index,
        dtype=bool,
    )

    # Short ECN series may remain only when their points obey the already-kept
    # neighbouring ECN curves.  No curve is fitted or added here.
    for (plot_group, unsat), indices in rows[short & keep].groupby(
        ["细类", "曲线不饱和度"], sort=False
    ).groups.items():
        reference_lines = valid_lines[valid_lines["细类"].eq(plot_group)].copy()
        group = rows.loc[indices].copy()
        unique_points = _median_points_by_x(group[["x碳数", "归一化保留时间"]])
        if len(unique_points) == 2:
            ok, _ = two_point_series_iup_status(group, reference_lines, config)
            if not ok:
                keep.loc[indices] = False
        else:
            for row_index in indices:
                ok, _ = point_iup_status(rows.loc[row_index], reference_lines, config)
                if not ok:
                    keep.at[row_index] = False

    # Direct same-carbon order check, with the locked two-feature drift and IUP
    # allowances.  It can only remove rows from the ECN subset.
    tolerance = relative_rt_tolerance_min(config)
    candidate_rows = rows[keep]
    for _, group in candidate_rows.groupby(["细类", "x碳数"], sort=False):
        medians = group.groupby("曲线不饱和度", sort=True)["归一化保留时间"].median().sort_index()
        last_kept_rt: float | None = None
        for unsat, rt_value in medians.items():
            current_rt = float(rt_value)
            if last_kept_rt is not None and current_rt > last_kept_rt + tolerance:
                bad = group["曲线不饱和度"].astype(float).eq(float(unsat))
                keep.loc[group.index[bad]] = False
                continue
            last_kept_rt = current_rt if last_kept_rt is None else min(last_kept_rt, current_rt)

    final_rows = rows[keep].copy()
    final_rows["最终保留"] = True
    final_rows["是否作图"] = True
    final_rows["保留原因"] = np.where(
        final_rows.get("点数不足保留", pd.Series(False, index=final_rows.index)).eq(True),
        "ECN点数不足保留且通过IUP顺序过滤",
        "ECN拟合通过并通过IUP顺序过滤",
    )
    plot_rows = final_rows.copy()

    remaining_curve_keys = set(
        zip(
            final_rows.loc[final_rows.get("拟合成功", pd.Series(False, index=final_rows.index)).eq(True), "细类"].astype(str),
            final_rows.loc[final_rows.get("拟合成功", pd.Series(False, index=final_rows.index)).eq(True), "曲线不饱和度"].astype(float),
        )
    )
    valid_lines = valid_lines[
        [key in remaining_curve_keys for key in zip(valid_lines["细类"].astype(str), valid_lines["不饱和度"].astype(float))]
    ].copy()

    ecn_keys = set(zip(ecn_lines["细类"].astype(str), ecn_lines["不饱和度"].astype(float)))
    iup_keys = set(zip(valid_lines["细类"].astype(str), valid_lines["不饱和度"].astype(float)))
    assert iup_keys <= ecn_keys
    if "合并后行ID" in ecn_rows.columns and "合并后行ID" in final_rows.columns:
        assert set(final_rows["合并后行ID"].astype(str)) <= set(ecn_rows["合并后行ID"].astype(str))
    else:
        assert set(final_rows.index) <= set(ecn_rows.index)

    stats = {
        "输入行数": int(len(ecn_rows)),
        "最终保留行数": int(len(final_rows)),
        "拟合内点保留行数": int(final_rows.get("同曲线0.2min内点", pd.Series(False, index=final_rows.index)).sum()),
        "点数不足保留行数": int(final_rows.get("点数不足保留", pd.Series(False, index=final_rows.index)).sum()),
        "作图点数": int(len(plot_rows)),
        "有效拟合曲线数": int(len(valid_lines)),
        "IUP从ECN删除行数": int(len(ecn_rows) - len(final_rows)),
        "IUP从ECN删除曲线数": int(len(ecn_lines) - len(valid_lines)),
    }
    return RTIUPResult(rows=final_rows, lines=valid_lines, plot_rows=plot_rows, stats=stats)


def _design_matrix(x: np.ndarray, fit_type: object) -> np.ndarray:
    if fit_type == "Quadratic":
        return np.column_stack([x * x, x, np.ones_like(x)])
    return np.column_stack([x, np.ones_like(x)])


def regression_ci(
    line: pd.Series,
    x: np.ndarray,
    y: np.ndarray,
    x_grid: np.ndarray,
    config: RTIUPConfig = RTIUPConfig(),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return fitted y-grid and 95% CI capped by the RT residual window."""

    y_grid = fit_prediction(line, x_grid)
    if y_grid is None:
        y_grid = np.full_like(x_grid, np.nan, dtype=float)
    y_grid = np.asarray(y_grid, dtype=float)
    y_fit = fit_prediction(line, x)
    if y_fit is None:
        return y_grid, np.full_like(y_grid, np.nan), np.full_like(y_grid, np.nan)
    y_fit = np.asarray(y_fit, dtype=float)

    x_design = _design_matrix(x, line.get("拟合类型"))
    x_grid_design = _design_matrix(x_grid, line.get("拟合类型"))
    n_obs, n_params = x_design.shape
    dof = n_obs - n_params
    if dof <= 0:
        return y_grid, np.full_like(y_grid, np.nan), np.full_like(y_grid, np.nan)
    residual = y - y_fit
    s2 = float(np.sum(residual * residual) / dof)
    if not math.isfinite(s2):
        return y_grid, np.full_like(y_grid, np.nan), np.full_like(y_grid, np.nan)
    cov = np.linalg.pinv(x_design.T @ x_design)
    leverage = np.sum((x_grid_design @ cov) * x_grid_design, axis=1)
    se_mean = np.sqrt(np.maximum(0, s2 * leverage))
    t_value = stats.t.ppf(0.975, dof)
    half_width = np.minimum(t_value * se_mean, config.rt_window_min)
    return y_grid, y_grid - half_width, y_grid + half_width


def set_plot_style(config: RTIUPConfig = RTIUPConfig()) -> None:
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["axes.grid"] = False
    plt.rcParams["axes.labelsize"] = config.axis_label_size
    plt.rcParams["xtick.labelsize"] = config.tick_label_size
    plt.rcParams["ytick.labelsize"] = config.tick_label_size


def legend_label(unsat: object, group: pd.DataFrame, line: pd.Series | None = None) -> str:
    """Return a compact R-squared legend label for one unsaturation series."""

    base = f"{int(float(unsat))}"
    if line is None:
        return base

    r_squared = finite(line.get("R²"))
    if r_squared is None:
        return base
    return f"{base}  R²={r_squared:.4f}"


def make_rt_iup_plot(
    plot_group: str,
    plot_rows: pd.DataFrame,
    valid_lines: pd.DataFrame,
    out_path: str | Path,
    config: RTIUPConfig = RTIUPConfig(),
) -> None:
    """Draw one RT/IUP plot with points, fit lines, capped 95% CIs and legends."""

    set_plot_style(config)
    fig, ax = plt.subplots(figsize=config.figure_size)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.tick_params(axis="both", which="both", direction="out", colors="black", length=7, width=2.0)

    clean_points = plot_rows.dropna(subset=["x碳数", "归一化保留时间", "曲线不饱和度"]).copy()
    legend_handles: list[Line2D] = []
    legend_labels: list[str] = []
    for unsat, group in sorted(clean_points.groupby("曲线不饱和度"), key=lambda item: float(item[0])):
        color = color_for_unsaturation(unsat, config)
        ax.scatter(
            group["x碳数"],
            group["归一化保留时间"],
            color=color,
            alpha=0.7,
            s=config.point_size,
            edgecolors="none",
            zorder=3,
        )
        two_point_line = group[["x碳数", "归一化保留时间"]].drop_duplicates().sort_values(["x碳数", "归一化保留时间"])
        if len(two_point_line) == 2:
            ax.plot(
                two_point_line["x碳数"],
                two_point_line["归一化保留时间"],
                color=color,
                linestyle="-",
                linewidth=config.two_point_line_width,
                alpha=0.82,
                zorder=2,
            )
        matching_lines = valid_lines[
            pd.to_numeric(valid_lines.get("不饱和度"), errors="coerce").eq(float(unsat))
        ] if not valid_lines.empty and "不饱和度" in valid_lines.columns else valid_lines.iloc[0:0]
        legend_line = matching_lines.iloc[0] if not matching_lines.empty else None
        legend_handles.append(
            Line2D([0], [0], marker="o", linestyle="None", markersize=9.5, markerfacecolor=color, markeredgecolor="none", alpha=0.7)
        )
        legend_labels.append(legend_label(unsat, group, legend_line))

    for _, line in valid_lines.iterrows():
        fit_unsat = finite(line.get("不饱和度"))
        if fit_unsat is None:
            continue
        group = clean_points[clean_points["曲线不饱和度"].astype(float).eq(fit_unsat)]
        if len(group) < 3:
            continue
        x = group["x碳数"].astype(float).to_numpy()
        y = group["归一化保留时间"].astype(float).to_numpy()
        if len(np.unique(x)) < 2:
            continue
        x_grid = np.linspace(float(np.min(x)), float(np.max(x)), 140)
        y_grid, y_low, y_high = regression_ci(line, x, y, x_grid, config)
        color = color_for_unsaturation(fit_unsat, config)
        if np.all(np.isfinite(y_low)) and np.all(np.isfinite(y_high)):
            ax.fill_between(x_grid, y_low, y_high, color=color, alpha=0.22, linewidth=0, zorder=1)
        ax.plot(x_grid, y_grid, color=color, linestyle="-", linewidth=config.fit_line_width, zorder=2)

    ax.set_xlabel("Carbon Number", fontname="Arial", fontweight="bold", labelpad=10)
    ax.set_ylabel("Normalized Retention Time(min)", fontname="Arial", fontweight="bold", labelpad=10)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontname("Arial")
    if legend_handles:
        ax.legend(
            legend_handles,
            legend_labels,
            loc="upper left",
            frameon=False,
            prop={"family": "Arial", "size": config.legend_font_size},
            handlelength=1.0,
            handletextpad=0.6,
            borderaxespad=0.8,
        )
    ax.margins(x=0.06, y=0.10)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def write_rt_iup_plots(
    plot_rows: pd.DataFrame,
    valid_lines: pd.DataFrame,
    out_dir: str | Path,
    config: RTIUPConfig = RTIUPConfig(),
) -> int:
    """Write one ``*_all.png`` plot per lipid subclass and return the count."""

    out_dir = Path(out_dir)
    plot_count = 0
    for plot_group, group_plot_rows in plot_rows.groupby("细类", sort=False):
        group_lines = valid_lines[valid_lines["细类"].eq(plot_group)] if not valid_lines.empty else valid_lines
        if group_plot_rows.empty:
            continue
        make_rt_iup_plot(plot_group, group_plot_rows, group_lines, out_dir / f"{safe_name(plot_group)}_all.png", config)
        plot_count += 1
    return plot_count


__all__ = [
    "RTIUPConfig",
    "RTIUPResult",
    "apply_fitted_point_iup_guard",
    "apply_parallel_checks",
    "apply_rescue_point_guard",
    "bool_value",
    "choose_iup_lines",
    "color_for_unsaturation",
    "equation_text",
    "derivative_values",
    "finite",
    "fit_line",
    "fit_ecn",
    "fit_rt_iup",
    "filter_iup_from_ecn",
    "iterative_refit",
    "legend_label",
    "make_rt_iup_plot",
    "normalize_json_value",
    "parse_series_unsaturation",
    "pair_parallel_status",
    "pearson_stats",
    "apply_point_level_iup_guard",
    "point_level_iup_status",
    "point_iup_status",
    "prepare_rank_table",
    "predict_values",
    "regression_ci",
    "rescue_missing_iup_lines",
    "safe_name",
    "series_to_plot_group",
    "set_plot_style",
    "write_rt_iup_plots",
]
