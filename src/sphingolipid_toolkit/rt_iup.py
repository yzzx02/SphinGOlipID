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
from scipy import stats
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
import warnings


warnings.filterwarnings("ignore", category=UndefinedMetricWarning)


DEFAULT_COLORS = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]


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


def _fit_candidate_model(points: pd.DataFrame, fit_type: str, config: RTIUPConfig) -> dict[str, object] | None:
    min_samples = 2 if fit_type == "Linear" else 3
    if len(points) < min_samples or points["x碳数"].nunique() < 2:
        return None

    x = points["x碳数"].astype(float).to_numpy().reshape(-1, 1)
    y = points["归一化保留时间"].astype(float).to_numpy()
    estimator = LinearRegression() if fit_type == "Linear" else make_pipeline(
        PolynomialFeatures(2, include_bias=False),
        LinearRegression(),
    )
    model = RANSACRegressor(
        estimator=estimator,
        min_samples=min_samples,
        residual_threshold=config.rt_window_min,
        random_state=config.random_state,
        max_trials=300,
    )
    try:
        model.fit(x, y)
    except Exception:
        return None

    pred = model.predict(x)
    inliers = np.abs(y - pred) <= config.rt_window_min
    if int(inliers.sum()) < 3:
        return None

    r2 = float(r2_score(y[inliers], pred[inliers])) if int(inliers.sum()) >= 2 else np.nan
    params = _model_params(model, fit_type)
    x_min = float(np.min(x[inliers]))
    x_max = float(np.max(x[inliers]))
    inlier_points = points.loc[inliers].copy()
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
        "点数": int(len(points)),
        "内点数": int(inliers.sum()),
        "x最小": x_min,
        "x最大": x_max,
    }


def fit_line(group: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> dict[str, object]:
    """Fit one unsaturation curve using linear/quadratic RANSAC candidates."""

    plot_group = str(group["细类"].iloc[0])
    line_unsat = float(group["曲线不饱和度"].iloc[0])
    points = (
        group[["x碳数", "归一化保留时间"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["x碳数", "归一化保留时间"])
        .reset_index(drop=True)
    )
    base = {
        "细类": plot_group,
        "不饱和度": line_unsat,
        "拟合成功": False,
        "拟合类型": None,
        "参数": None,
        "R²": None,
        "点数": int(len(points)),
        "内点数": 0,
        "x最小": None,
        "x最大": None,
        "失败原因": "",
        "是否有效曲线": False,
        "去除原因": "",
    }
    if len(points) < 3:
        base["失败原因"] = "少于3个点"
        return base
    if points["x碳数"].nunique() < 2:
        base["失败原因"] = "有效x少于2个"
        return base

    candidates = [_fit_candidate_model(points, "Linear", config), _fit_candidate_model(points, "Quadratic", config)]
    candidates = [item for item in candidates if item is not None]
    if not candidates:
        base["失败原因"] = "未达到R²或单调ECN要求"
        return base

    def score(item: dict[str, object]) -> tuple[float, int, int]:
        complexity_bonus = 0 if item["拟合类型"] == "Linear" else -1
        return (float(item["R²"]), int(item["内点数"]), complexity_bonus)

    best = sorted(candidates, key=score, reverse=True)[0]
    base.update(best)
    base["拟合成功"] = True
    return base


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
    bad_fraction = float(np.mean(diff < -config.order_tolerance_min))
    return bool(
        bad_fraction >= 0.75
        and float(np.median(diff)) < -config.order_tolerance_min
        and float(np.mean(diff)) < -config.order_tolerance_min
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


def point_iup_status(point: pd.Series, valid_lines: pd.DataFrame, config: RTIUPConfig = RTIUPConfig()) -> tuple[bool, str]:
    """Return whether a 1-2 point series can be shown by IUP order alone."""

    x = finite(point.get("x碳数"))
    y = finite(point.get("归一化保留时间"))
    unsat = finite(point.get("曲线不饱和度"))
    if x is None or y is None or unsat is None:
        return False, "缺少x/保留时间/不饱和度"

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
            if pred + config.iup_tolerance_min < y:
                violations.append(f"低不饱和度{other_unsat:g}曲线未在上方")
            else:
                supports.append(f"低不饱和度{other_unsat:g}曲线在上方")
        else:
            if y + config.iup_tolerance_min < pred:
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
    guarded["点数不足过滤原因"] = ""
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
    if float(upper_y) + config.order_tolerance_min < float(lower_y):
        return False, "相邻曲线在该x处交叉", None
    if y > float(upper_y) + config.iup_tolerance_min:
        return False, "高于低不饱和度相邻曲线", None
    if y < float(lower_y) - config.iup_tolerance_min:
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


def _is_judgeable_iup_point(note: str) -> bool:
    skipped = (
        "缺少x或保留时间",
        "相邻曲线缺少x范围",
        "x不在相邻有效曲线重叠范围",
        "相邻曲线无法预测",
        "相邻曲线在该x处交叉",
    )
    return not any(text in note for text in skipped)


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
    if np.mean((upper_y[ok] + config.order_tolerance_min) < lower_y[ok]) > 0.10:
        return False
    inside = (y[ok] <= upper_y[ok] + config.iup_tolerance_min) & (y[ok] >= lower_y[ok] - config.iup_tolerance_min)
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
    x = points["x碳数"].astype(float).to_numpy()
    y = points["归一化保留时间"].astype(float).to_numpy()
    if len(points) < config.rescue_min_distinct_x or len(np.unique(x)) < config.rescue_min_distinct_x:
        return None
    degree = 1 if fit_type == "Linear" else 2
    if len(np.unique(x)) < degree + 1:
        return None
    try:
        coeff = np.polyfit(x, y, degree)
    except Exception:
        return None
    params = [float(coeff[0]), float(coeff[1])] if fit_type == "Linear" else [float(coeff[0]), float(coeff[1]), float(coeff[2])]
    pred = predict_values(fit_type, params, x)
    if pred is None:
        return None
    pred = np.asarray(pred, dtype=float)
    residual = np.abs(y - pred)
    if not np.all(np.isfinite(residual)) or float(np.max(residual)) > config.rt_window_min:
        return None
    r2 = float(r2_score(y, pred)) if len(points) >= 2 else np.nan
    if not math.isfinite(r2) or r2 < config.r2_threshold:
        return None
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    if not _is_positive_monotonic(fit_type, params, x_min, x_max, config):
        return None
    return {
        "拟合类型": fit_type,
        "参数": params,
        "R²": r2,
        "点数": int(len(points)),
        "内点数": int(len(points)),
        "x最小": x_min,
        "x最大": x_max,
    }


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

    def score(item: dict[str, object]) -> tuple[float, int, int]:
        complexity_bonus = 0 if item["拟合类型"] == "Linear" else -1
        return (float(item["R²"]), int(item["内点数"]), complexity_bonus)

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


def fit_rt_iup(
    data: pd.DataFrame,
    config: RTIUPConfig = RTIUPConfig(),
    keep_short_series: bool = True,
) -> RTIUPResult:
    """Run RT filtering using RANSAC fits, IUP ordering, and strict rescue."""

    required = {"细类", "曲线不饱和度", "x碳数", "归一化保留时间"}
    missing = required - set(data.columns)
    if missing:
        raise KeyError(f"Missing required column(s): {sorted(missing)}")

    data = data.dropna(subset=list(required)).copy()
    line_records = [fit_line(group, config) for _, group in data.groupby(["细类", "曲线不饱和度"], sort=False)]
    lines = pd.DataFrame(line_records)
    chosen_lines = [choose_iup_lines(group, config) for _, group in lines.groupby("细类", sort=False)]
    lines = pd.concat(chosen_lines, ignore_index=True) if chosen_lines else ensure_rescue_columns(lines)
    lines = rescue_missing_iup_lines(data, lines, config)

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
    ]
    line_key = lines[[col for col in line_key_columns if col in lines.columns]].copy()
    line_key = line_key.rename(columns={"不饱和度": "曲线不饱和度"})
    rows = data.merge(line_key, on=["细类", "曲线不饱和度"], how="left")

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
    rows, lines = apply_fitted_point_iup_guard(rows, lines, config)

    if keep_short_series:
        for plot_group, idx in rows[rows["点数不足保留"]].groupby("细类").groups.items():
            valid = lines[lines["细类"].eq(plot_group) & lines["是否有效曲线"].eq(True)].copy()
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
        r_value, p_value = pearson_stats(x, y)
        record = line.to_dict()
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
    }
    return RTIUPResult(rows=final_rows, lines=valid_lines, plot_rows=plot_rows, stats=stats_payload)


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


def legend_label(unsat: object, group: pd.DataFrame) -> str:
    base = f"{int(float(unsat))}"
    x = group["x碳数"].astype(float).to_numpy()
    y = group["归一化保留时间"].astype(float).to_numpy()
    r_value, p_value = pearson_stats(x, y)
    if r_value is None:
        return base
    return f"{base}  r={r_value:.4f}, p={p_value}"


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
        legend_handles.append(
            Line2D([0], [0], marker="o", linestyle="None", markersize=9.5, markerfacecolor=color, markeredgecolor="none", alpha=0.7)
        )
        legend_labels.append(legend_label(unsat, group))

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
    "apply_rescue_point_guard",
    "bool_value",
    "choose_iup_lines",
    "color_for_unsaturation",
    "equation_text",
    "finite",
    "fit_line",
    "fit_rt_iup",
    "legend_label",
    "make_rt_iup_plot",
    "normalize_json_value",
    "parse_series_unsaturation",
    "pearson_stats",
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
