"""Regenerate TOP3 ECN/IUP plots without overwriting archived figures.

The script writes every single-group plot plus six-panel summary figures for
three grouping schemes: exact chain, LCB type, and total carbon/unsaturation.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sphingolipid_toolkit.rt_iup import (  # noqa: E402
    RTIUPConfig,
    RTIUPResult,
    color_for_unsaturation,
    equation_text,
    fit_ecn as fit_ecn_core,
    fit_line,
    fit_rt_iup,
    pearson_stats,
    predict_values,
    relative_rt_tolerance_min,
    safe_name,
    write_rt_iup_plots,
)


C_CLASS = "细类"
C_NOTE = "注释"
C_PRECURSOR = "母离子"
C_X = "x碳数"
C_UNSAT = "曲线不饱和度"
C_RT = "归一化保留时间"
C_SCORE = "匹配度分数"

PLOT_STYLE = {
    "figure_size": (4.72, 4.72),
    "point_size": 26,
    "two_point_line_width": 1.2,
    "fit_line_width": 1.5,
    "legend_font_size": 8,
    "axis_label_size": 10,
    "tick_label_size": 9,
}

TOTAL_RE = re.compile(r"\(([mdt])\s*(\d+):(\d+)\)", re.IGNORECASE)
BARE_TOTAL_RE = re.compile(r"\b([mdt])\s*(\d+):(\d+)\b", re.IGNORECASE)
CHAIN_RE = re.compile(r"(?:(?P<base>[mdt])\s*)?(?P<c>\d{1,3}):(?P<u>\d{1,2})", re.IGNORECASE)
OH_RE = re.compile(r"\((?P<oh>(?:\d+)?OH)\)", re.IGNORECASE)


def text_value(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def parse_row(row: pd.Series) -> tuple[str, str, int | None, int | None, str]:
    note = text_value(row.get(C_NOTE))
    precursor = text_value(row.get(C_PRECURSOR))
    source = note or precursor
    head = re.sub(r"\s+", " ", source.split("(", 1)[0].strip()) if source else ""

    base = ""
    total_c: int | None = None
    total_u: int | None = None
    for pattern in (TOTAL_RE, BARE_TOTAL_RE):
        match = pattern.search(precursor)
        if match:
            base = match.group(1).lower()
            total_c = int(match.group(2))
            total_u = int(match.group(3))
            break
    if total_c is None or total_u is None:
        matches = list(CHAIN_RE.finditer(note or precursor))
        if matches:
            base = base or next((m.group("base").lower() for m in matches if m.group("base")), "")
            total_c = sum(int(m.group("c")) for m in matches)
            total_u = sum(int(m.group("u")) for m in matches)

    markers: list[str] = []
    for item in (precursor, note):
        for match in OH_RE.finditer(item):
            marker = match.group("oh").upper()
            if marker not in markers:
                markers.append(marker)
    precise = [marker for marker in markers if marker != "OH"]
    oh = "+".join(precise or markers)
    return head, base, total_c, total_u, oh


def prepare_groupings(raw: pd.DataFrame, score_threshold: float) -> dict[str, pd.DataFrame]:
    score = pd.to_numeric(raw.get(C_SCORE), errors="coerce")
    filtered = raw.loc[score.ge(score_threshold)].copy()
    for column in (C_X, C_UNSAT, C_RT):
        filtered[column] = pd.to_numeric(filtered[column], errors="coerce")

    detailed = filtered.dropna(subset=[C_CLASS, C_X, C_UNSAT, C_RT]).copy()
    parsed = pd.DataFrame(
        [parse_row(row) for _, row in filtered.iterrows()],
        columns=["脂类头", "长链碱类型", "总碳数", "总不饱和度", "OH标记"],
        index=filtered.index,
    )
    enriched = pd.concat([filtered, parsed], axis=1)
    valid = (
        enriched["脂类头"].astype(str).str.len().gt(0)
        & enriched["总碳数"].notna()
        & enriched["总不饱和度"].notna()
        & enriched[C_RT].notna()
    )
    enriched = enriched.loc[valid].copy()

    def with_group(frame: pd.DataFrame, group: pd.Series) -> pd.DataFrame:
        out = frame.copy()
        out[C_CLASS] = group
        out[C_X] = pd.to_numeric(out["总碳数"], errors="coerce")
        out[C_UNSAT] = pd.to_numeric(out["总不饱和度"], errors="coerce")
        return out.dropna(subset=[C_CLASS, C_X, C_UNSAT, C_RT]).reset_index(drop=True)

    oh_suffix = enriched["OH标记"].map(lambda value: f" ({value})" if text_value(value) else "")
    lcb_group = enriched["脂类头"].astype(str) + enriched["长链碱类型"].map(lambda value: f"({value})" if text_value(value) else "") + oh_suffix
    total_group = enriched["脂类头"].astype(str) + oh_suffix
    return {
        "详细链分类": detailed.reset_index(drop=True),
        "长链碱系列分类": with_group(enriched, lcb_group),
        "总碳数总不饱和度分类": with_group(enriched, total_group),
    }


def fit_ecn(data: pd.DataFrame, config: RTIUPConfig) -> RTIUPResult:
    """Fit each ECN curve independently, without cross-curve IUP ordering."""

    required = [C_CLASS, C_UNSAT, C_X, C_RT]
    source = data.dropna(subset=required).copy()
    records = [fit_line(group, config) for _, group in source.groupby([C_CLASS, C_UNSAT], sort=False)]
    lines = pd.DataFrame(records)
    if lines.empty:
        return RTIUPResult(rows=source.iloc[0:0], lines=lines, plot_rows=source.iloc[0:0], stats={})
    lines["是否有效曲线"] = lines["拟合成功"].eq(True)
    key = lines.rename(columns={"不饱和度": C_UNSAT})
    keep_cols = [C_CLASS, C_UNSAT, "拟合成功", "拟合类型", "参数", "R²", "点数", "内点数", "x最小", "x最大", "是否有效曲线"]
    rows = source.merge(key[keep_cols], on=[C_CLASS, C_UNSAT], how="left")

    predictions = []
    for _, row in rows.iterrows():
        predictions.append(predict_values(row.get("拟合类型"), row.get("参数"), row.get(C_X)))
    rows["预测保留时间"] = pd.to_numeric(pd.Series(predictions, index=rows.index), errors="coerce")
    rows["保留时间绝对偏差"] = (rows[C_RT] - rows["预测保留时间"]).abs()
    rows["同曲线0.2min内点"] = rows["拟合成功"].eq(True) & rows["保留时间绝对偏差"].le(config.rt_window_min)
    rows["点数不足保留"] = False

    short = rows[rows["拟合成功"].eq(False) & rows["点数"].between(1, 2, inclusive="both")]
    for _, indices in short.groupby([C_CLASS, C_UNSAT], sort=False).groups.items():
        pts = rows.loc[indices, [C_X, C_RT]].dropna().drop_duplicates().sort_values([C_X, C_RT])
        keep = True
        if len(pts) == 2 and pts[C_X].nunique() == 2:
            keep = bool(float(pts.iloc[1][C_RT]) + relative_rt_tolerance_min(config) >= float(pts.iloc[0][C_RT]))
        rows.loc[indices, "点数不足保留"] = keep

    rows["最终保留"] = rows["同曲线0.2min内点"] | rows["点数不足保留"]
    rows["是否作图"] = rows["最终保留"]
    rows["保留原因"] = np.where(rows["点数不足保留"], "点数不足但ECN未明显下降，保留", "拟合成功且在0.2min内")
    final_rows = rows[rows["最终保留"]].copy()

    line_stats: list[dict[str, object]] = []
    for _, line in lines[lines["是否有效曲线"]].sort_values([C_CLASS, "不饱和度"]).iterrows():
        pts = final_rows[
            final_rows[C_CLASS].eq(line[C_CLASS])
            & final_rows[C_UNSAT].astype(float).eq(float(line["不饱和度"]))
            & final_rows["同曲线0.2min内点"]
        ]
        x = pts[C_X].astype(float).to_numpy()
        y = pts[C_RT].astype(float).to_numpy()
        if len(pts) < 3 or len(np.unique(x)) < 2:
            continue
        r_value, p_value = pearson_stats(x, y)
        record = line.to_dict()
        record["最终IUP点数"] = int(len(pts))
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
    }
    return RTIUPResult(rows=final_rows, lines=valid_lines, plot_rows=final_rows, stats=stats)


def apply_same_carbon_iup_guard(result: RTIUPResult, config: RTIUPConfig) -> RTIUPResult:
    """Remove higher-unsaturation short-series points that invert at equal carbon.

    Curve-level IUP selection cannot judge two unrelated one-point series.  This
    final guard enforces the directly observable same-carbon ordering while
    allowing the two features to drift independently inside the locked 0.20 min
    RT window plus the locked 0.05 min IUP tolerance.
    """

    rows = result.rows.copy()
    remove_index: set[int] = set()
    tolerance = relative_rt_tolerance_min(config)
    for _, group in rows.groupby([C_CLASS, C_X], sort=False):
        medians = group.groupby(C_UNSAT, sort=True)[C_RT].median().sort_index()
        last_kept_rt: float | None = None
        for unsat, rt_value in medians.items():
            current_rt = float(rt_value)
            if last_kept_rt is not None and current_rt > last_kept_rt + tolerance:
                mask = group[C_UNSAT].astype(float).eq(float(unsat))
                remove_index.update(group.index[mask].tolist())
                continue
            last_kept_rt = current_rt if last_kept_rt is None else min(last_kept_rt, current_rt)

    if not remove_index:
        stats = dict(result.stats)
        stats["同碳数IUP顺序过滤行数"] = 0
        return RTIUPResult(rows=rows, lines=result.lines, plot_rows=result.plot_rows, stats=stats)

    kept_rows = rows.drop(index=sorted(remove_index)).copy()
    kept_plot_rows = result.plot_rows.drop(index=result.plot_rows.index.intersection(remove_index)).copy()
    stats = dict(result.stats)
    stats["同碳数IUP顺序过滤行数"] = int(len(remove_index))
    stats["最终保留行数"] = int(len(kept_rows))
    stats["拟合内点保留行数"] = int(kept_rows.get("同曲线0.2min内点", pd.Series(False, index=kept_rows.index)).sum())
    stats["点数不足保留行数"] = int(kept_rows.get("点数不足保留", pd.Series(False, index=kept_rows.index)).sum())
    stats["二次捞点保留行数"] = int(
        kept_rows.get("二次捞点", pd.Series(False, index=kept_rows.index)).fillna(False).astype(bool).sum()
    )
    stats["作图点数"] = int(len(kept_plot_rows))
    return RTIUPResult(rows=kept_rows, lines=result.lines, plot_rows=kept_plot_rows, stats=stats)


def title_suffix(grouping: str, method: str) -> str:
    grouping_code = {
        "详细链分类": "fixed_lcb_chain",
        "长链碱系列分类": "lcb_type",
        "总碳数总不饱和度分类": "total_cdb",
    }[grouping]
    return f"{grouping_code}_{method.lower()}_fit_ms2_only"


def make_contact_sheet(
    plot_dir: Path,
    result: RTIUPResult,
    out_path: Path,
    grouping: str,
    method: str,
    preferred_groups: list[str] | None = None,
) -> None:
    counts = result.plot_rows.groupby(C_CLASS).size().sort_values(ascending=False)
    groups: list[str] = []
    for group in preferred_groups or []:
        if group in counts.index and group not in groups:
            groups.append(group)
    for group in counts.index:
        if group not in groups:
            groups.append(str(group))
        if len(groups) >= 6:
            break

    fig, axes = plt.subplots(2, 3, figsize=(18, 12), facecolor="white")
    for ax, group in zip(axes.ravel(), groups):
        image_path = plot_dir / f"{safe_name(group)}_all.png"
        ax.imshow(plt.imread(image_path))
        ax.set_title(f"{group}_{title_suffix(grouping, method)}", fontsize=12, fontname="Arial", pad=8)
        ax.axis("off")
    for ax in axes.ravel()[len(groups):]:
        ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.01, wspace=0.06, hspace=0.10)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, facecolor="white")
    plt.close(fig)


def make_glccer_comparison(root: Path) -> Path:
    grouping = "详细链分类"
    target = "GlcCer(m15:2/x:y)"
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), facecolor="white")
    config = replace(RTIUPConfig(), **PLOT_STYLE)
    for ax, method in zip(np.atleast_1d(axes), ("ECN", "IUP")):
        table_path = root / grouping / f"{method}_表" / "最终保留明细.csv"
        retained = pd.read_csv(table_path)
        retained = retained[retained[C_CLASS].astype(str).eq(target)].copy()
        for unsat, group in retained.groupby(C_UNSAT, sort=True):
            color = color_for_unsaturation(unsat, config)
            if method == "ECN":
                ax.scatter(group[C_X], group[C_RT], s=54, color=color, alpha=0.8, edgecolors="none", label=f"{int(float(unsat))}")
            else:
                ax.scatter(
                    group[C_X],
                    group[C_RT],
                    s=62,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=1.8,
                    label=f"{int(float(unsat))} retained",
                )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(2.0)
        ax.spines["bottom"].set_linewidth(2.0)
        ax.tick_params(axis="both", direction="out", length=7, width=2.0)
        ax.set_xlabel("Carbon Number", fontname="Arial", fontweight="bold")
        ax.set_ylabel("Normalized Retention Time(min)", fontname="Arial", fontweight="bold")
        title_note = "" if method == "ECN" else " retained, unvalidated"
        ax.set_title(f"GlcCer(m15:2) - {method}{title_note}", fontsize=14, fontname="Arial")
        ax.legend(loc="lower left", frameon=False, prop={"family": "Arial", "size": 9})
        ax.margins(x=0.10, y=0.16)
        if method == "IUP":
            ax.text(
                0.98,
                0.98,
                f"n={len(retained)} retained\n0 fitted / plot-supported",
                ha="right",
                va="top",
                fontsize=10,
                color="#4A5568",
                transform=ax.transAxes,
            )
    fig.tight_layout()
    out_path = root / "汇总图" / "TOP3_GlcCer(m15_2)_ECN_IUP_对照.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, facecolor="white")
    plt.close(fig)
    return out_path


def build_three_status_tables(prepared: pd.DataFrame, result: RTIUPResult) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Classify every score-filtered candidate into exactly three fit states."""

    point_counts = (
        prepared[[C_CLASS, C_UNSAT, C_X, C_RT]]
        .dropna()
        .drop_duplicates()
        .groupby([C_CLASS, C_UNSAT], sort=False)
        .size()
        .rename("拟合独立点数")
        .reset_index()
    )
    audit = prepared.merge(point_counts, on=[C_CLASS, C_UNSAT], how="left")
    audit["拟合独立点数"] = pd.to_numeric(audit["拟合独立点数"], errors="coerce").fillna(0).astype(int)

    id_column = "合并后行ID"
    if id_column not in audit.columns or id_column not in result.rows.columns:
        raise KeyError("三状态明细需要合并后行ID来追踪候选结果")
    retained_ids = set(result.rows[id_column].dropna().tolist())
    audit["是否最终保留"] = audit[id_column].isin(retained_ids)
    audit["判定状态"] = np.select(
        [audit["拟合独立点数"].lt(3), audit["是否最终保留"]],
        ["点数不足未拟合", "拟合通过"],
        default="拟合未通过",
    )
    audit["状态说明"] = np.select(
        [
            audit["判定状态"].eq("点数不足未拟合") & audit["是否最终保留"],
            audit["判定状态"].eq("点数不足未拟合"),
            audit["判定状态"].eq("拟合通过"),
        ],
        [
            "独立点少于3，未拟合；按短序列保护规则保留",
            "独立点少于3，未拟合；未通过短序列保护规则",
            "拟合成功并通过当前ECN/IUP与0.20 min筛选",
        ],
        default="点数达到拟合要求，但未通过拟合、ECN/IUP或0.20 min筛选",
    )

    annotated_rows = result.rows.merge(point_counts, on=[C_CLASS, C_UNSAT], how="left")
    annotated_rows["拟合独立点数"] = pd.to_numeric(annotated_rows["拟合独立点数"], errors="coerce").fillna(0).astype(int)
    annotated_rows["判定状态"] = np.where(annotated_rows["拟合独立点数"].lt(3), "点数不足未拟合", "拟合通过")

    status_summary = (
        audit.groupby("判定状态", sort=False)
        .agg(候选行数=(id_column, "size"), 最终保留行数=("是否最终保留", "sum"))
        .reindex(["拟合通过", "拟合未通过", "点数不足未拟合"], fill_value=0)
        .reset_index()
    )
    return audit, annotated_rows, status_summary


def add_limited_plot_rt_shift(
    frame: pd.DataFrame,
    grouping: str,
    method: str,
    config: RTIUPConfig,
) -> pd.DataFrame:
    """Expose the bounded IUP whole-series correction already chosen by fitting."""

    out = frame.copy()
    if method == "IUP" and "IUP_RT漂移校正(min)" in out.columns:
        out["作图RT漂移校正(min)"] = pd.to_numeric(
            out["IUP_RT漂移校正(min)"], errors="coerce"
        ).fillna(0.0)
    else:
        out["作图RT漂移校正(min)"] = 0.0
    # IUP result rows already contain the corrected RT and retain the original
    # value in 原始归一化保留时间. Do not add a second display-only shift.
    out["作图归一化保留时间"] = pd.to_numeric(out[C_RT], errors="coerce")
    return out


def prepare_plot_with_limited_rt_shift(
    result: RTIUPResult,
    grouping: str,
    method: str,
    config: RTIUPConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Prepare already-corrected rows and lines for plotting and audit tables."""

    plot_rows = add_limited_plot_rt_shift(result.plot_rows, grouping, method, config)
    plot_rows[C_RT] = plot_rows["作图归一化保留时间"]
    plot_lines = result.lines.copy()
    line_table = result.lines.copy()
    line_table["作图RT漂移校正(min)"] = 0.0
    line_table["作图参数"] = line_table.get("参数")

    if method == "IUP" and "IUP_RT漂移校正(min)" in line_table.columns:
        line_table["作图RT漂移校正(min)"] = pd.to_numeric(
            line_table["IUP_RT漂移校正(min)"], errors="coerce"
        ).fillna(0.0)
        line_table["作图参数"] = line_table.get("参数")
    return plot_rows, plot_lines, line_table


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--input-sheet", default=None, help="Worksheet name for Excel input; CSV input ignores this option.")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--score-threshold", type=float, default=0.50)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_root.exists() and any(args.output_root.iterdir()):
        raise RuntimeError(f"Preview output directory is not empty: {args.output_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)

    if args.input.suffix.lower() == ".csv":
        raw = pd.read_csv(args.input, low_memory=False)
    else:
        raw = pd.read_excel(args.input, sheet_name=args.input_sheet or 0)
    groupings = prepare_groupings(raw, args.score_threshold)
    config = replace(RTIUPConfig(), **PLOT_STYLE)
    results: dict[tuple[str, str], RTIUPResult] = {}
    summary_records: list[dict[str, object]] = []

    preferred = {
        "详细链分类": ["Cer(m18:2/x:y)", "GlcCer(d18:1/x:y)", "SM(d18:1/x:y)", "Cer(d18:1/x:y)", "GlcCer(m18:2/x:y)", "Cer(d18:0/x:y)"],
        "长链碱系列分类": ["Cer(d)", "GlcCer(m)", "SM(d)", "GlcCer(t)", "DHCerP(d)", "So(m)"],
        "总碳数总不饱和度分类": ["Cer", "GlcCer", "SM", "LacCer", "DHCerP", "So"],
    }

    for grouping, prepared in groupings.items():
        print(f"[{grouping}] score-filtered input rows: {len(prepared)}", flush=True)
        print(f"[{grouping}] fitting ECN...", flush=True)
        ecn_result = fit_ecn_core(prepared, config)
        print(f"[{grouping}] independently fitting and coordinating IUP...", flush=True)
        iup_result = apply_same_carbon_iup_guard(fit_rt_iup(prepared, config), config)
        for method in ("ECN", "IUP"):
            result = ecn_result if method == "ECN" else iup_result
            results[(grouping, method)] = result
            audit, annotated_rows, status_summary = build_three_status_tables(prepared, result)
            audit = add_limited_plot_rt_shift(audit, grouping, method, config)
            annotated_rows = add_limited_plot_rt_shift(annotated_rows, grouping, method, config)
            plot_rows, plot_lines, line_table = prepare_plot_with_limited_rt_shift(result, grouping, method, config)
            plot_dir = args.output_root / grouping / f"{method}_单图"
            count = write_rt_iup_plots(plot_rows, plot_lines, plot_dir, config)
            table_dir = args.output_root / grouping / f"{method}_表"
            write_csv(annotated_rows, table_dir / "最终保留明细.csv")
            write_csv(line_table, table_dir / "有效拟合曲线.csv")
            write_csv(audit, table_dir / "全部候选三状态明细.csv")
            write_csv(status_summary, table_dir / "三状态汇总.csv")
            status_counts = audit["判定状态"].value_counts()
            summary_records.append(
                {
                    "分类方法": grouping,
                    "规则": method,
                    "单图数": count,
                    **result.stats,
                    "拟合通过行数": int(status_counts.get("拟合通过", 0)),
                    "拟合未通过行数": int(status_counts.get("拟合未通过", 0)),
                    "点数不足未拟合行数": int(status_counts.get("点数不足未拟合", 0)),
                }
            )
            summary_path = args.output_root / "汇总图" / f"TOP3_{grouping}_{method}_MS2_汇总图.png"
            make_contact_sheet(plot_dir, result, summary_path, grouping, method, preferred[grouping])
            print(f"[{grouping}] {method}: kept={len(result.rows)}, plots={count}", flush=True)

    comparison = make_glccer_comparison(args.output_root)
    summary = pd.DataFrame(summary_records)
    write_csv(summary, args.output_root / "TOP3_ECN_IUP_重新绘图汇总.csv")
    print(f"Output: {args.output_root}", flush=True)
    print(f"GlcCer comparison: {comparison}", flush=True)


if __name__ == "__main__":
    main()
