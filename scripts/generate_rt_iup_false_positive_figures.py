"""Create overview figures for RT/IUP false-positive filtering effects."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle


BASE_DIR = Path("E:/李奕晓数据备份/20240422-Agilent/原始数据")
OLD_WORKBOOKS = {
    "top5": BASE_DIR / "SphinGOlipID_top5_归一化保留时间_IUP最终拟合结果.xlsx",
    "top10": BASE_DIR / "SphinGOlipID_top10_归一化保留时间_IUP最终拟合结果.xlsx",
}
REFILTERED_WORKBOOK = BASE_DIR / "SphinGOlipID_RT_IUP_refiltered_top5_top10.xlsx"
OUT_DIR = BASE_DIR / "RT_IUP_false_positive_figures"

COL_ANNOTATION = "注释"
COL_CLASS = "细类"
COL_RANK = "候选排名"
COL_MATCH_SCORE = "匹配度分数"
COL_TOTAL_SCORE = "总分数"
COL_PLOT = "是否作图"
COL_KEEP_REASON = "保留原因"

COLORS = {
    "top5": "#4C78A8",
    "top10": "#F58518",
    "retained": "#4C78A8",
    "removed": "#E45756",
    "plot": "#54A24B",
    "neutral": "#B9C0C9",
}


@dataclass
class DatasetPayload:
    old_rows: pd.DataFrame
    new_rows: pd.DataFrame
    annotations: pd.DataFrame


def truthy_plot(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "是"}


def load_dataset(dataset: str) -> DatasetPayload:
    old_rows = pd.read_excel(OLD_WORKBOOKS[dataset], sheet_name="最终拟合结果")
    new_rows = pd.read_excel(REFILTERED_WORKBOOK, sheet_name=f"{dataset}重新过滤结果")
    new_annotations = set(new_rows[COL_ANNOTATION].astype(str))

    old_rows["__kept"] = old_rows[COL_ANNOTATION].astype(str).isin(new_annotations)
    annotations = (
        old_rows.groupby(COL_ANNOTATION)
        .agg(
            rows=(COL_ANNOTATION, "size"),
            lipid_class=(COL_CLASS, "first"),
            best_rank=(COL_RANK, "min"),
            best_match_score=(COL_MATCH_SCORE, "max"),
            best_total_score=(COL_TOTAL_SCORE, "max"),
            plotted_any=(COL_PLOT, lambda values: any(truthy_plot(value) for value in values)),
            fitted_rows=(COL_KEEP_REASON, lambda values: values.astype(str).str.contains("拟合成功", regex=False).sum()),
            short_rows=(COL_KEEP_REASON, lambda values: values.astype(str).str.contains("点数不足", regex=False).sum()),
            rescued_rows=(COL_KEEP_REASON, lambda values: values.astype(str).str.contains("捞点", regex=False).sum()),
            kept=("__kept", "any"),
        )
        .reset_index()
        .rename(columns={COL_ANNOTATION: "annotation"})
    )
    annotations["removed"] = ~annotations["kept"]
    annotations["dataset"] = dataset
    return DatasetPayload(old_rows=old_rows, new_rows=new_rows, annotations=annotations)


def save_figure(fig: plt.Figure, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.4,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def annotate_barh(ax: plt.Axes, bars, fmt: str = "{:,.0f}") -> None:
    for bar in bars:
        width = bar.get_width()
        ax.text(width, bar.get_y() + bar.get_height() / 2, "  " + fmt.format(width), va="center", ha="left", fontsize=8)


def figure_workflow(payloads: dict[str, DatasetPayload]) -> Path:
    rows = []
    for dataset, payload in payloads.items():
        rows.extend(
            [
                {"dataset": dataset, "stage": "Candidate rows", "count": len(payload.old_rows), "panel": "raw"},
                {"dataset": dataset, "stage": "Unique annotations", "count": payload.annotations.shape[0], "panel": "final"},
                {"dataset": dataset, "stage": "RT/IUP retained", "count": int(payload.annotations["kept"].sum()), "panel": "final"},
                {"dataset": dataset, "stage": "Plot-supported", "count": int(payload.new_rows[COL_PLOT].map(truthy_plot).sum()), "panel": "final"},
            ]
        )
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), gridspec_kw={"width_ratios": [1.05, 1.65]})

    raw = df[df["panel"].eq("raw")]
    y = np.arange(len(raw))
    bars = axes[0].barh(y, raw["count"], color=[COLORS[d] for d in raw["dataset"]], alpha=0.92)
    axes[0].set_yticks(y, [f"{d} candidate rows" for d in raw["dataset"]])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Rows")
    axes[0].set_title("Raw candidate records")
    annotate_barh(axes[0], bars)
    axes[0].set_xlim(0, raw["count"].max() * 1.22)

    final = df[df["panel"].eq("final")]
    stages = ["Unique annotations", "RT/IUP retained", "Plot-supported"]
    x = np.arange(len(stages))
    width = 0.36
    for offset, dataset in [(-width / 2, "top5"), (width / 2, "top10")]:
        vals = [int(final[(final["dataset"].eq(dataset)) & (final["stage"].eq(stage))]["count"].iloc[0]) for stage in stages]
        bars = axes[1].bar(x + offset, vals, width=width, color=COLORS[dataset], label=dataset, alpha=0.92)
        for bar, value in zip(bars, vals):
            axes[1].text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}", ha="center", va="bottom", fontsize=8, rotation=90)
    axes[1].set_xticks(x, stages, rotation=18, ha="right")
    axes[1].set_ylabel("Annotations / plot points")
    axes[1].set_title("RT/IUP filtering workflow")
    axes[1].legend(frameon=False)
    axes[1].set_ylim(0, final["count"].max() * 1.26)
    fig.suptitle("Filtering overview", y=1.03, fontsize=12, fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, "01_filtering_workflow.png")


def figure_retention_stacked(payloads: dict[str, DatasetPayload]) -> Path:
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    datasets = ["top5", "top10"]
    x = np.arange(len(datasets))
    retained = np.array([payloads[d].annotations["kept"].sum() for d in datasets], dtype=float)
    removed = np.array([payloads[d].annotations["removed"].sum() for d in datasets], dtype=float)
    total = retained + removed
    ax.bar(x, retained / total * 100, color=COLORS["retained"], label="Retained")
    ax.bar(x, removed / total * 100, bottom=retained / total * 100, color=COLORS["removed"], label="RT/IUP removed")
    for i, dataset in enumerate(datasets):
        ax.text(i, retained[i] / total[i] * 50, f"{int(retained[i]):,}\nretained", ha="center", va="center", color="white", fontsize=8)
        ax.text(
            i,
            101.2,
            f"{int(removed[i]):,}\n({removed[i] / total[i] * 100:.1f}%)",
            ha="center",
            va="bottom",
            color=COLORS["removed"],
            fontsize=8,
        )
    ax.set_xticks(x, datasets)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Unique annotations (%)")
    ax.set_title("Putative false positives removed by RT/IUP")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2)
    fig.tight_layout()
    return save_figure(fig, "02_rt_iup_removed_fraction.png")


def class_stats(payloads: dict[str, DatasetPayload]) -> pd.DataFrame:
    frames = []
    for dataset, payload in payloads.items():
        grouped = (
            payload.annotations.groupby("lipid_class")["removed"]
            .agg(total_annotations="count", removed_annotations="sum")
            .reset_index()
        )
        grouped["retained_annotations"] = grouped["total_annotations"] - grouped["removed_annotations"]
        grouped["removed_rate"] = grouped["removed_annotations"] / grouped["total_annotations"]
        grouped["dataset"] = dataset
        frames.append(grouped)
    return pd.concat(frames, ignore_index=True)


def figure_class_removed_top20(stats: pd.DataFrame) -> Path:
    top_classes = (
        stats.groupby("lipid_class")["removed_annotations"]
        .sum()
        .sort_values(ascending=False)
        .head(20)
        .index.tolist()
    )
    plot = stats[stats["lipid_class"].isin(top_classes)].pivot(index="lipid_class", columns="dataset", values="removed_annotations").fillna(0)
    plot["combined"] = plot.sum(axis=1)
    plot = plot.sort_values("combined")
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    y = np.arange(len(plot))
    width = 0.38
    ax.barh(y - width / 2, plot.get("top5", 0), height=width, color=COLORS["top5"], label="top5")
    ax.barh(y + width / 2, plot.get("top10", 0), height=width, color=COLORS["top10"], label="top10")
    ax.set_yticks(y, plot.index)
    ax.set_xlabel("Removed annotations")
    ax.set_title("Top lipid subclasses affected by RT/IUP filtering")
    ax.legend(frameon=False)
    fig.tight_layout()
    return save_figure(fig, "03_removed_annotations_by_class_top20.png")


def figure_class_bubble(stats: pd.DataFrame) -> Path:
    candidates = stats[stats["total_annotations"].ge(5)].copy()
    top_classes = (
        candidates.assign(score=candidates["removed_annotations"] * (0.5 + candidates["removed_rate"]))
        .groupby("lipid_class")["score"]
        .max()
        .sort_values(ascending=False)
        .head(22)
        .index.tolist()
    )
    plot = candidates[candidates["lipid_class"].isin(top_classes)].copy()
    class_order = (
        plot.groupby("lipid_class")["removed_rate"]
        .max()
        .sort_values()
        .index.tolist()
    )
    y_lookup = {name: idx for idx, name in enumerate(class_order)}
    x_lookup = {"top5": 0, "top10": 1}
    fig, ax = plt.subplots(figsize=(5.4, 6.2))
    sizes = 45 + plot["removed_annotations"].astype(float) * 28
    scatter = ax.scatter(
        plot["dataset"].map(x_lookup),
        plot["lipid_class"].map(y_lookup),
        s=sizes,
        c=plot["removed_rate"] * 100,
        cmap="Reds",
        edgecolors="#333333",
        linewidths=0.5,
        alpha=0.88,
    )
    ax.set_xticks([0, 1], ["top5", "top10"])
    ax.set_yticks(range(len(class_order)), class_order)
    ax.set_xlim(-0.45, 1.45)
    ax.set_xlabel("Dataset")
    ax.set_title("Subclass-level removal rate")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Removed (%)")
    for removed in [2, 5, 10]:
        ax.scatter([], [], s=45 + removed * 28, c="white", edgecolors="#333333", label=f"{removed} removed")
    ax.legend(frameon=False, loc="lower right", title="Bubble size")
    fig.tight_layout()
    return save_figure(fig, "04_removed_rate_bubble_by_class.png")


def figure_score_distribution(payloads: dict[str, DatasetPayload]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.8), sharey=False)
    metrics = [("best_total_score", "Best total score"), ("best_match_score", "Best match score")]
    for ax, (metric, label) in zip(axes, metrics):
        data = []
        labels = []
        colors = []
        for dataset in ["top5", "top10"]:
            ann = payloads[dataset].annotations
            for status, status_label, color in [(False, "Retained", COLORS[dataset]), (True, "Removed", COLORS["removed"])]:
                data.append(ann.loc[ann["removed"].eq(status), metric].dropna().to_numpy())
                labels.append(f"{dataset}\n{status_label}")
                colors.append(color)
        box = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
            patch.set_edgecolor("#333333")
        for median in box["medians"]:
            median.set_color("#111111")
            median.set_linewidth(1.4)
        ax.set_title(label)
        ax.tick_params(axis="x", rotation=0)
        ax.grid(axis="y", color="#E8E8E8", linewidth=0.8)
    fig.suptitle("MS2 score distributions before vs after RT/IUP filtering", y=1.04, fontsize=12, fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, "05_score_distribution_retained_vs_removed.png")


def figure_overlap(payloads: dict[str, DatasetPayload]) -> Path:
    removed = {dataset: set(payload.annotations.loc[payload.annotations["removed"], "annotation"].astype(str)) for dataset, payload in payloads.items()}
    overlap = len(removed["top5"] & removed["top10"])
    top5_only = len(removed["top5"] - removed["top10"])
    top10_only = len(removed["top10"] - removed["top5"])

    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    ax.add_patch(Circle((0.43, 0.5), 0.28, color=COLORS["top5"], alpha=0.42, lw=2, ec=COLORS["top5"]))
    ax.add_patch(Circle((0.62, 0.5), 0.28, color=COLORS["top10"], alpha=0.42, lw=2, ec=COLORS["top10"]))
    ax.text(0.31, 0.50, f"{top5_only}", ha="center", va="center", fontsize=20, fontweight="bold")
    ax.text(0.525, 0.50, f"{overlap}", ha="center", va="center", fontsize=20, fontweight="bold")
    ax.text(0.74, 0.50, f"{top10_only}", ha="center", va="center", fontsize=20, fontweight="bold")
    ax.text(0.35, 0.82, f"top5 removed\n{len(removed['top5'])}", ha="center", va="center", fontsize=10, color=COLORS["top5"])
    ax.text(0.70, 0.82, f"top10 removed\n{len(removed['top10'])}", ha="center", va="center", fontsize=10, color=COLORS["top10"])
    ax.text(0.525, 0.18, "Shared RT/IUP-filtered annotations", ha="center", va="center", fontsize=10)
    ax.set_xlim(0.1, 0.95)
    ax.set_ylim(0.12, 0.9)
    ax.axis("off")
    ax.set_title("Overlap of removed annotations", fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, "06_removed_annotation_overlap.png")


def write_summary_tables(payloads: dict[str, DatasetPayload], stats: pd.DataFrame) -> None:
    summary_rows = []
    for dataset, payload in payloads.items():
        ann = payload.annotations
        summary_rows.append(
            {
                "dataset": dataset,
                "candidate_rows": len(payload.old_rows),
                "unique_annotations": len(ann),
                "rt_iup_retained": int(ann["kept"].sum()),
                "rt_iup_removed": int(ann["removed"].sum()),
                "rt_iup_removed_rate": float(ann["removed"].mean()),
                "plot_supported": int(payload.new_rows[COL_PLOT].map(truthy_plot).sum()),
            }
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "summary_counts.csv", index=False, encoding="utf-8-sig")
    stats.to_csv(OUT_DIR / "class_removal_stats.csv", index=False, encoding="utf-8-sig")
    pd.concat([payload.annotations for payload in payloads.values()], ignore_index=True).to_csv(
        OUT_DIR / "annotation_removal_stats.csv", index=False, encoding="utf-8-sig"
    )


def main() -> None:
    setup_style()
    payloads = {dataset: load_dataset(dataset) for dataset in ["top5", "top10"]}
    stats = class_stats(payloads)
    paths = [
        figure_workflow(payloads),
        figure_retention_stacked(payloads),
        figure_class_removed_top20(stats),
        figure_class_bubble(stats),
        figure_score_distribution(payloads),
        figure_overlap(payloads),
    ]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
