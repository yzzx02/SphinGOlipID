"""Compare RT/IUP filtering effects across top3, top5, and top10 outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path("E:/李奕晓数据备份/20240422-Agilent/原始数据")
CHAIN_DEDUP_WORKBOOK = BASE_DIR / "SphinGOlipID_最终鉴定结果_链信息去重.xlsx"
TOP3_WORKBOOK = BASE_DIR / "SphinGOlipID_top3_归一化保留时间_IUP最终拟合结果.xlsx"
REFILTERED_WORKBOOK = BASE_DIR / "SphinGOlipID_RT_IUP_refiltered_top5_top10.xlsx"
OUT_DIR = BASE_DIR / "RT_IUP_top3_top5_top10_comparison_figures"

DATASETS = {
    "top3": {
        "source_workbook": TOP3_WORKBOOK,
        "source_sheet": "top3过滤前唯一注释",
        "result_workbook": TOP3_WORKBOOK,
        "result_sheet": "最终拟合结果",
    },
    "top5": {
        "source_workbook": CHAIN_DEDUP_WORKBOOK,
        "source_sheet": "top5最终鉴定",
        "result_workbook": REFILTERED_WORKBOOK,
        "result_sheet": "top5重新过滤结果",
    },
    "top10": {
        "source_workbook": CHAIN_DEDUP_WORKBOOK,
        "source_sheet": "top10最终鉴定",
        "result_workbook": REFILTERED_WORKBOOK,
        "result_sheet": "top10重新过滤结果",
    },
}

COL_ANNOTATION = "注释"
COL_CLASS = "细类"
COL_RANK = "候选排名"
COL_MATCH_SCORE = "匹配度分数"
COL_TOTAL_SCORE = "总分数"
COL_PLOT = "是否作图"

ORDER = ["top3", "top5", "top10"]
COLORS = {
    "top3": "#54A24B",
    "top5": "#4C78A8",
    "top10": "#F58518",
    "removed": "#E45756",
    "neutral": "#B9C0C9",
}


@dataclass
class DatasetPayload:
    source_rows: pd.DataFrame
    result_rows: pd.DataFrame
    annotations: pd.DataFrame


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "是"}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.3,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def load_dataset(dataset: str) -> DatasetPayload:
    spec = DATASETS[dataset]
    source = pd.read_excel(spec["source_workbook"], sheet_name=spec["source_sheet"])
    result = pd.read_excel(spec["result_workbook"], sheet_name=spec["result_sheet"])
    retained = set(result[COL_ANNOTATION].astype(str))
    source["__kept"] = source[COL_ANNOTATION].astype(str).isin(retained)

    annotations = (
        source.groupby(COL_ANNOTATION)
        .agg(
            rows=(COL_ANNOTATION, "size"),
            lipid_class=(COL_CLASS, "first"),
            best_rank=(COL_RANK, "min"),
            best_match_score=(COL_MATCH_SCORE, "max"),
            best_total_score=(COL_TOTAL_SCORE, "max"),
            kept=("__kept", "any"),
        )
        .reset_index()
        .rename(columns={COL_ANNOTATION: "annotation"})
    )
    annotations["removed"] = ~annotations["kept"]
    annotations["dataset"] = dataset
    return DatasetPayload(source_rows=source, result_rows=result, annotations=annotations)


def save_figure(fig: plt.Figure, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def figure_workflow(payloads: dict[str, DatasetPayload]) -> Path:
    stages = ["Input annotations", "RT/IUP retained", "Plot-supported"]
    values = {
        dataset: [
            len(payload.annotations),
            int(payload.annotations["kept"].sum()),
            int(payload.result_rows[COL_PLOT].map(truthy).sum()),
        ]
        for dataset, payload in payloads.items()
    }
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    x = np.arange(len(stages))
    width = 0.24
    offsets = [-width, 0, width]
    for offset, dataset in zip(offsets, ORDER):
        bars = ax.bar(x + offset, values[dataset], width=width, color=COLORS[dataset], label=dataset, alpha=0.92)
        for bar, value in zip(bars, values[dataset]):
            ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}", ha="center", va="bottom", fontsize=8, rotation=90)
    ax.set_xticks(x, stages, rotation=12, ha="right")
    ax.set_ylabel("Annotations / plot points")
    ax.set_title("TopN RT/IUP filtering workflow")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    ax.set_ylim(0, max(max(v) for v in values.values()) * 1.25)
    fig.tight_layout()
    return save_figure(fig, "01_topn_filtering_workflow.png")


def figure_removed_fraction(payloads: dict[str, DatasetPayload]) -> Path:
    retained = np.array([payloads[d].annotations["kept"].sum() for d in ORDER], dtype=float)
    removed = np.array([payloads[d].annotations["removed"].sum() for d in ORDER], dtype=float)
    total = retained + removed
    x = np.arange(len(ORDER))

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.bar(x, retained / total * 100, color=[COLORS[d] for d in ORDER], label="Retained", alpha=0.9)
    ax.bar(x, removed / total * 100, bottom=retained / total * 100, color=COLORS["removed"], label="RT/IUP removed")
    for i, dataset in enumerate(ORDER):
        ax.text(i, 48, f"{int(retained[i]):,}\nretained", ha="center", va="center", color="white", fontsize=8)
        ax.text(
            i,
            101.2,
            f"{int(removed[i]):,}\n({removed[i] / total[i] * 100:.1f}%)",
            ha="center",
            va="bottom",
            color=COLORS["removed"],
            fontsize=8,
        )
    ax.set_xticks(x, ORDER)
    ax.set_ylim(0, 113)
    ax.set_ylabel("Input annotations (%)")
    ax.set_title("RT/IUP-filtered annotations by topN")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2)
    fig.tight_layout()
    return save_figure(fig, "02_topn_removed_fraction.png")


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
        .head(22)
        .index.tolist()
    )
    plot = stats[stats["lipid_class"].isin(top_classes)].pivot(index="lipid_class", columns="dataset", values="removed_annotations").fillna(0)
    plot["combined"] = plot.sum(axis=1)
    plot = plot.sort_values("combined")
    fig, ax = plt.subplots(figsize=(6.8, 6.3))
    y = np.arange(len(plot))
    width = 0.25
    for offset, dataset in zip([-width, 0, width], ORDER):
        ax.barh(y + offset, plot.get(dataset, 0), height=width, color=COLORS[dataset], label=dataset)
    ax.set_yticks(y, plot.index)
    ax.set_xlabel("Removed annotations")
    ax.set_title("Lipid subclasses most affected by RT/IUP filtering")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    return save_figure(fig, "03_topn_removed_annotations_by_class.png")


def figure_class_bubble(stats: pd.DataFrame) -> Path:
    candidates = stats[stats["total_annotations"].ge(5)].copy()
    top_classes = (
        candidates.assign(score=candidates["removed_annotations"] * (0.5 + candidates["removed_rate"]))
        .groupby("lipid_class")["score"]
        .max()
        .sort_values(ascending=False)
        .head(24)
        .index.tolist()
    )
    plot = candidates[candidates["lipid_class"].isin(top_classes)].copy()
    class_order = plot.groupby("lipid_class")["removed_rate"].max().sort_values().index.tolist()
    x_lookup = {dataset: idx for idx, dataset in enumerate(ORDER)}
    y_lookup = {name: idx for idx, name in enumerate(class_order)}

    fig, ax = plt.subplots(figsize=(5.8, 6.4))
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
    ax.set_xticks(range(len(ORDER)), ORDER)
    ax.set_yticks(range(len(class_order)), class_order)
    ax.set_xlim(-0.5, len(ORDER) - 0.5)
    ax.set_xlabel("Dataset")
    ax.set_title("Subclass removal rate by topN")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Removed (%)")
    for removed in [2, 5, 10]:
        ax.scatter([], [], s=45 + removed * 28, c="white", edgecolors="#333333", label=f"{removed} removed")
    ax.legend(frameon=False, loc="lower right", title="Bubble size")
    fig.tight_layout()
    return save_figure(fig, "04_topn_removed_rate_bubble_by_class.png")


def figure_score_distribution(payloads: dict[str, DatasetPayload]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
    metrics = [("best_total_score", "Best total score"), ("best_match_score", "Best match score")]
    for ax, (metric, title) in zip(axes, metrics):
        data = []
        labels = []
        colors = []
        for dataset in ORDER:
            ann = payloads[dataset].annotations
            for status, suffix, color in [(False, "retained", COLORS[dataset]), (True, "removed", COLORS["removed"])]:
                data.append(ann.loc[ann["removed"].eq(status), metric].dropna().to_numpy())
                labels.append(f"{dataset}\n{suffix}")
                colors.append(color)
        box = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
            patch.set_edgecolor("#333333")
        for median in box["medians"]:
            median.set_color("#111111")
            median.set_linewidth(1.4)
        ax.set_title(title)
        ax.grid(axis="y", color="#E8E8E8", linewidth=0.8)
    fig.suptitle("Score distributions of retained and RT/IUP-filtered annotations", y=1.04, fontsize=12, fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, "05_topn_score_distribution.png")


def figure_removed_overlap(payloads: dict[str, DatasetPayload]) -> Path:
    removed_sets = {dataset: set(payload.annotations.loc[payload.annotations["removed"], "annotation"].astype(str)) for dataset, payload in payloads.items()}
    union = sorted(set().union(*removed_sets.values()))
    combo_counts: dict[tuple[str, ...], int] = {}
    for annotation in union:
        combo = tuple(dataset for dataset in ORDER if annotation in removed_sets[dataset])
        combo_counts[combo] = combo_counts.get(combo, 0) + 1
    combos = sorted(combo_counts, key=lambda combo: (combo_counts[combo], len(combo)), reverse=True)
    counts = [combo_counts[combo] for combo in combos]

    fig = plt.figure(figsize=(7.0, 4.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 3.0], height_ratios=[2.6, 1.2], hspace=0.05, wspace=0.18)
    ax_sets = fig.add_subplot(gs[:, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_matrix = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    set_sizes = [len(removed_sets[d]) for d in ORDER]
    y = np.arange(len(ORDER))
    ax_sets.barh(y, set_sizes, color=[COLORS[d] for d in ORDER])
    ax_sets.set_yticks(y, ORDER)
    ax_sets.invert_yaxis()
    ax_sets.set_xlabel("Set size")
    ax_sets.set_title("Removed")
    for idx, value in enumerate(set_sizes):
        ax_sets.text(value, idx, f" {value}", va="center", fontsize=8)

    x = np.arange(len(combos))
    ax_bar.bar(x, counts, color="#555555")
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_title("Overlap of RT/IUP-filtered annotations")
    for xi, value in zip(x, counts):
        ax_bar.text(xi, value, str(value), ha="center", va="bottom", fontsize=8)
    ax_bar.tick_params(axis="x", labelbottom=False)

    ax_matrix.set_ylim(-0.5, len(ORDER) - 0.5)
    ax_matrix.set_yticks(range(len(ORDER)), ORDER)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xticks(x, ["+".join(combo) for combo in combos], rotation=45, ha="right")
    for xi, combo in zip(x, combos):
        active_y = []
        for yi, dataset in enumerate(ORDER):
            active = dataset in combo
            ax_matrix.scatter(xi, yi, s=58 if active else 22, color=COLORS[dataset] if active else "#D0D0D0", zorder=3)
            if active:
                active_y.append(yi)
        if len(active_y) > 1:
            ax_matrix.plot([xi, xi], [min(active_y), max(active_y)], color="#333333", linewidth=1.2, zorder=2)
    ax_matrix.spines["left"].set_visible(False)
    ax_matrix.spines["bottom"].set_visible(False)
    ax_matrix.tick_params(axis="y", length=0)
    ax_matrix.tick_params(axis="x", length=0)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.24, top=0.88, wspace=0.22, hspace=0.05)
    return save_figure(fig, "06_topn_removed_annotation_overlap.png")


def write_summary_tables(payloads: dict[str, DatasetPayload], stats: pd.DataFrame) -> None:
    summary = []
    for dataset, payload in payloads.items():
        ann = payload.annotations
        summary.append(
            {
                "dataset": dataset,
                "input_annotations": len(ann),
                "rt_iup_retained": int(ann["kept"].sum()),
                "rt_iup_removed": int(ann["removed"].sum()),
                "rt_iup_removed_rate": float(ann["removed"].mean()),
                "plot_supported": int(payload.result_rows[COL_PLOT].map(truthy).sum()),
            }
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary).to_csv(OUT_DIR / "summary_counts.csv", index=False, encoding="utf-8-sig")
    stats.to_csv(OUT_DIR / "class_removal_stats.csv", index=False, encoding="utf-8-sig")
    pd.concat([payload.annotations for payload in payloads.values()], ignore_index=True).to_csv(
        OUT_DIR / "annotation_removal_stats.csv", index=False, encoding="utf-8-sig"
    )


def main() -> None:
    setup_style()
    payloads = {dataset: load_dataset(dataset) for dataset in ORDER}
    stats = class_stats(payloads)
    paths = [
        figure_workflow(payloads),
        figure_removed_fraction(payloads),
        figure_class_removed_top20(stats),
        figure_class_bubble(stats),
        figure_score_distribution(payloads),
        figure_removed_overlap(payloads),
    ]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
