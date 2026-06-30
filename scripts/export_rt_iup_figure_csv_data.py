"""Export CSV data behind all RT/IUP overview and comparison figures."""

from __future__ import annotations

from itertools import chain
from pathlib import Path

import pandas as pd


BASE_DIR = Path("E:/李奕晓数据备份/20240422-Agilent/原始数据")
OUT_DIR = BASE_DIR / "RT_IUP_all_figure_data_csv"

CHAIN_DEDUP_WORKBOOK = BASE_DIR / "SphinGOlipID_最终鉴定结果_链信息去重.xlsx"
TOP3_WORKBOOK = BASE_DIR / "SphinGOlipID_top3_归一化保留时间_IUP最终拟合结果.xlsx"
TOP5_TOP10_WORKBOOK = BASE_DIR / "SphinGOlipID_RT_IUP_refiltered_top5_top10.xlsx"
OLD_WORKBOOKS = {
    "top5": BASE_DIR / "SphinGOlipID_top5_归一化保留时间_IUP最终拟合结果.xlsx",
    "top10": BASE_DIR / "SphinGOlipID_top10_归一化保留时间_IUP最终拟合结果.xlsx",
}

COL_ANNOTATION = "注释"
COL_CLASS = "细类"
COL_RANK = "候选排名"
COL_MATCH_SCORE = "匹配度分数"
COL_TOTAL_SCORE = "总分数"
COL_PLOT = "是否作图"
COL_KEEP_REASON = "保留原因"

TOPN_ORDER = ["top3", "top5", "top10"]


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "是"}


def write_csv(df: pd.DataFrame, relative_path: str) -> Path:
    path = OUT_DIR / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def load_topn_source_and_result(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if dataset == "top3":
        source = pd.read_excel(TOP3_WORKBOOK, sheet_name="top3过滤前唯一注释")
        result = pd.read_excel(TOP3_WORKBOOK, sheet_name="最终拟合结果")
    elif dataset == "top5":
        source = pd.read_excel(CHAIN_DEDUP_WORKBOOK, sheet_name="top5最终鉴定")
        result = pd.read_excel(TOP5_TOP10_WORKBOOK, sheet_name="top5重新过滤结果")
    elif dataset == "top10":
        source = pd.read_excel(CHAIN_DEDUP_WORKBOOK, sheet_name="top10最终鉴定")
        result = pd.read_excel(TOP5_TOP10_WORKBOOK, sheet_name="top10重新过滤结果")
    else:
        raise ValueError(f"Unknown dataset: {dataset}")
    return source, result


def annotation_stats_from_source(source: pd.DataFrame, retained_annotations: set[str], dataset: str) -> pd.DataFrame:
    data = source.copy()
    data["__kept"] = data[COL_ANNOTATION].astype(str).isin(retained_annotations)
    out = (
        data.groupby(COL_ANNOTATION)
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
    out["removed"] = ~out["kept"]
    out["dataset"] = dataset
    return out


def load_topn_payloads() -> dict[str, dict[str, pd.DataFrame]]:
    payloads: dict[str, dict[str, pd.DataFrame]] = {}
    for dataset in TOPN_ORDER:
        source, result = load_topn_source_and_result(dataset)
        retained = set(result[COL_ANNOTATION].astype(str))
        annotations = annotation_stats_from_source(source, retained, dataset)
        payloads[dataset] = {"source": source, "result": result, "annotations": annotations}
    return payloads


def load_false_positive_payloads() -> dict[str, dict[str, pd.DataFrame]]:
    payloads: dict[str, dict[str, pd.DataFrame]] = {}
    for dataset in ["top5", "top10"]:
        old_rows = pd.read_excel(OLD_WORKBOOKS[dataset], sheet_name="最终拟合结果")
        _, result = load_topn_source_and_result(dataset)
        retained = set(result[COL_ANNOTATION].astype(str))
        old_rows["__kept"] = old_rows[COL_ANNOTATION].astype(str).isin(retained)
        annotations = (
            old_rows.groupby(COL_ANNOTATION)
            .agg(
                rows=(COL_ANNOTATION, "size"),
                lipid_class=(COL_CLASS, "first"),
                best_rank=(COL_RANK, "min"),
                best_match_score=(COL_MATCH_SCORE, "max"),
                best_total_score=(COL_TOTAL_SCORE, "max"),
                plotted_any=(COL_PLOT, lambda values: any(truthy(value) for value in values)),
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
        payloads[dataset] = {"old_rows": old_rows, "result": result, "annotations": annotations}
    return payloads


def class_stats(payloads: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    frames = []
    for dataset, payload in payloads.items():
        grouped = (
            payload["annotations"].groupby("lipid_class")["removed"]
            .agg(total_annotations="count", removed_annotations="sum")
            .reset_index()
        )
        grouped["retained_annotations"] = grouped["total_annotations"] - grouped["removed_annotations"]
        grouped["removed_rate"] = grouped["removed_annotations"] / grouped["total_annotations"]
        grouped["dataset"] = dataset
        frames.append(grouped)
    return pd.concat(frames, ignore_index=True)


def workflow_counts(payloads: dict[str, dict[str, pd.DataFrame]], include_candidate_rows: bool = False) -> pd.DataFrame:
    rows = []
    for dataset, payload in payloads.items():
        if include_candidate_rows:
            rows.append({"dataset": dataset, "stage": "Candidate rows", "count": len(payload["old_rows"])})
        rows.extend(
            [
                {"dataset": dataset, "stage": "Input annotations", "count": len(payload["annotations"])},
                {"dataset": dataset, "stage": "RT/IUP retained", "count": int(payload["annotations"]["kept"].sum())},
                {"dataset": dataset, "stage": "Plot-supported", "count": int(payload["result"][COL_PLOT].map(truthy).sum())},
            ]
        )
    return pd.DataFrame(rows)


def removed_fraction(payloads: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    rows = []
    for dataset, payload in payloads.items():
        annotations = payload["annotations"]
        retained = int(annotations["kept"].sum())
        removed = int(annotations["removed"].sum())
        total = retained + removed
        rows.append(
            {
                "dataset": dataset,
                "input_annotations": total,
                "retained": retained,
                "removed": removed,
                "retained_rate": retained / total if total else 0,
                "removed_rate": removed / total if total else 0,
            }
        )
    return pd.DataFrame(rows)


def class_top_data(stats: pd.DataFrame, n: int = 25) -> pd.DataFrame:
    top_classes = (
        stats.groupby("lipid_class")["removed_annotations"]
        .sum()
        .sort_values(ascending=False)
        .head(n)
        .index
    )
    return stats[stats["lipid_class"].isin(top_classes)].sort_values(["lipid_class", "dataset"]).reset_index(drop=True)


def class_bubble_data(stats: pd.DataFrame, n: int = 25) -> pd.DataFrame:
    candidates = stats[stats["total_annotations"].ge(5)].copy()
    top_classes = (
        candidates.assign(score=candidates["removed_annotations"] * (0.5 + candidates["removed_rate"]))
        .groupby("lipid_class")["score"]
        .max()
        .sort_values(ascending=False)
        .head(n)
        .index
    )
    return candidates[candidates["lipid_class"].isin(top_classes)].sort_values(["lipid_class", "dataset"]).reset_index(drop=True)


def score_distribution_data(payloads: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    frames = []
    for dataset, payload in payloads.items():
        data = payload["annotations"][["dataset", "annotation", "removed", "best_total_score", "best_match_score", "best_rank", "lipid_class"]].copy()
        data["status"] = data["removed"].map({True: "removed", False: "retained"})
        frames.append(data)
    return pd.concat(frames, ignore_index=True)


def removed_overlap_data(payloads: dict[str, dict[str, pd.DataFrame]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    removed_sets = {
        dataset: set(payload["annotations"].loc[payload["annotations"]["removed"], "annotation"].astype(str))
        for dataset, payload in payloads.items()
    }
    all_annotations = sorted(set(chain.from_iterable(removed_sets.values())))
    membership_rows = []
    for annotation in all_annotations:
        row = {"annotation": annotation}
        for dataset in payloads:
            row[f"removed_in_{dataset}"] = annotation in removed_sets[dataset]
        row["combination"] = "+".join(dataset for dataset in payloads if row[f"removed_in_{dataset}"])
        membership_rows.append(row)
    membership = pd.DataFrame(membership_rows)
    intersections = (
        membership.groupby("combination")
        .size()
        .reset_index(name="intersection_size")
        .sort_values("intersection_size", ascending=False)
        .reset_index(drop=True)
    )
    return intersections, membership


def export_rt_iup_tables() -> list[Path]:
    written: list[Path] = []
    table_specs = {
        "top3": (TOP3_WORKBOOK, "top3过滤前唯一注释", "最终拟合结果", "拟合曲线结果"),
        "top5": (TOP5_TOP10_WORKBOOK, None, "top5重新过滤结果", "top5重新拟合曲线"),
        "top10": (TOP5_TOP10_WORKBOOK, None, "top10重新过滤结果", "top10重新拟合曲线"),
    }
    for dataset, (workbook, source_sheet, result_sheet, lines_sheet) in table_specs.items():
        if source_sheet is not None:
            written.append(write_csv(pd.read_excel(workbook, sheet_name=source_sheet), f"rt_iup_tables/{dataset}_source_annotations.csv"))
        result = pd.read_excel(workbook, sheet_name=result_sheet)
        lines = pd.read_excel(workbook, sheet_name=lines_sheet)
        written.append(write_csv(result, f"rt_iup_tables/{dataset}_final_results.csv"))
        written.append(write_csv(lines, f"rt_iup_tables/{dataset}_fit_lines.csv"))
        if COL_PLOT in result.columns:
            written.append(write_csv(result[result[COL_PLOT].map(truthy)].copy(), f"rt_iup_tables/{dataset}_plot_points.csv"))
    return written


def export_topn_comparison() -> list[Path]:
    payloads = load_topn_payloads()
    stats = class_stats(payloads)
    intersections, membership = removed_overlap_data(payloads)
    written = [
        write_csv(workflow_counts(payloads), "top3_top5_top10_comparison/01_workflow_counts.csv"),
        write_csv(removed_fraction(payloads), "top3_top5_top10_comparison/02_removed_fraction.csv"),
        write_csv(class_top_data(stats, 25), "top3_top5_top10_comparison/03_class_removed_top_data.csv"),
        write_csv(class_bubble_data(stats, 25), "top3_top5_top10_comparison/04_class_bubble_data.csv"),
        write_csv(score_distribution_data(payloads), "top3_top5_top10_comparison/05_score_distribution_data.csv"),
        write_csv(intersections, "top3_top5_top10_comparison/06_removed_overlap_intersections.csv"),
        write_csv(membership, "top3_top5_top10_comparison/06_removed_overlap_membership.csv"),
        write_csv(removed_fraction(payloads), "top3_top5_top10_comparison/summary_counts.csv"),
        write_csv(stats, "top3_top5_top10_comparison/class_removal_stats.csv"),
        write_csv(pd.concat([payload["annotations"] for payload in payloads.values()], ignore_index=True), "top3_top5_top10_comparison/annotation_removal_stats.csv"),
    ]
    return written


def export_false_positive_figures() -> list[Path]:
    payloads = load_false_positive_payloads()
    stats = class_stats(payloads)
    intersections, membership = removed_overlap_data(payloads)
    written = [
        write_csv(workflow_counts(payloads, include_candidate_rows=True), "top5_top10_false_positive/01_workflow_counts.csv"),
        write_csv(removed_fraction(payloads), "top5_top10_false_positive/02_removed_fraction.csv"),
        write_csv(class_top_data(stats, 25), "top5_top10_false_positive/03_class_removed_top_data.csv"),
        write_csv(class_bubble_data(stats, 25), "top5_top10_false_positive/04_class_bubble_data.csv"),
        write_csv(score_distribution_data(payloads), "top5_top10_false_positive/05_score_distribution_data.csv"),
        write_csv(intersections, "top5_top10_false_positive/06_removed_overlap_intersections.csv"),
        write_csv(membership, "top5_top10_false_positive/06_removed_overlap_membership.csv"),
        write_csv(removed_fraction(payloads), "top5_top10_false_positive/summary_counts.csv"),
        write_csv(stats, "top5_top10_false_positive/class_removal_stats.csv"),
        write_csv(pd.concat([payload["annotations"] for payload in payloads.values()], ignore_index=True), "top5_top10_false_positive/annotation_removal_stats.csv"),
    ]
    return written


def main() -> None:
    written = []
    written.extend(export_rt_iup_tables())
    written.extend(export_topn_comparison())
    written.extend(export_false_positive_figures())
    for path in written:
        print(path)
    print(f"csv_count={len(written)}")


if __name__ == "__main__":
    main()
