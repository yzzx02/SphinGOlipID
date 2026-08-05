"""Build the canonical TOP3 RT/IUP input table from the fixed MS2 table.

The output is a CSV so downstream scientific calculations do not depend on an
intermediate Excel workbook.  Existing rows can inherit stable metadata from a
previous final workbook; newly added HILIC-fraction rows are derived from the
same fixed-table columns and matched to the archived MS1 feature clusters.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


CHAIN_RE = re.compile(r"(?P<base>[mdt])?(?P<c>\d{1,3}):(?P<u>\d{1,2})", re.IGNORECASE)
WINDOW_END = {
    "0-5": 5.0,
    "5-15": 15.0,
    "15-17.5": 17.5,
    "17.5-20.5": 20.5,
    "20.5-22": 22.0,
    "22-30": 30.0,
}


def text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def parse_detailed_chain(name: object) -> tuple[str, float | None, float | None]:
    """Return plot group, variable-chain carbon and unsaturation."""

    value = text(name)
    matches = list(CHAIN_RE.finditer(value))
    if not matches:
        return "", None, None
    variable = matches[-1]
    replacement = f"{variable.group('base') or ''}x:y"
    group = value[: variable.start()] + replacement + value[variable.end() :]
    return group, float(variable.group("c")), float(variable.group("u"))


def choose_cluster(
    mz: float,
    rt: float,
    clusters: pd.DataFrame,
) -> pd.Series | None:
    ppm = (clusters["feature_mz"] - mz).abs() / mz * 1e6
    rt_delta = (clusters["feature_rt"] - rt).abs()
    candidates = clusters.loc[(ppm <= 10.0) & (rt_delta <= 0.20)].copy()
    if candidates.empty:
        return None
    candidates["__ppm"] = ppm.loc[candidates.index]
    candidates["__rt"] = rt_delta.loc[candidates.index]
    return candidates.sort_values(["__rt", "__ppm", "ms1_feature_cluster_id"], kind="mergesort").iloc[0]


def previous_lookup(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    previous = pd.read_excel(path, sheet_name="详细链IUP_全部候选")
    keys = ["特征ID", "注释"]
    return previous.sort_values("合并后行ID", kind="mergesort").drop_duplicates(keys, keep="first").set_index(keys)


def build_table(fixed: pd.DataFrame, clusters: pd.DataFrame, previous: pd.DataFrame) -> pd.DataFrame:
    numeric_cluster_columns = ["feature_mz", "feature_rt", "normalized_rt", "abund_max"]
    for column in numeric_cluster_columns:
        clusters[column] = pd.to_numeric(clusters[column], errors="coerce")
    clusters = clusters.dropna(subset=["feature_mz", "feature_rt", "normalized_rt"]).copy()

    records: list[dict[str, object]] = []
    used_ids: set[str] = set()
    for row_number, row in enumerate(fixed.itertuples(index=False), start=1):
        source = row._asdict()
        lipid_name = text(source.get("lipid_name") or source.get("repr_注释"))
        feature_id = text(source.get("feature_id") or source.get("repr_feature_id"))
        window = text(source.get("source_sheets") or source.get("repr_source_sheet"))
        if ";" in window:
            window = window.split(";", 1)[0].strip()
        plot_group, x_carbon, curve_unsaturation = parse_detailed_chain(lipid_name)
        feature_mz = pd.to_numeric(pd.Series([source.get("feature_mz")]), errors="coerce").iloc[0]
        feature_rt = pd.to_numeric(pd.Series([source.get("ms1_rt")]), errors="coerce").iloc[0]
        cluster = None
        if pd.notna(feature_mz) and pd.notna(feature_rt):
            cluster = choose_cluster(float(feature_mz), float(feature_rt), clusters)

        previous_row: pd.Series | None = None
        if not previous.empty and (feature_id, lipid_name) in previous.index:
            candidate = previous.loc[(feature_id, lipid_name)]
            previous_row = candidate.iloc[0] if isinstance(candidate, pd.DataFrame) else candidate

        old_id = text(previous_row.get("合并后行ID")) if previous_row is not None else ""
        if old_id and old_id not in used_ids:
            merged_id = old_id
        else:
            merged_id = f"MR6_{row_number:07d}"
        used_ids.add(merged_id)

        if previous_row is not None:
            match_flag = text(previous_row.get("匹配到特征"))
            cluster_id = text(previous_row.get("特征簇ID"))
            out_feature_mz = previous_row.get("特征mz")
            out_feature_rt = previous_row.get("特征RT")
            feature_rt_delta = previous_row.get("特征RT差")
            feature_mz_ppm = previous_row.get("特征mz误差ppm")
            feature_abundance = previous_row.get("特征丰度")
            feature_source = previous_row.get("特征来源时间目录")
            normalized_rt = previous_row.get("归一化保留时间")
        elif cluster is not None:
            match_flag = "是"
            cluster_id = text(cluster.get("ms1_feature_cluster_id"))
            out_feature_mz = cluster.get("feature_mz")
            out_feature_rt = cluster.get("feature_rt")
            feature_rt_delta = abs(float(out_feature_rt) - float(feature_rt))
            feature_mz_ppm = abs(float(out_feature_mz) - float(feature_mz)) / float(feature_mz) * 1e6
            feature_abundance = cluster.get("abund_max")
            feature_source = text(cluster.get("best_source_file")).replace("-gradient1.csv", "").strip("() ").upper()
            normalized_rt = cluster.get("normalized_rt")
        else:
            match_flag = "否"
            cluster_id = ""
            out_feature_mz = np.nan
            out_feature_rt = np.nan
            feature_rt_delta = np.nan
            feature_mz_ppm = np.nan
            feature_abundance = np.nan
            feature_source = ""
            normalized_rt = float(feature_rt) - WINDOW_END.get(window, np.nan) if pd.notna(feature_rt) else np.nan

        record = {
            "合并后行ID": merged_id,
            "数据集": "top3",
            "时间窗口": window,
            "细类": plot_group,
            "注释": lipid_name,
            "母离子": source.get("repr_母离子"),
            "候选排名": source.get("best_rank"),
            "匹配度分数": source.get("best_ms2_match_score"),
            "总分数": source.get("repr_总分数"),
            "丰度": source.get("repr_相对强度"),
            "x碳数": x_carbon,
            "曲线不饱和度": curve_unsaturation,
            "归一化保留时间": normalized_rt,
            "特征簇ID": cluster_id,
            "匹配到特征": match_flag,
            "特征mz": out_feature_mz,
            "特征RT": out_feature_rt,
            "特征RT差": feature_rt_delta,
            "特征mz误差ppm": feature_mz_ppm,
            "特征丰度": feature_abundance,
            "特征来源时间目录": feature_source,
            "二级谱图RT均值": source.get("ms2_rt_mean"),
            "二级谱图RT列表": source.get("ms2_rt_list"),
            "合并谱图数": source.get("unique_ms2_spectrum_count"),
            "源行ID": source.get("repr_identification_id"),
            "特征ID": feature_id,
            "MS2扫描ID": source.get("ms2_scan_ids"),
        }
        records.append(record)
    output = pd.DataFrame(records)
    output["合并后行ID"] = output["合并后行ID"].astype(str)
    if output["合并后行ID"].duplicated().any():
        raise RuntimeError("合并后行ID is not unique")
    return output


def validate_overlap(output: pd.DataFrame, previous: pd.DataFrame) -> None:
    if previous.empty:
        return
    old = previous.reset_index()
    overlap = old.merge(output, on=["特征ID", "注释"], suffixes=("_old", "_new"), how="inner")
    print(f"Previous score-filtered rows: {len(old):,}; overlap: {len(overlap):,}")
    for column in ("细类", "x碳数", "曲线不饱和度", "归一化保留时间"):
        left = overlap[f"{column}_old"]
        right = overlap[f"{column}_new"]
        if column == "细类":
            matches = left.fillna("").astype(str).eq(right.fillna("").astype(str))
        else:
            matches = np.isclose(pd.to_numeric(left, errors="coerce"), pd.to_numeric(right, errors="coerce"), equal_nan=True)
        print(f"{column}: {int(matches.sum()):,}/{len(matches):,} exact")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-csv", type=Path, required=True)
    parser.add_argument("--ms1-clusters", type=Path, required=True)
    parser.add_argument("--previous-final", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    fixed = pd.read_csv(args.fixed_csv, low_memory=False)
    clusters = pd.read_csv(args.ms1_clusters, low_memory=False)
    previous = previous_lookup(args.previous_final)
    output = build_table(fixed, clusters, previous)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Input fixed rows: {len(fixed):,}")
    print(f"RT input rows: {len(output):,}")
    print(f"Parsed detailed-chain rows: {output['细类'].astype(str).str.len().gt(0).sum():,}")
    print(f"Matched MS1 rows: {output['匹配到特征'].eq('是').sum():,}")
    print(f"Output: {args.output}")
    validate_overlap(output, previous)


if __name__ == "__main__":
    main()
