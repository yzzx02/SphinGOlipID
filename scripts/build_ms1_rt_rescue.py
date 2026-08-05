"""Rebuild RT-supported MS1-only sphingolipid rescue candidates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sphingolipid_toolkit.rt_iup import predict_values  # noqa: E402


RT_TOLERANCE_MIN = 0.20
MZ_TOLERANCE_PPM = 10.0


def no_ms2_feature_ids(ms1: pd.DataFrame, top3: pd.DataFrame) -> set[str]:
    features = ms1[["ms1_feature_cluster_id", "feature_mz", "normalized_rt"]].drop_duplicates("ms1_feature_cluster_id")
    ms2 = top3[top3["匹配到特征"].astype(str).str.strip().eq("是")][["特征mz", "归一化保留时间"]].dropna().drop_duplicates()
    ms2_mz = pd.to_numeric(ms2["特征mz"], errors="coerce").to_numpy(dtype=float)
    ms2_rt = pd.to_numeric(ms2["归一化保留时间"], errors="coerce").to_numpy(dtype=float)
    no_ms2: set[str] = set()
    for row in features.itertuples(index=False):
        mz = float(row.feature_mz)
        rt = float(row.normalized_rt)
        ppm = np.abs(ms2_mz - mz) / mz * 1e6
        matched = bool(np.any((ppm <= MZ_TOLERANCE_PPM) & (np.abs(ms2_rt - rt) <= RT_TOLERANCE_MIN)))
        if not matched:
            no_ms2.add(str(row.ms1_feature_cluster_id))
    return no_ms2


def line_lookup(path: Path) -> dict[tuple[str, float], pd.Series]:
    lines = pd.read_csv(path)
    lookup: dict[tuple[str, float], pd.Series] = {}
    for _, line in lines.iterrows():
        unsat = pd.to_numeric(pd.Series([line.get("不饱和度")]), errors="coerce").iloc[0]
        if pd.isna(unsat):
            continue
        lookup[(str(line.get("细类")), float(unsat))] = line
    return lookup


def line_support(line: pd.Series | None, carbon: float, rt: float) -> tuple[bool, float | None]:
    if line is None:
        return False, None
    x_min = pd.to_numeric(pd.Series([line.get("x最小")]), errors="coerce").iloc[0]
    x_max = pd.to_numeric(pd.Series([line.get("x最大")]), errors="coerce").iloc[0]
    if pd.isna(x_min) or pd.isna(x_max) or not (float(x_min) <= carbon <= float(x_max)):
        return False, None
    # IUP stores the bounded whole-series display/order correction in 参数 and
    # keeps the data-derived fit in IUP原始参数. An MS1-only candidate belongs
    # to the same series and would receive the same constant correction, so its
    # point-to-curve residual must be evaluated against the original fit (the
    # shift cancels on both sides) rather than being penalized by the display
    # offset a second time.
    original_params = line.get("IUP原始参数")
    params = original_params if isinstance(original_params, str) and original_params.strip() not in {"", "[]", "nan"} else line.get("参数")
    predicted = predict_values(line.get("拟合类型"), params, carbon)
    if predicted is None:
        return False, None
    residual = abs(rt - float(predicted))
    return residual <= RT_TOLERANCE_MIN, residual


def evaluate_method(
    candidates: pd.DataFrame,
    total_lines: dict[tuple[str, float], pd.Series],
    lcb_lines: dict[tuple[str, float], pd.Series],
    method: str,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for _, row in candidates.iterrows():
        base_class = str(row.get("base_class") or "").strip()
        lcb_type = str(row.get("lcb_type") or "").strip().lower()
        carbon = pd.to_numeric(pd.Series([row.get("总碳数")]), errors="coerce").iloc[0]
        unsat = pd.to_numeric(pd.Series([row.get("总不饱和度")]), errors="coerce").iloc[0]
        rt = pd.to_numeric(pd.Series([row.get("normalized_rt")]), errors="coerce").iloc[0]
        if not base_class or lcb_type not in {"d", "m", "t"} or pd.isna(carbon) or pd.isna(unsat) or pd.isna(rt):
            continue
        total_ok, total_residual = line_support(total_lines.get((base_class, float(unsat))), float(carbon), float(rt))
        lcb_group = f"{base_class}({lcb_type})"
        lcb_ok, lcb_residual = line_support(lcb_lines.get((lcb_group, float(unsat))), float(carbon), float(rt))
        if not (total_ok or lcb_ok):
            continue
        record = row.to_dict()
        record["RT规则"] = method
        record["总碳数RT支持"] = total_ok
        record["长链碱系列RT支持"] = lcb_ok
        record["总碳数RT绝对偏差"] = total_residual
        record["长链碱系列RT绝对偏差"] = lcb_residual
        record["MS2匹配状态"] = "无匹配MS2，MS1-only"
        record["MS1捞回判定"] = "RT支持，捞回"
        records.append(record)
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ms1-candidates", type=Path, required=True)
    parser.add_argument("--top3", type=Path, required=True)
    parser.add_argument("--plots-root", type=Path, required=True)
    parser.add_argument("--output-detail", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    args = parser.parse_args()

    ms1 = pd.read_csv(args.ms1_candidates, low_memory=False)
    if args.top3.suffix.lower() == ".csv":
        top3 = pd.read_csv(args.top3, low_memory=False)
    else:
        top3 = pd.read_excel(args.top3, sheet_name=0)
    for column in ("feature_mz", "normalized_rt", "总碳数", "总不饱和度"):
        ms1[column] = pd.to_numeric(ms1[column], errors="coerce")
    ms1 = ms1.dropna(subset=["ms1_feature_cluster_id", "feature_mz", "normalized_rt", "总碳数", "总不饱和度"])
    no_ms2 = no_ms2_feature_ids(ms1, top3)
    candidates = ms1[ms1["ms1_feature_cluster_id"].astype(str).isin(no_ms2)].copy()

    outputs: dict[str, pd.DataFrame] = {}
    for method in ("ECN", "IUP"):
        total_lines = line_lookup(args.plots_root / "总碳数总不饱和度分类" / f"{method}_表" / "有效拟合曲线.csv")
        lcb_lines = line_lookup(args.plots_root / "长链碱系列分类" / f"{method}_表" / "有效拟合曲线.csv")
        outputs[method] = evaluate_method(candidates, total_lines, lcb_lines, method)

    key_columns = ["ms1_feature_cluster_id", "candidate_total_name"]
    def supported_keys(frame: pd.DataFrame, support_column: str | None = None) -> set[tuple[str, str]]:
        if frame.empty:
            return set()
        selected = frame if support_column is None else frame[frame[support_column].fillna(False).astype(bool)]
        return set(map(tuple, selected[key_columns].astype(str).to_numpy()))

    ecn_keys = supported_keys(outputs["ECN"])
    iup_keys = supported_keys(outputs["IUP"])
    key_sets = {
        "总碳数总不饱和度分类": {
            "ECN": supported_keys(outputs["ECN"], "总碳数RT支持"),
            "IUP": supported_keys(outputs["IUP"], "总碳数RT支持"),
        },
        "长链碱系列分类": {
            "ECN": supported_keys(outputs["ECN"], "长链碱系列RT支持"),
            "IUP": supported_keys(outputs["IUP"], "长链碱系列RT支持"),
        },
    }
    raw_combined = pd.concat(outputs.values(), ignore_index=True).drop_duplicates(key_columns + ["RT规则"])
    consolidated_records: list[dict[str, object]] = []
    for _, group in raw_combined.groupby(key_columns, sort=False):
        record = group.iloc[0].to_dict()
        pair_key = tuple(group.iloc[0][key_columns].astype(str))
        record["RT规则"] = "ECN且IUP" if pair_key in ecn_keys & iup_keys else ("ECN" if pair_key in ecn_keys else "IUP")
        record["ECN支持"] = pair_key in ecn_keys
        record["IUP支持"] = pair_key in iup_keys
        ecn_group = group[group["RT规则"].eq("ECN")]
        iup_group = group[group["RT规则"].eq("IUP")]
        record["ECN总碳数RT支持"] = bool(ecn_group["总碳数RT支持"].fillna(False).astype(bool).any())
        record["ECN长链碱系列RT支持"] = bool(ecn_group["长链碱系列RT支持"].fillna(False).astype(bool).any())
        record["IUP总碳数RT支持"] = bool(iup_group["总碳数RT支持"].fillna(False).astype(bool).any())
        record["IUP长链碱系列RT支持"] = bool(iup_group["长链碱系列RT支持"].fillna(False).astype(bool).any())
        record["总碳数RT支持"] = bool(group["总碳数RT支持"].fillna(False).astype(bool).any())
        record["长链碱系列RT支持"] = bool(group["长链碱系列RT支持"].fillna(False).astype(bool).any())
        record["总碳数RT绝对偏差"] = pd.to_numeric(group.loc[group["总碳数RT支持"].fillna(False).astype(bool), "总碳数RT绝对偏差"], errors="coerce").min()
        record["长链碱系列RT绝对偏差"] = pd.to_numeric(group.loc[group["长链碱系列RT支持"].fillna(False).astype(bool), "长链碱系列RT绝对偏差"], errors="coerce").min()
        consolidated_records.append(record)
    combined = pd.DataFrame(consolidated_records)

    summary_rows: list[dict[str, object]] = []
    for classification, method_keys in key_sets.items():
        ecn_class_keys = method_keys["ECN"]
        iup_class_keys = method_keys["IUP"]
        for label, keys in (
            ("ECN", ecn_class_keys),
            ("IUP", iup_class_keys),
            ("ECN且IUP", ecn_class_keys & iup_class_keys),
            ("ECN或IUP", ecn_class_keys | iup_class_keys),
        ):
            subset = combined[combined[key_columns].astype(str).apply(tuple, axis=1).isin(keys)].drop_duplicates(key_columns)
            summary_rows.append(
                {
                    "分类模式": classification,
                    "RT口径": label,
                    "MS1-only候选配对": int(len(keys)),
                    "MS1-only feature": int(subset["ms1_feature_cluster_id"].nunique()),
                    "MS1捞回脂质名称": int(subset["candidate_total_name"].nunique()),
                }
            )

    args.output_detail.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output_detail, index=False, encoding="utf-8-sig")
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.output_summary, index=False, encoding="utf-8-sig")
    print(f"Original MS1 candidates: {len(ms1)}")
    print(f"MS1-only features: {len(no_ms2)}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
