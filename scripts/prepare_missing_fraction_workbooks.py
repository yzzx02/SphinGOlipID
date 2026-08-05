"""Prepare CSV datasets for rebuilt initial SphinGOlipID workbooks.

Excel files are authored separately with the bundled spreadsheet artifact
runtime.  This module performs only the scientific tabular transformations.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


WINDOWS = ["0-5", "5-15", "15-17.5", "17.5-20.5", "20.5-22", "22-30"]
REPRESENTATIVE_RULE = (
    "lowest rank, then highest ms2_match_score, matched_fragment_count, "
    "matched_intensity_sum, then lowest fragment_error_ppm_mean"
)


def unique_text(values: pd.Series, separator: str = "; ") -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value)
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return separator.join(ordered)


def fmt_float(value: object, digits: int = 6) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return ""


def select_representative(group: pd.DataFrame) -> pd.Series:
    ranked = group.copy()
    ranked["_rank"] = pd.to_numeric(ranked["rank"], errors="coerce").fillna(np.inf)
    ranked["_score"] = pd.to_numeric(ranked["ms2_match_score"], errors="coerce").fillna(-np.inf)
    ranked["_fragments"] = pd.to_numeric(ranked["matched_fragment_count"], errors="coerce").fillna(-np.inf)
    ranked["_intensity"] = pd.to_numeric(ranked["matched_intensity_sum"], errors="coerce").fillna(-np.inf)
    ranked["_ppm"] = pd.to_numeric(ranked["fragment_error_ppm_mean"], errors="coerce").fillna(np.inf)
    return ranked.sort_values(
        ["_rank", "_score", "_fragments", "_intensity", "_ppm"],
        ascending=[True, False, False, False, True],
        kind="stable",
    ).iloc[0]


def deduplicate(raw: pd.DataFrame, top_n: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (feature_id, lipid_name), group in raw.groupby(["feature_id", "lipid_name"], sort=True, dropna=False):
        representative = select_representative(group)
        ms2_rt = pd.to_numeric(group["ms2_rt"], errors="coerce")
        spectrum_keys = group["ms2_file_name"].astype(str) + "|" + group["ms2_scan_id"].astype(str)
        rank_values = unique_text(group["rank"], separator="; ")
        rt_list = "; ".join(fmt_float(value) for value in ms2_rt if pd.notna(value))
        file_scan_rt = "; ".join(
            f"{file_name}|{scan_id}|rt={fmt_float(rt)}"
            for file_name, scan_id, rt in zip(group["ms2_file_name"], group["ms2_scan_id"], ms2_rt)
            if pd.notna(rt)
        )
        item: dict[str, object] = {
            "topN": top_n,
            "feature_id": feature_id,
            "lipid_name": lipid_name,
            "raw_candidate_row_count": len(group),
            "unique_ms2_spectrum_count": spectrum_keys.nunique(),
            "rank_values": rank_values,
            "best_rank": representative.get("rank"),
            "best_ms2_match_score": representative.get("ms2_match_score"),
            "best_matched_fragment_count": representative.get("matched_fragment_count"),
            "best_matched_intensity_sum": representative.get("matched_intensity_sum"),
            "ms1_rt": representative.get("ms1_rt"),
            "feature_mz": representative.get("feature_mz"),
            "feature_mz_error_ppm": representative.get("feature_mz_error_ppm"),
            "ms1_library_mz_error_ppm": representative.get("ms1_library_mz_error_ppm"),
            "ms2_rt_min": ms2_rt.min(),
            "ms2_rt_mean": ms2_rt.mean(),
            "ms2_rt_max": ms2_rt.max(),
            "ms2_rt_list": rt_list,
            "ms2_file_scan_rt_list": file_scan_rt,
            "collision_energies": unique_text(group["collision_energy"]),
            "source_sheets": unique_text(group["source_sheet"]),
            "targetlist_groups": unique_text(group["targetlist_groups"]),
            "ms2_file_names": unique_text(group["ms2_file_name"]),
            "ms2_scan_ids": unique_text(group["ms2_scan_id"]),
            "representative_selection_rule": REPRESENTATIVE_RULE,
        }
        for column in raw.columns:
            item[f"repr_{column}"] = representative.get(column)
        rows.append(item)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-rank10", required=True, type=Path)
    parser.add_argument("--new-rank10-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    excel = pd.ExcelFile(args.source_rank10)
    by_window: dict[str, pd.DataFrame] = {}
    base_columns: list[str] | None = None
    for window in WINDOWS[:-1]:
        frame = pd.read_excel(args.source_rank10, sheet_name=window)
        if base_columns is None:
            base_columns = [str(column).lstrip("\ufeff") for column in frame.columns]
        frame.columns = [str(column).lstrip("\ufeff") for column in frame.columns]
        frame.insert(0, "source_sheet", window)
        by_window[window] = frame
    assert base_columns is not None

    new = pd.read_csv(args.new_rank10_csv)
    new.columns = [str(column).lstrip("\ufeff") for column in new.columns]
    for column in base_columns:
        if column not in new.columns:
            new[column] = np.nan
    by_window["22-30"] = new[["source_sheet", *base_columns]].copy()

    combined = pd.concat([by_window[window] for window in WINDOWS], ignore_index=True, sort=False)
    old_mask = ~combined["source_sheet"].eq("22-30")
    old_numbers = pd.to_numeric(
        combined.loc[old_mask, "identification_id"].astype(str).str.extract(r"(\d+)", expand=False),
        errors="coerce",
    )
    next_id = int(old_numbers.max()) + 1
    new_count = int((~old_mask).sum())
    combined.loc[~old_mask, "identification_id"] = [
        f"ID{i:08d}" for i in range(next_id, next_id + new_count)
    ]

    for window in WINDOWS:
        source = combined[combined["source_sheet"].eq(window)][base_columns]
        source.to_csv(args.output_dir / f"source_{window}.csv", index=False, encoding="utf-8-sig")

    counts: list[dict[str, int]] = []
    dedup_frames: dict[int, pd.DataFrame] = {}
    for top_n in (3, 5, 10):
        rank = pd.to_numeric(combined["rank"], errors="coerce")
        raw = combined[rank <= top_n].copy().reset_index(drop=True)
        raw.to_csv(args.output_dir / f"TOP{top_n}_raw_rows.csv", index=False, encoding="utf-8-sig")
        dedup = deduplicate(raw, top_n)
        dedup.to_csv(args.output_dir / f"TOP{top_n}_feature_lipid.csv", index=False, encoding="utf-8-sig")
        dedup_frames[top_n] = dedup
        counts.append(
            {
                "top_n": top_n,
                "raw_candidate_rows": len(raw),
                "feature_lipid_rows": len(dedup),
            }
        )

    pd.DataFrame(counts).to_csv(args.output_dir / "counts.csv", index=False, encoding="utf-8-sig")
    examples = dedup_frames[3].sort_values(
        ["raw_candidate_row_count", "feature_id", "lipid_name"],
        ascending=[False, True, True],
    ).head(10)
    example_table = pd.DataFrame(
        {
            "topN": "TOP3",
            "example_no": range(1, len(examples) + 1),
            "feature_id": examples["feature_id"].values,
            "lipid_name": examples["lipid_name"].values,
            "raw_rows_before_dedup": examples["raw_candidate_row_count"].values,
            "unique_ms2_spectra": examples["unique_ms2_spectrum_count"].values,
            "rank_values": examples["rank_values"].values,
            "ms2_rt_list": examples["ms2_rt_list"].values,
            "ms2_file_scan_rt_list": examples["ms2_file_scan_rt_list"].values,
            "explanation": (
                "Repeated MS2 evidence for the same MS1 feature and same lipid candidate; "
                "dedup keeps one feature-level candidate while retaining all MS2 RT evidence."
            ),
        }
    )
    example_table.to_csv(args.output_dir / "dedup_examples.csv", index=False, encoding="utf-8-sig")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
