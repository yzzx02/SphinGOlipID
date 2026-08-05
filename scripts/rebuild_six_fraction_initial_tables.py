"""Rebuild six-fraction initial tables with the corrected fifth fraction.

The first four fractions and the sixth fraction are read from the previously
validated rank-10 CSV export.  The fifth fraction (20.5--22 min) is replaced by
the independently re-identified rank-10 CSV.  Feature and identification IDs
for the replacement fraction are namespaced to prevent collisions.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from prepare_missing_fraction_workbooks import WINDOWS, deduplicate


def numeric_id_max(values: pd.Series) -> int:
    numbers = pd.to_numeric(
        values.astype(str).str.extract(r"(\d+)", expand=False), errors="coerce"
    )
    return int(numbers.max()) if numbers.notna().any() else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validated-source-dir", required=True, type=Path)
    parser.add_argument("--corrected-fifth-rank10", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    by_window: dict[str, pd.DataFrame] = {}
    base_columns: list[str] | None = None
    for window in WINDOWS:
        path = args.validated_source_dir / f"source_{window}.csv"
        frame = pd.read_csv(path, low_memory=False)
        frame.columns = [str(column).lstrip("\ufeff") for column in frame.columns]
        if base_columns is None:
            base_columns = frame.columns.tolist()
        elif frame.columns.tolist() != base_columns:
            raise RuntimeError(f"Column mismatch in {path}")
        by_window[window] = frame
    assert base_columns is not None

    fifth = pd.read_csv(args.corrected_fifth_rank10, low_memory=False)
    fifth.columns = [str(column).lstrip("\ufeff") for column in fifth.columns]
    if "source_sheet" in fifth.columns:
        if not fifth["source_sheet"].astype(str).eq("20.5-22").all():
            raise RuntimeError("Corrected fifth CSV contains a non-20.5-22 source_sheet")
        fifth = fifth.drop(columns=["source_sheet"])
    missing = [column for column in base_columns if column not in fifth.columns]
    if missing:
        raise RuntimeError(f"Corrected fifth CSV is missing columns: {missing}")
    fifth = fifth[base_columns].copy()

    feature_ids = sorted(fifth["feature_id"].astype(str).unique())
    feature_map = {old: f"F205{index:06d}" for index, old in enumerate(feature_ids, start=1)}
    fifth["feature_id"] = fifth["feature_id"].astype(str).map(feature_map)

    existing = pd.concat(
        [by_window[window] for window in WINDOWS if window != "20.5-22"],
        ignore_index=True,
    )
    by_window["20.5-22"] = fifth

    # Legacy sheets contain repeated identification_id values across fractions.
    # Reissue this row-level identifier globally so the rebuilt six-fraction
    # source table has an unambiguous primary key.
    next_identification = 1
    for window in WINDOWS:
        row_count = len(by_window[window])
        by_window[window]["identification_id"] = [
            f"ID{value:08d}"
            for value in range(next_identification, next_identification + row_count)
        ]
        next_identification += row_count

    for window in WINDOWS:
        by_window[window].to_csv(
            args.output_dir / f"source_{window}.csv", index=False, encoding="utf-8-sig"
        )

    combined_frames: list[pd.DataFrame] = []
    for window in WINDOWS:
        frame = by_window[window].copy()
        frame.insert(0, "source_sheet", window)
        combined_frames.append(frame)
    combined = pd.concat(combined_frames, ignore_index=True, sort=False)

    counts: list[dict[str, int]] = []
    dedup_frames: dict[int, pd.DataFrame] = {}
    numeric_rank = pd.to_numeric(combined["rank"], errors="coerce")
    for top_n in (3, 5, 10):
        raw = combined[numeric_rank <= top_n].copy().reset_index(drop=True)
        raw.to_csv(
            args.output_dir / f"TOP{top_n}_raw_rows.csv", index=False, encoding="utf-8-sig"
        )
        dedup = deduplicate(raw, top_n)
        dedup.to_csv(
            args.output_dir / f"TOP{top_n}_feature_lipid.csv", index=False, encoding="utf-8-sig"
        )
        dedup_frames[top_n] = dedup
        counts.append(
            {
                "top_n": top_n,
                "raw_candidate_rows": len(raw),
                "feature_lipid_rows": len(dedup),
                "unique_feature_ids": raw["feature_id"].nunique(),
            }
        )

    pd.DataFrame(counts).to_csv(
        args.output_dir / "counts.csv", index=False, encoding="utf-8-sig"
    )
    examples = dedup_frames[3].sort_values(
        ["raw_candidate_row_count", "feature_id", "lipid_name"],
        ascending=[False, True, True],
    ).head(10)
    pd.DataFrame(
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
    ).to_csv(args.output_dir / "dedup_examples.csv", index=False, encoding="utf-8-sig")

    fraction_rows = []
    for window in WINDOWS:
        frame = by_window[window]
        rank = pd.to_numeric(frame["rank"], errors="coerce")
        fraction_rows.append(
            {
                "fraction": window,
                "rank10_candidate_rows": len(frame),
                "top3_candidate_rows": int((rank <= 3).sum()),
                "top3_unique_features": int(frame.loc[rank <= 3, "feature_id"].nunique()),
                "top3_feature_lipid_rows": int(
                    frame.loc[rank <= 3, ["feature_id", "lipid_name"]].drop_duplicates().shape[0]
                ),
            }
        )
    pd.DataFrame(fraction_rows).to_csv(
        args.output_dir / "fraction_counts.csv", index=False, encoding="utf-8-sig"
    )

    if combined["identification_id"].astype(str).duplicated().any():
        raise RuntimeError("identification_id collision after rebuild")
    fifth_ids = set(fifth["feature_id"].astype(str))
    other_ids = set(existing["feature_id"].astype(str))
    if fifth_ids & other_ids:
        raise RuntimeError("feature_id collision after fifth-fraction namespacing")

    print(pd.DataFrame(counts).to_string(index=False))
    print(pd.DataFrame(fraction_rows).to_string(index=False))
    print(f"Corrected fifth feature IDs: {len(feature_map):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
