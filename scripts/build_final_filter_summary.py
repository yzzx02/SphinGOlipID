"""Calculate the final ECN/IUP filter-summary table on the full TOP3 universe."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from regenerate_top3_rt_ecn_iup_preview import (  # noqa: E402
    PLOT_STYLE,
    prepare_groupings,
)
from sphingolipid_toolkit.rt_iup import (  # noqa: E402
    RTIUPConfig,
    fit_ecn,
    fit_rt_iup,
)


ID_COLUMN = "合并后行ID"
FEATURE_COLUMN = "特征ID"
NAME_COLUMN = "注释"
SCORE_COLUMN = "匹配度分数"


def normalized_values(series: pd.Series) -> set[str]:
    return set(series.dropna().astype(str).str.strip().loc[lambda values: values.ne("")])


def summarize(label: str, universe: pd.DataFrame, ids: set[str]) -> dict[str, object]:
    selected = universe[universe[ID_COLUMN].astype(str).isin(ids)].drop_duplicates(ID_COLUMN)
    return {
        "过滤口径": label,
        "初始总候选": int(universe[ID_COLUMN].nunique()),
        "过滤后剩余": int(selected[ID_COLUMN].nunique()),
        "MS1 feature": int(selected[FEATURE_COLUMN].dropna().astype(str).str.strip().loc[lambda values: values.ne("")].nunique()),
    }


def rescue_counts_for_group(rescue: pd.DataFrame, grouping: str) -> dict[str, int | str]:
    if grouping == "详细链分类":
        return {"ECN": "不适用", "IUP": "不适用", "ECN且IUP": "不适用", "ECN或IUP": "不适用"}
    subset = rescue[rescue["分类模式"].eq(grouping)].set_index("RT口径")["MS1-only候选配对"].astype(int)
    return subset.to_dict()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--ms1-rescue-summary", type=Path, required=True)
    parser.add_argument(
        "--scored-results-root",
        type=Path,
        default=None,
        help="Optional already-generated score>=0.50 ECN/IUP result root.",
    )
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    if args.input.suffix.lower() == ".csv":
        raw = pd.read_csv(args.input, low_memory=False)
    else:
        raw = pd.read_excel(args.input, sheet_name=0)
    groupings = prepare_groupings(raw, -np.inf)
    config = replace(RTIUPConfig(), **PLOT_STYLE)
    rescue = pd.read_csv(args.ms1_rescue_summary)
    for grouping, universe in groupings.items():
        ecn = fit_ecn(universe, config)
        iup = fit_rt_iup(universe, config)
        ecn_ids = normalized_values(ecn.rows[ID_COLUMN])
        iup_ids = normalized_values(iup.rows[ID_COLUMN])

        if args.scored_results_root is not None:
            ecn_scored_rows = pd.read_csv(
                args.scored_results_root / grouping / "ECN_表" / "最终保留明细.csv",
                low_memory=False,
            )
            iup_scored_rows = pd.read_csv(
                args.scored_results_root / grouping / "IUP_表" / "最终保留明细.csv",
                low_memory=False,
            )
            ecn_scored_ids = normalized_values(ecn_scored_rows[ID_COLUMN])
            iup_scored_ids = normalized_values(iup_scored_rows[ID_COLUMN])
        else:
            scored_universe = universe.loc[
                pd.to_numeric(universe[SCORE_COLUMN], errors="coerce").ge(0.50)
            ].copy()
            ecn_scored = fit_ecn(scored_universe, config)
            iup_scored = fit_rt_iup(scored_universe, config)
            ecn_scored_ids = normalized_values(ecn_scored.rows[ID_COLUMN])
            iup_scored_ids = normalized_values(iup_scored.rows[ID_COLUMN])
        rows = [
            summarize("ECN RT筛选后保留", universe, ecn_ids),
            summarize("IUP RT筛选后保留", universe, iup_ids),
            summarize("ECN RT保留且 score≥0.50", universe, ecn_scored_ids),
            summarize("IUP RT保留且 score≥0.50", universe, iup_scored_ids),
            summarize("ECN、IUP同时保留且 score≥0.50", universe, ecn_scored_ids & iup_scored_ids),
            summarize("ECN或IUP任一保留且 score≥0.50", universe, ecn_scored_ids | iup_scored_ids),
        ]
        output = pd.DataFrame(rows)
        grouping_rescue = rescue_counts_for_group(rescue, grouping)
        rescue_by_row = {
            "ECN RT筛选后保留": grouping_rescue["ECN"],
            "IUP RT筛选后保留": grouping_rescue["IUP"],
            "ECN RT保留且 score≥0.50": grouping_rescue["ECN"],
            "IUP RT保留且 score≥0.50": grouping_rescue["IUP"],
            "ECN、IUP同时保留且 score≥0.50": grouping_rescue["ECN且IUP"],
            "ECN或IUP任一保留且 score≥0.50": grouping_rescue["ECN或IUP"],
        }
        output["MS1 RT捞回"] = output["过滤口径"].map(rescue_by_row)
        output_path = args.output_root / grouping / "过滤口径统计.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output.to_csv(output_path, index=False, encoding="utf-8-sig")
        print(f"\n[{grouping}]")
        print(output.to_string(index=False))


if __name__ == "__main__":
    main()
