"""Regenerate filtered RT/IUP plots from a top5/top10 SphinGOlipID workbook."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_rt_iup, write_rt_iup_plots


DEFAULT_INPUT = Path(
    "E:/李奕晓数据备份/20240422-Agilent/原始数据/SphinGOlipID_最终鉴定结果_链信息去重.xlsx"
)
DEFAULT_OUTPUT_DIR = Path("E:/李奕晓数据备份/20240422-Agilent/原始数据")

DATASETS = {
    "top5": "top5最终鉴定",
    "top10": "top10最终鉴定",
}

OLD_RT_IUP_COLUMNS = {
    "合并后预测保留时间",
    "合并后保留时间偏差",
    "合并后保留时间绝对偏差",
    "合并后是否作图",
    "合并后IUP说明",
    "预测保留时间",
    "保留时间偏差",
    "保留时间绝对偏差",
    "同曲线0.2min内点",
    "点数不足保留",
    "点数不足过滤原因",
    "拟合点IUP支持数",
    "拟合点IUP判断数",
    "IUP说明",
    "是否作图",
    "最终保留",
    "保留原因",
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
}

REQUIRED_COLUMNS = ["细类", "曲线不饱和度", "x碳数", "归一化保留时间"]
PAPER_B_PLOT_STYLE = {
    "figure_size": (4.72, 4.72),
    "point_size": 26,
    "two_point_line_width": 1.2,
    "fit_line_width": 1.5,
    "legend_font_size": 8,
    "axis_label_size": 10,
    "tick_label_size": 9,
}


def clean_input_table(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.drop(columns=[col for col in OLD_RT_IUP_COLUMNS if col in df.columns], errors="ignore").copy()
    missing = [col for col in REQUIRED_COLUMNS if col not in cleaned.columns]
    if missing:
        raise KeyError(f"Missing required column(s): {missing}")
    for col in ["曲线不饱和度", "x碳数", "归一化保留时间"]:
        cleaned[col] = pd.to_numeric(cleaned[col], errors="coerce")
    return cleaned.dropna(subset=REQUIRED_COLUMNS).reset_index(drop=True)


def excel_safe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "参数" in out.columns:
        out["参数"] = out["参数"].map(
            lambda value: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, tuple)) else value
        )
    return out


def write_result_workbook(results: dict[str, object], output_path: Path) -> None:
    summary_rows = []
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for dataset, payload in results.items():
            result = payload["result"]
            summary = {"数据集": dataset, **result.stats, "图数量": payload["plot_count"]}
            summary_rows.append(summary)
            excel_safe(result.rows).to_excel(writer, sheet_name=f"{dataset}重新过滤结果", index=False)
            excel_safe(result.lines).to_excel(writer, sheet_name=f"{dataset}重新拟合曲线", index=False)
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name="summary", index=False)


def clear_old_plots(plot_dir: Path) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    for path in plot_dir.glob("*_all.png"):
        path.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-xlsx", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--r2-threshold", type=float, default=0.99)
    parser.add_argument("--rt-window-min", type=float, default=0.2)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = RTIUPConfig(r2_threshold=args.r2_threshold, rt_window_min=args.rt_window_min, **PAPER_B_PLOT_STYLE)
    results: dict[str, object] = {}
    for dataset, sheet_name in DATASETS.items():
        raw = pd.read_excel(args.input_xlsx, sheet_name=sheet_name)
        prepared = clean_input_table(raw)
        result = fit_rt_iup(prepared, config)
        plot_dir = args.output_dir / f"{dataset}_RT_IUP_filtered_plots"
        clear_old_plots(plot_dir)
        plot_count = write_rt_iup_plots(result.plot_rows, result.lines, plot_dir, config)
        results[dataset] = {"result": result, "plot_count": plot_count, "plot_dir": plot_dir}
        print(f"{dataset}: input={len(prepared)}, kept={len(result.rows)}, plotted={len(result.plot_rows)}, plots={plot_count}")
        print(f"{dataset}: plot_dir={plot_dir}")

    workbook_path = args.output_dir / "SphinGOlipID_RT_IUP_refiltered_top5_top10.xlsx"
    write_result_workbook(results, workbook_path)
    print(f"summary_workbook={workbook_path}")


if __name__ == "__main__":
    main()
