"""Build top3 RT/IUP results from the chain-deduplicated top5 table."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_rt_iup, write_rt_iup_plots


BASE_DIR = Path("E:/李奕晓数据备份/20240422-Agilent/原始数据")
SOURCE_WORKBOOK = BASE_DIR / "SphinGOlipID_最终鉴定结果_链信息去重.xlsx"
OUTPUT_WORKBOOK = BASE_DIR / "SphinGOlipID_top3_归一化保留时间_IUP最终拟合结果.xlsx"
PLOT_DIR = BASE_DIR / "top3_RT_IUP_filtered_plots"

SOURCE_SHEET = "top5最终鉴定"
RANK_LIMIT = 3
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


def excel_safe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "参数" in out.columns:
        out["参数"] = out["参数"].map(
            lambda value: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, tuple)) else value
        )
    return out


def prepare_top3_source() -> pd.DataFrame:
    raw = pd.read_excel(SOURCE_WORKBOOK, sheet_name=SOURCE_SHEET)
    raw["候选排名"] = pd.to_numeric(raw["候选排名"], errors="coerce")
    top3 = raw[raw["候选排名"].le(RANK_LIMIT)].copy()
    top3 = top3.drop(columns=[col for col in OLD_RT_IUP_COLUMNS if col in top3.columns], errors="ignore")
    missing = [col for col in REQUIRED_COLUMNS if col not in top3.columns]
    if missing:
        raise KeyError(f"Missing required column(s): {missing}")
    for col in ["曲线不饱和度", "x碳数", "归一化保留时间"]:
        top3[col] = pd.to_numeric(top3[col], errors="coerce")
    if "数据集" in top3.columns:
        top3["数据集"] = "top3"
    if "最终鉴定ID" in top3.columns:
        top3["最终鉴定ID"] = [f"TOP3_ID{i:06d}" for i in range(1, len(top3) + 1)]
    return top3.dropna(subset=REQUIRED_COLUMNS).reset_index(drop=True)


def clear_old_plots() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    for path in PLOT_DIR.glob("*_all.png"):
        path.unlink()


def main() -> None:
    config = RTIUPConfig(r2_threshold=0.99, rt_window_min=0.2, **PAPER_B_PLOT_STYLE)
    source = prepare_top3_source()
    result = fit_rt_iup(source, config)
    clear_old_plots()
    plot_count = write_rt_iup_plots(result.plot_rows, result.lines, PLOT_DIR, config)

    summary = pd.DataFrame(
        [
            {
                "数据集": "top3",
                **result.stats,
                "图数量": plot_count,
                "过滤前唯一注释数": len(source),
            }
        ]
    )

    with pd.ExcelWriter(OUTPUT_WORKBOOK, engine="openpyxl") as writer:
        source.to_excel(writer, sheet_name="top3过滤前唯一注释", index=False)
        excel_safe(result.rows).to_excel(writer, sheet_name="最终拟合结果", index=False)
        excel_safe(result.lines).to_excel(writer, sheet_name="拟合曲线结果", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)

    print(f"top3: source={len(source)}, kept={len(result.rows)}, plotted={len(result.plot_rows)}, plots={plot_count}")
    print(f"workbook={OUTPUT_WORKBOOK}")
    print(f"plot_dir={PLOT_DIR}")


if __name__ == "__main__":
    main()
