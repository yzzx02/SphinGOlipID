from pathlib import Path

import pandas as pd


def test_fit_rt_iup_and_plot(tmp_path: Path):
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_rt_iup, make_rt_iup_plot

    rows = []
    for unsat, intercept in [(0, 1.2), (1, 0.55), (4, -0.65)]:
        for x in [16, 18, 20, 22, 24]:
            rows.append(
                {
                    "细类": "SM(d18:0/x:y)",
                    "曲线不饱和度": float(unsat),
                    "x碳数": float(x),
                    "归一化保留时间": 0.4 * x + intercept,
                    "总分数": 40,
                    "匹配度分数": 0.5,
                    "丰度": 1.0,
                }
            )
    # Two-point series should be retained for plotting by IUP and drawn as a
    # visual connector, but it is not treated as a fitted curve.
    for x, y in [(20, 7.35), (22, 8.15)]:
        rows.append(
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": float(x),
                "归一化保留时间": y,
                "总分数": 40,
                "匹配度分数": 0.5,
                "丰度": 1.0,
            }
        )

    result = fit_rt_iup(pd.DataFrame(rows), RTIUPConfig(r2_threshold=0.98))

    assert result.stats["有效拟合曲线数"] == 3
    assert set(result.lines["不饱和度"]) == {0.0, 1.0, 4.0}
    assert result.plot_rows[result.plot_rows["曲线不饱和度"].eq(2.0)].shape[0] == 2

    plot_path = tmp_path / "sm.png"
    make_rt_iup_plot("SM(d18:0/x:y)", result.plot_rows, result.lines, plot_path)
    assert plot_path.exists()
    assert plot_path.stat().st_size > 0
