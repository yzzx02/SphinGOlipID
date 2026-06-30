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


def test_decreasing_tail_guard_keeps_recovering_third_point():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, _has_clear_decreasing_tail

    config = RTIUPConfig()
    recovering = pd.DataFrame(
        {
            "x碳数": [20.0, 22.0, 24.0],
            "归一化保留时间": [8.0, 7.8, 7.92],
        }
    )
    falling_tail = pd.DataFrame(
        {
            "x碳数": [20.0, 22.0, 24.0],
            "归一化保留时间": [7.6, 8.0, 7.85],
        }
    )

    assert not _has_clear_decreasing_tail(recovering, config)
    assert _has_clear_decreasing_tail(falling_tail, config)


def test_negative_short_series_is_removed():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_rt_iup

    rows = []
    for unsat, intercept in [(0, 1.2), (4, -0.65)]:
        for x in [16, 18, 20, 22, 24]:
            rows.append(
                {
                    "细类": "SM(d18:0/x:y)",
                    "曲线不饱和度": float(unsat),
                    "x碳数": float(x),
                    "归一化保留时间": 0.4 * x + intercept,
                }
            )
    for x, y in [(20, 8.45), (22, 8.30)]:
        rows.append(
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": float(x),
                "归一化保留时间": y,
            }
        )

    result = fit_rt_iup(pd.DataFrame(rows), RTIUPConfig(r2_threshold=0.98))

    assert result.rows[result.rows["曲线不饱和度"].eq(2.0)].empty
    assert result.plot_rows[result.plot_rows["曲线不饱和度"].eq(2.0)].empty
    assert result.stats["点数不足过滤行数"] == 2


def test_fitted_curve_requires_two_iup_supported_points():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, apply_fitted_point_iup_guard

    lines = pd.DataFrame(
        [
            {
                "细类": "SM(d18:0/x:y)",
                "不饱和度": 1.0,
                "拟合成功": True,
                "是否有效曲线": True,
                "拟合类型": "Linear",
                "参数": [0.4, 0.8],
                "x最小": 18.0,
                "x最大": 22.0,
                "去除原因": "",
            },
            {
                "细类": "SM(d18:0/x:y)",
                "不饱和度": 2.0,
                "拟合成功": True,
                "是否有效曲线": True,
                "拟合类型": "Linear",
                "参数": [0.4, 0.3],
                "x最小": 18.0,
                "x最大": 22.0,
                "去除原因": "",
            },
            {
                "细类": "SM(d18:0/x:y)",
                "不饱和度": 3.0,
                "拟合成功": True,
                "是否有效曲线": True,
                "拟合类型": "Linear",
                "参数": [0.4, -0.2],
                "x最小": 18.0,
                "x最大": 22.0,
                "去除原因": "",
            },
        ]
    )
    rows = pd.DataFrame(
        [
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": 18.0,
                "归一化保留时间": 7.5,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
                "去除原因": "",
            },
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": 20.0,
                "归一化保留时间": 9.2,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
                "去除原因": "",
            },
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": 22.0,
                "归一化保留时间": 8.1,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
                "去除原因": "",
            },
        ]
    )

    guarded_rows, guarded_lines = apply_fitted_point_iup_guard(rows, lines, RTIUPConfig())

    removed = guarded_lines[guarded_lines["不饱和度"].eq(2.0)].iloc[0]
    assert not bool(removed["是否有效曲线"])
    assert "拟合曲线IUP支持点不足" in removed["去除原因"]
    assert guarded_rows["是否作图"].sum() == 0
    assert guarded_rows["拟合点IUP支持数"].iloc[0] == 1
    assert guarded_rows["拟合点IUP判断数"].iloc[0] == 3
