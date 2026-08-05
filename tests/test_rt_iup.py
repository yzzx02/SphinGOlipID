from pathlib import Path

import pandas as pd


def _sm_d18_1_regression_fixture() -> pd.DataFrame:
    """Representative SM(d18:1) points retained from the real analysis shape."""

    values = {
        0.0: {
            14: [6.2795], 16: [7.044], 17: [7.362], 18: [7.791],
            20: [8.621], 22: [9.293], 23: [9.895], 24: [10.0625], 25: [10.729],
        },
        # Fourteen raw records at six distinct carbon numbers.  The ECN fit
        # must use the six per-carbon medians, not treat replicates as six
        # additional independent fit points.
        1.0: {
            16: [6.445], 19: [7.613], 22: [8.764, 8.879, 8.897],
            23: [9.156, 9.180, 9.295], 24: [9.395, 9.413, 9.580],
            25: [9.841, 9.962, 9.978],
        },
        # Three distinct points: linear is valid and quadratic is forbidden.
        2.0: {23: [8.579], 24: [8.9905], 25: [9.372]},
        3.0: {
            17: [6.6325], 18: [7.0345], 20: [7.736], 22: [8.612],
            23: [9.089], 24: [9.439], 25: [9.897], 26: [10.2735],
        },
    }
    records = []
    row_id = 1
    for unsaturation, by_carbon in values.items():
        for carbon, retention_times in by_carbon.items():
            for retention_time in retention_times:
                records.append(
                    {
                        "合并后行ID": f"SM_FIX_{row_id:03d}",
                        "细类": "SM(d18:1/x:y)",
                        "曲线不饱和度": unsaturation,
                        "x碳数": float(carbon),
                        "归一化保留时间": float(retention_time),
                    }
                )
                row_id += 1
    return pd.DataFrame(records)


def test_sm_d18_1_expected_ecn_and_iup_lines():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_ecn, fit_rt_iup

    fixture = _sm_d18_1_regression_fixture()
    config = RTIUPConfig()
    ecn = fit_ecn(fixture, config)
    iup = fit_rt_iup(fixture, config)

    assert set(ecn.lines["不饱和度"]) == {0.0, 1.0, 2.0, 3.0}
    assert set(iup.lines["不饱和度"]) == {0.0, 1.0, 3.0}

    unsat_1 = ecn.lines[ecn.lines["不饱和度"].eq(1.0)].iloc[0]
    assert unsat_1["点数"] == 14
    assert unsat_1["不同碳数"] == 6
    assert unsat_1["R²"] >= 0.99

    unsat_2 = ecn.lines[ecn.lines["不饱和度"].eq(2.0)].iloc[0]
    assert unsat_2["拟合类型"] == "Linear"
    assert unsat_2["不同碳数"] == 3

    parallel = ecn.lines[ecn.lines["不饱和度"].isin([1.0, 2.0, 3.0])]
    assert parallel["平行性是否通过"].eq(True).all()


def test_ecn_keeps_individually_valid_nonparallel_curves():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_ecn

    rows = []
    for unsat, slope, intercept in [(0.0, 0.40, 1.0), (1.0, 0.70, -5.5)]:
        for carbon in [18, 20, 22, 24, 26]:
            rows.append(
                {
                    "细类": "Cer",
                    "曲线不饱和度": unsat,
                    "x碳数": float(carbon),
                    "归一化保留时间": slope * carbon + intercept,
                }
            )

    result = fit_ecn(pd.DataFrame(rows), RTIUPConfig())

    assert set(result.lines["不饱和度"]) == {0.0, 1.0}
    assert result.lines["平行性是否通过"].eq(False).sum() == 1
    assert result.lines["是否有效曲线"].eq(True).all()


def test_iup_can_rescue_curve_that_ecn_could_not_fit():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_ecn, fit_rt_iup

    rows = []
    for unsat, intercept in [(0.0, 1.0), (2.0, 0.0)]:
        for carbon in [18, 20, 22]:
            rows.append(
                {
                    "合并后行ID": f"boundary-{unsat}-{carbon}",
                    "细类": "Cer",
                    "曲线不饱和度": unsat,
                    "x碳数": float(carbon),
                    "归一化保留时间": 0.4 * carbon + intercept,
                }
            )
    good = {18: 7.7, 20: 8.5, 22: 9.3}
    bad = {18: 4.0, 20: 13.0, 22: 5.0}
    for carbon in [18, 20, 22]:
        for label, rt in [("good", good[carbon]), ("bad", bad[carbon])]:
            rows.append(
                {
                    "合并后行ID": f"middle-{label}-{carbon}",
                    "细类": "Cer",
                    "曲线不饱和度": 1.0,
                    "x碳数": float(carbon),
                    "归一化保留时间": rt,
                }
            )

    frame = pd.DataFrame(rows)
    config = RTIUPConfig()
    ecn = fit_ecn(frame, config)
    iup = fit_rt_iup(frame, config)

    assert 1.0 not in set(ecn.lines["不饱和度"])
    rescued = iup.lines[iup.lines["不饱和度"].eq(1.0)].iloc[0]
    assert bool(rescued["IUP重新寻优"])
    assert rescued["救回不同碳数"] == 3


def test_iterative_refit_recomputes_final_model_after_outlier_removal():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_line

    points = pd.DataFrame(
        {
            "细类": ["SM(d18:1/x:y)"] * 6,
            "曲线不饱和度": [1.0] * 6,
            "x碳数": [16.0, 18.0, 20.0, 22.0, 24.0, 26.0],
            "归一化保留时间": [6.4, 7.2, 10.5, 8.8, 9.6, 10.4],
        }
    )
    line = fit_line(points, RTIUPConfig(r2_threshold=0.99))

    assert line["拟合成功"]
    assert line["拟合类型"] == "Linear"
    assert line["内点数"] == 5
    assert abs(line["参数"][0] - 0.4) < 1e-10
    assert abs(line["参数"][1]) < 1e-10
    assert line["R²"] == 1.0


def test_legend_label_reports_only_fit_r_squared():
    from sphingolipid_toolkit.rt_iup import legend_label

    group = pd.DataFrame(
        {
            "x碳数": [32.0, 34.0, 36.0],
            "归一化保留时间": [7.0, 7.9, 8.8],
        }
    )
    line = pd.Series({"拟合类型": "Linear", "R²": 0.996714, "内点数": 3})

    assert legend_label(1, group, line) == "1  R²=0.9967"
    assert legend_label(1, group.iloc[:2]) == "1"


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
                "归一化保留时间": 9.31,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
                "去除原因": "",
            },
            {
                "细类": "SM(d18:0/x:y)",
                "曲线不饱和度": 2.0,
                "x碳数": 22.0,
                "归一化保留时间": 8.0,
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


def _short_series_fixture(higher_unsat: float) -> pd.DataFrame:
    rows = []
    for unsat, intercept in [(0, 2.0), (higher_unsat, 0.5)]:
        for x in [18, 20, 22, 24, 26]:
            rows.append(
                {
                    "细类": "SM d",
                    "曲线不饱和度": float(unsat),
                    "x碳数": float(x),
                    "归一化保留时间": 0.4 * x + intercept,
                }
            )
    target_unsat = 2.0
    fraction = target_unsat / higher_unsat
    target_intercept = 2.0 - fraction * (2.0 - 0.5)
    for x in [20, 22]:
        rows.append(
            {
                "细类": "SM d",
                "曲线不饱和度": target_unsat,
                "x碳数": float(x),
                "归一化保留时间": 0.4 * x + target_intercept,
            }
        )
    return pd.DataFrame(rows)


def test_two_point_short_series_may_cross_one_unsaturation_only():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, fit_rt_iup

    allowed = fit_rt_iup(_short_series_fixture(3.0), RTIUPConfig(r2_threshold=0.98))
    blocked = fit_rt_iup(_short_series_fixture(4.0), RTIUPConfig(r2_threshold=0.98))

    assert allowed.plot_rows[allowed.plot_rows["曲线不饱和度"].eq(2.0)].shape[0] == 2
    assert blocked.rows[blocked.rows["曲线不饱和度"].eq(2.0)].empty
    assert blocked.stats["点数不足过滤行数"] == 2


def test_unsaturation_7_color_differs_from_1():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, color_for_unsaturation

    config = RTIUPConfig()
    assert color_for_unsaturation(1, config) != color_for_unsaturation(7, config)


def test_relative_rt_tolerance_keeps_locked_thresholds_and_allows_opposite_drift():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, relative_rt_tolerance_min

    config = RTIUPConfig()
    assert config.rt_window_min == 0.2
    assert config.iup_tolerance_min == 0.05
    assert relative_rt_tolerance_min(config) == 0.45

    assert 8.486 - 8.047 <= relative_rt_tolerance_min(config)
    assert 10.900 - 9.984 > relative_rt_tolerance_min(config)


def test_point_level_iup_guard_keeps_good_boundary_points():
    from sphingolipid_toolkit.rt_iup import RTIUPConfig, apply_point_level_iup_guard

    lines = pd.DataFrame(
        [
            {
                "细类": "SM d",
                "不饱和度": 6.0,
                "是否有效曲线": True,
                "拟合类型": "Linear",
                "参数": [0.28, -3.32],
                "x最小": 37.0,
                "x最大": 44.0,
            },
            {
                "细类": "SM d",
                "不饱和度": 7.0,
                "是否有效曲线": True,
                "拟合类型": "Linear",
                "参数": [0.09, 4.4],
                "x最小": 37.0,
                "x最大": 44.0,
            },
        ]
    )
    rows = pd.DataFrame(
        [
            {
                "细类": "SM d",
                "曲线不饱和度": 7.0,
                "x碳数": 37.0,
                "归一化保留时间": 7.733,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
            },
            {
                "细类": "SM d",
                "曲线不饱和度": 7.0,
                "x碳数": 44.0,
                "归一化保留时间": 8.361,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
            },
            {
                "细类": "SM d",
                "曲线不饱和度": 7.0,
                "x碳数": 44.0,
                "归一化保留时间": 8.380,
                "同曲线0.2min内点": True,
                "是否作图": True,
                "IUP说明": "有效拟合曲线，0.2min内点",
            },
        ]
    )

    guarded = apply_point_level_iup_guard(rows, lines, RTIUPConfig())

    assert not bool(guarded.loc[0, "同曲线0.2min内点"])
    assert not bool(guarded.loc[0, "是否作图"])
    assert guarded.loc[1:, "同曲线0.2min内点"].sum() == 2
    assert guarded.loc[1:, "x碳数"].nunique() == 1
