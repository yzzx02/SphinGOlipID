import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LinearRegression

from sphingolipid_toolkit import rt_iup as rt


def points(x, y):
    return pd.DataFrame({"细类":"synthetic Cer", "曲线不饱和度":1.,
                         "x碳数":x, "归一化保留时间":y})


def test_ransac_excludes_outlier_and_final_model_is_new_ols(monkeypatch):
    real_ransac = rt.RANSACRegressor
    seeds = []

    class ObservedSeed(real_ransac):
        def fit(self, x, y, **kwargs):
            super().fit(x,y,**kwargs)
            seeds.append(self)
            # Canary: only the real RANSAC inlier mask should be consumed.
            self.estimator_.intercept_ = 1000.
            return self

    monkeypatch.setattr(rt,"RANSACRegressor",ObservedSeed)
    data = points([16.,18.,20.,22.,24.,26.], [6.4,7.2,10.5,8.8,9.6,10.4])
    model, stable = rt.iterative_refit(data,"Linear")
    assert not seeds[0].inlier_mask_[2]
    assert stable["x碳数"].tolist() == [16.,18.,22.,24.,26.]
    expected = LinearRegression().fit(stable[["x碳数"]].to_numpy(),stable["归一化保留时间"].to_numpy())
    assert model is not seeds[0].estimator_
    assert model.coef_ == pytest.approx(expected.coef_)
    assert model.intercept_ == pytest.approx(expected.intercept_)
    assert model.coef_[0] == pytest.approx(.4)


def test_replicates_are_medians_not_independent_fit_points():
    data = points([16,16,16,18,20,22], [6.0,6.4,9.0,7.2,8.,8.8])
    model, stable = rt.iterative_refit(data,"Linear")
    assert len(stable) == 4
    assert stable.loc[stable["x碳数"].eq(16),"归一化保留时间"].iloc[0] == 6.4
    assert model.coef_[0] == pytest.approx(.4)


@pytest.mark.parametrize("curvature,expected", [(0.,"Linear"),(.004,"Linear"),(.01,"Quadratic")])
def test_linear_quadratic_selection_with_real_synthetic_fits(curvature, expected):
    x = np.array([16.,18.,20.,22.,24.])
    data = points(x,.4*x + curvature*(x-20)**2)
    config = rt.RTIUPConfig()
    linear = rt._fit_candidate_model(data,"Linear",config)
    quadratic = rt._fit_candidate_model(data,"Quadratic",config)
    assert linear is not None and quadratic is not None
    gain = quadratic["R²"]-linear["R²"]
    assert (gain >= .002) == (expected == "Quadratic")
    line = rt.fit_line(data)
    assert line["拟合类型"] == expected
    assert line["拟合成功"]


def test_quadratic_requires_four_distinct_carbons_not_four_rows():
    data = points([16,16,18,20], [6.4,6.4,7.2,8.])
    assert rt._fit_candidate_model(data,"Quadratic",rt.RTIUPConfig()) is None
    assert rt.fit_line(data)["拟合类型"] == "Linear"


@pytest.mark.parametrize("linear_r2,expected", [(.9980001,"Linear"),(.998,"Quadratic")])
def test_model_selection_gain_boundary(monkeypatch,linear_r2,expected):
    # Isolate the decision boundary from numerical fitting noise. Real model
    # construction and validation are exercised in the tests above.
    def candidate(data,fit_type,config):
        return {"拟合类型":fit_type,"R²":linear_r2 if fit_type == "Linear" else 1.}
    monkeypatch.setattr(rt,"_fit_candidate_model",candidate)
    assert rt.fit_line(points([16,18,20,22],[6.4,7.2,8.,8.8]))["拟合类型"] == expected


@pytest.mark.parametrize("fit_type", ["Linear","Quadratic"])
def test_r2_below_point_99_rejected(fit_type):
    x = np.array([16.,18.,20.,22.,24.])
    data = points(x,3+.04*(x-16)+np.array([0,.05,-.05,.05,0]))
    assert rt.iterative_refit(data,fit_type) is not None
    assert rt._fit_candidate_model(data,fit_type,rt.RTIUPConfig()) is None


def test_decreasing_rt_is_rejected_despite_perfect_r2():
    x = np.array([16.,18.,20.,22.,24.])
    line = rt.fit_line(points(x,12-.4*x))
    assert not line["拟合成功"]


def test_defaults_and_summary_aliases_preserve_chinese_fields():
    config = rt.RTIUPConfig()
    assert (config.r2_threshold,config.rt_window_min,config.quadratic_min_distinct_x,config.quadratic_min_r2_gain) == (.99,.2,4,.002)
    data = points([16.,18.,20.,22.],[6.4,7.2,8.,8.8])
    summaries = [rt.fit_line(data),rt.fit_line(data.iloc[:2])]
    summaries.extend(rt.fit_ecn(data).lines.to_dict("records"))
    summaries.extend(rt.fit_rt_iup(data).lines.to_dict("records"))
    for line in summaries:
        assert line["fit_strategy"] == "RANSAC-seeded OLS refit"
        for alias, original in {
            "fit_type":"拟合类型","r2":"R²", "representative_point_count":"代表点数",
            "inlier_count":"内点数","x_min":"x最小","x_max":"x最大",
        }.items():
            assert line[alias] == line[original]


def test_existing_iup_rescue_ranking_is_not_the_fit_line_gain_rule(monkeypatch):
    """Known scope gap: rescue compares parallelism, then R2, not a 0.002 gate.

    Do not silently change this existing RT selection path while adding aliases.
    """
    data = points([16,18,20,22],[6.4,7.2,8.,8.8])
    monkeypatch.setattr(rt,"pick_rescue_points",lambda *args:data)
    monkeypatch.setattr(rt,"is_near_horizontal",lambda *args:False)
    monkeypatch.setattr(rt,"line_between_bracket",lambda *args:True)
    monkeypatch.setattr(rt,"pair_parallel_status",lambda *args:(True,{"斜率相对差":0.}))
    def candidate(points,fit_type,config):
        return {"拟合类型":fit_type,"R²":.999 if fit_type == "Linear" else .9995,
                "参数":[.4,0.] if fit_type == "Linear" else [0.,.4,0.],
                "内点数":4,"代表点数":4,"x最小":16.,"x最大":22.}
    monkeypatch.setattr(rt,"_fit_rescue_candidate",candidate)
    result = rt.build_rescue_line("synthetic",1.,data,pd.Series({"不饱和度":0.}),pd.Series({"不饱和度":2.}))
    assert result["拟合类型"] == "Quadratic"
