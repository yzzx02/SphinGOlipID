import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from sphingolipid_toolkit.scoring import (
    _sequence_length, add_standard_score_columns, parse_fragment_mz_values,
    match_fragments, score_fragment_matches,
)


@pytest.mark.parametrize("value,expected", [
    ([184.0733,264.2686,282.2792],3), ((100.,200.),2), (np.array([100.,200.]),2),
    (np.array(100.),1), (100.,1), (np.float64(100),1), (np.nan,0), (None,0),
    (pd.NA,0), ("",0), ("   ",0), ([],0), (np.array([]),0),
    ("[184.0733, 264.2686, 282.2792]",3), ("184.0733,264.2686,282.2792",3),
    ("(184.0733, 264.2686)",2), (" [ 100 ; 200 ; 300 ] ",3),
    ("[100. 200. 300.]",3), ("array([100., 200.], dtype=float64)",2),
    ("np.array([1e2, 2e2])",2), ("100; 200",2), ("100",1),
    ("[]",0), ("None",0), ("NaN",0), ([100, None, np.nan],1),
    ([100.,100.],2),  # The same observed m/z can be recorded twice by legacy.
])
def test_fragment_cell_parser_counts_real_entries(value, expected):
    assert len(parse_fragment_mz_values(value)) == expected
    assert _sequence_length(value) == expected
    out = add_standard_score_columns(pd.DataFrame({"实际mz":[value]}))
    assert out.matched_fragment_count.iloc[0] == expected


@pytest.mark.parametrize("value", ["100,garbage", "[100, broken]", "__import__('os').system('bad')", True, [[100,200]]])
def test_malformed_fragment_cell_is_not_silently_counted(value):
    with pytest.raises(ValueError):
        parse_fragment_mz_values(value)


def test_normalized_count_does_not_change_legacy_scores_or_rank_order(tmp_path):
    original = pd.DataFrame({
        "注释":["A","B"], "实际mz":["[100, 200, 300]","100; 200"],
        "target":[500.,500.], "匹配度分数":[.75,.5], "总分数":[75.,45.],
        "强度总和":[300.,200.],
    })
    # Exercise Excel string round-trip using synthetic cells only.
    path = tmp_path / "synthetic.xlsx"
    original.to_excel(path,index=False)
    loaded = pd.read_excel(path)
    result = add_standard_score_columns(loaded)
    assert_frame_equal(result[list(loaded)],loaded)
    assert result.matched_fragment_count.tolist() == [3,2]
    assert result["rank"].tolist() == [1,2]


def test_api_assigns_each_observed_peak_once():
    observed = pd.DataFrame({"fragment_mz":[100.,100.00005],"fragment_intensity":[50.,80.]})
    theoretical = pd.DataFrame({"theoretical_mz":[100.0001,100.0002,100.0001]})
    matches = match_fragments(observed,theoretical,20)
    summary = score_fragment_matches(matches,total_fragment_intensity=130)
    assert summary["matched_fragment_count"] == 2
    assert summary["matched_intensity_sum"] == 130.
    assert sorted(matches.observed_mz.tolist()) == [100.,100.00005]


def test_api_keeps_candidate_identity_when_masses_are_shared():
    observed = pd.DataFrame({"fragment_mz":[100.],"fragment_intensity":[50.]})
    theory = pd.DataFrame({"theoretical_mz":[100.,100.],"anno":["A","B"]})
    matches = match_fragments(observed,theory,20)
    assert set(matches.anno) == {"A","B"}
    assert len(matches) == 2
