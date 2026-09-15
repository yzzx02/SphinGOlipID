from itertools import permutations

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from sphingolipid_toolkit.result_cleaning import deduplicate_by_rt_score_intensity


def frame(rts, scores, intensities=None):
    return pd.DataFrame({
        "series": ["Cer"] * len(rts), "x": [34] * len(rts), "y": rts,
        "匹配度分数": scores, "强度总和": intensities or [100] * len(rts),
        "label": [chr(65+i) for i in range(len(rts))],
    })


@pytest.mark.parametrize("scores,expected", [([.9, .7], "A"), ([.7, .9], "B")])
def test_best_evidence_wins_in_either_rt_direction(scores, expected):
    data = frame([10., 10.05], scores)
    original = data.copy(deep=True)
    for order in permutations(range(len(data))):
        result = deduplicate_by_rt_score_intensity(data.iloc[list(order)])
        assert result.label.tolist() == [expected]
    assert_frame_equal(data, original)


def test_separate_chromatographic_features_are_retained():
    result = deduplicate_by_rt_score_intensity(frame([10., 10.2], [.9, .7]))
    assert set(result.label) == {"A", "B"}


@pytest.mark.parametrize("rts,scores,expected", [
    ([10., 10.04, 10.08], [.7, .9, .8], {"B"}),
    # Proximity is not transitive: two endpoints may be distinct features.
    ([10., 10.08, 10.16], [.9, .8, .7], {"A", "C"}),
    ([10., 10.08, 10.16], [.8, .9, .7], {"B"}),
])
def test_cluster_selection_checks_every_retained_rt(rts, scores, expected):
    data = frame(rts, scores)
    for order in permutations(range(len(data))):
        result = deduplicate_by_rt_score_intensity(data.iloc[list(order)])
        assert set(result.label) == expected
        assert all(abs(a-b) > .1 for a,b in permutations(result.y, 2))


def test_intensity_then_rt_then_content_break_ties_independent_of_row_order():
    data = frame([10.04, 10., 10.], [.9]*3, [90, 100, 100])
    for order in permutations(range(3)):
        result = deduplicate_by_rt_score_intensity(data.iloc[list(order)].reset_index(drop=True))
        assert result.label.tolist() == ["B"]
    result = deduplicate_by_rt_score_intensity(frame([10.05, 10.], [.9, .9]))
    assert result.label.tolist() == ["B"]


def test_groups_and_inclusive_tolerance_boundary():
    data = frame([10., 10.125, 10., 10.], [.9, .8, .7, .6])
    data.loc[2,"series"] = "SM"
    data.loc[3,"x"] = 36
    reference = deduplicate_by_rt_score_intensity(data, rt_tolerance=.125)
    assert set(reference.label) == {"A", "C", "D"}
    assert_frame_equal(reference, deduplicate_by_rt_score_intensity(data.iloc[::-1], rt_tolerance=.125))


def test_content_tiebreak_does_not_round_away_distinct_metadata():
    data = frame([10.,10.],[.9,.9])
    data["label"] = "same"
    data["extra"] = [1.,1.0000000000000002]
    first = deduplicate_by_rt_score_intensity(data)
    reversed_rows = deduplicate_by_rt_score_intensity(data.iloc[::-1].reset_index(drop=True))
    assert_frame_equal(first,reversed_rows)


def test_empty_missing_rt_and_invalid_tolerance():
    data = frame([None, None], [.9,.8])
    assert len(deduplicate_by_rt_score_intensity(data)) == 2
    assert list(deduplicate_by_rt_score_intensity(data.iloc[:0])) == list(data)
    with pytest.raises(ValueError, match="non-negative"):
        deduplicate_by_rt_score_intensity(data, rt_tolerance=-.1)
