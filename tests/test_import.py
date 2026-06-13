def test_import_package():
    import sphingolipid_toolkit
    from sphingolipid_toolkit import SphinGOlipIDConfig

    assert sphingolipid_toolkit.__version__
    assert SphinGOlipIDConfig


def test_parse_file_indices():
    from sphingolipid_toolkit.ms2_pipeline import parse_file_indices
    assert parse_file_indices("1-3") == (1, 2, 3)
    assert parse_file_indices("1,3,5") == (1, 3, 5)


def test_rt_series_parser():
    from sphingolipid_toolkit.rt_validation import parse_lipid_series

    series, x = parse_lipid_series("Cer(d18:1/24:0)")
    assert series == "Cer(d18:1/x:0)"
    assert x == 24


def test_result_cleaning_deduplicate():
    import pandas as pd
    from sphingolipid_toolkit.result_cleaning import deduplicate_by_rt_score_intensity

    df = pd.DataFrame(
        {
            "series": ["A", "A", "A"],
            "x": [20, 20, 21],
            "y": [1.00, 1.05, 2.00],
            "匹配度分数": [0.5, 0.6, 0.7],
            "强度总和": [10, 8, 5],
        }
    )
    out = deduplicate_by_rt_score_intensity(df)
    assert len(out) == 2
    assert out.loc[out["x"] == 20, "匹配度分数"].iloc[0] == 0.6


def test_no_notebooks_in_project_tree():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    assert not list(root.rglob("*.ipynb"))
