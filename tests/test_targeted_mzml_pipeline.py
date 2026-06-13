from pathlib import Path

import numpy as np
import pandas as pd

from sphingolipid_toolkit.targeted_mzml_pipeline import (
    MS2Spectrum,
    TargetedMzMLConfig,
    apply_rt_fit_filter,
    build_ms1_feature_table,
    finalize_spectrum_matches,
    find_matching_feature,
    match_observed_to_theoretical,
    read_targetlist_csv,
    search_ms1_candidates,
)


def test_read_targetlist_csv_and_build_feature_table(tmp_path):
    csv_path = tmp_path / "0-5min-1.csv"
    csv_path.write_text(
        "\n".join(
            [
                "TargetedMSMSTable,,,,,,,",
                "On,Prec. m/z,Z,Ret. Time (min),Delta Ret. Time (min),Iso. Width,Collision Energy,Acquisition Time (ms/spec)",
                "TRUE,329.3278,1,5.039,0.4,Medium (~4 m/z),25,",
                "TRUE,329.3280,1,5.041,0.4,Medium (~4 m/z),40,",
                "TRUE,500.0000,1,6.000,0.4,Medium (~4 m/z),25,",
            ]
        ),
        encoding="utf-8",
    )

    raw = read_targetlist_csv(csv_path)
    features = build_ms1_feature_table(raw, mz_tolerance_da=0.005, rt_tolerance_min=0.05)

    assert len(raw) == 3
    assert len(features) == 2
    first = features.sort_values("feature_mz").iloc[0]
    assert first["targetlist_row_count"] == 2
    assert abs(first["feature_mz"] - 329.3279) < 1e-5


def test_match_observed_to_theoretical_and_finalize(tmp_path):
    spectrum = MS2Spectrum(
        mzml_file=Path("sample.mzML"),
        scan_id="scanId=1",
        precursor_mz=329.3278,
        ms2_rt=5.04,
        collision_energy=25.0,
        fragment_mz=np.array([100.0, 150.0, 200.0]),
        fragment_intensity=np.array([30.0, 5.0, 40.0]),
    )
    theoretical = pd.DataFrame(
        {
            "idf": ["frag_a", "frag_b"],
            "mz": [100.0, 200.0],
            "anno": ["Cer(d18:1/10:0)", "Cer(d18:1/10:0)"],
            "total": ["Cer(d28:1)", "Cer(d28:1)"],
            "target": [329.3278, 329.3278],
            "lenDB": [2, 2],
            "RT": [5.039, 5.039],
            "Abund": [1.0, 1.0],
        }
    )

    matched = match_observed_to_theoretical(
        spectrum,
        theoretical,
        fragment_tolerance_fraction=20 / 1_000_000,
        min_fragment_intensity=20,
    )
    result = finalize_spectrum_matches(matched, min_matched_fragments=2, min_match_score=0.35, top_n=3)

    assert len(matched) == 2
    assert len(result) == 1
    assert result["注释"].iloc[0] == "Cer(d18:1/10:0)"
    assert result["匹配度分数"].iloc[0] == 1.0
    assert result["fragment_error_ppm_max"].iloc[0] == 0.0


def test_ppm_filters_for_feature_and_ms1_candidate():
    features = pd.DataFrame(
        {
            "feature_id": ["F1", "F2"],
            "feature_mz": [500.004, 500.006],
            "ms1_rt": [5.0, 5.0],
        }
    )
    matched = find_matching_feature(features, precursor_mz=500.0, ms2_rt=5.0, mz_tolerance_ppm=10, rt_tolerance_min=0.1)
    assert matched is not None
    assert matched["feature_id"] == "F1"
    assert matched["feature_mz_error_ppm"] <= 10

    library = pd.DataFrame(
        {
            "理论值": [500.004, 500.006],
            "name": ["in", "out"],
            "classy": ["Cer", "Cer"],
            "structure": ["", ""],
        }
    )
    candidates = search_ms1_candidates(library, observed_mz=500.0, tolerance_ppm=10, ms1_rt=5.0)
    assert candidates["name"].tolist() == ["in"]
    assert candidates["ms1_library_mz_error_ppm"].iloc[0] <= 10


def test_rt_fit_filter_keeps_not_evaluated_rows():
    data = pd.DataFrame(
        {
            "identification_id": ["ID1", "ID2"],
            "lipid_name": ["Cer(d18:1/20:0)", "Cer(d18:1/21:0)"],
            "ms1_rt": [10.0, 10.3],
        }
    )
    config = TargetedMzMLConfig(
        data_dir=Path("."),
        ms1_library_path=Path("library.xlsx"),
        output_dir=Path("out"),
        rt_fit_min_points=3,
    )

    filtered = apply_rt_fit_filter(data, config)

    assert len(filtered) == 2
    assert set(filtered["rt_validation_status"]) == {"not_evaluated"}
