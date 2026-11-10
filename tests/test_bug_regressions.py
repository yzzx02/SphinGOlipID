import logging
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from pyteomics.auxiliary import unitfloat

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.io_utils import InputValidationError, ensure_output_dir, parse_ms2_spectrum_blocks
from sphingolipid_toolkit.ms2_pipeline import _finalize_results, _write_db_sheets, run_batch
from sphingolipid_toolkit.result_cleaning import deduplicate_by_rt_score_intensity, filter_target_abundance_duplicates
from sphingolipid_toolkit.rt_validation import fit_ransac_models, match_predicted_to_targets, predict_rt
from sphingolipid_toolkit.targeted_mzml_pipeline import TargetedMzMLConfig, _extract_scan_rt, apply_rt_fit_filter


def test_theoretical_fragments_and_scores_stay_with_each_precursor(tmp_path):
    rows = []
    for target, rt, abundance, intensity in ((500., 1., 100., 20.), (501., 1., 100., 30.), (500., 2., 200., 40.)):
        for mz in (100., 200.):
            rows.append(["fragment", mz, "Cer(d18:1/16:0)", "Cer(d34:1)", target, 2, rt, abundance, intensity])
    data = pd.DataFrame(rows)
    db_csv = tmp_path / "db.csv"
    data.iloc[:, :8].to_csv(db_csv, header=False, index=False)
    assert _write_db_sheets(db_csv, tmp_path / "db.xlsx", .00002)
    with pd.ExcelFile(tmp_path / "db.xlsx") as book:
        assert set(book.sheet_names) == {"500.0", "501.0"}
        assert len(book.parse("500.0")) == 4
    matches = [row[:8] + [row[1] * .99998, row[1] * 1.00002, 1, row[1], row[8]] for row in rows]
    query = tmp_path / "query.csv"
    pd.DataFrame(matches).to_csv(query, header=False, index=False)
    result = _finalize_results(query, tmp_path / "result.xlsx", 2, .35, 1)
    assert len(result) == 3
    assert result["matched_fragment_count"].tolist() == [2, 2, 2]
    assert result["rank"].tolist() == [1, 1, 1]
    assert set(zip(result.target, result.RT, result["强度总和"])) == {(500., 1., 40.), (501., 1., 60.), (500., 2., 80.)}


def test_assignment_preserves_same_annotation_at_different_precursors(tmp_path, monkeypatch):
    monkeypatch.setattr(core, "folder", str(tmp_path) + "/")
    observed = pd.DataFrame({"idx": [1], "A": [100.], "B": [50.]})
    theory = pd.DataFrame({"idf": ["first", "first"], "mz": [100., 100.], "anno": ["same", "same"],
                          "total": ["total", "total"], "target": [500., 501.], "lenDB": [1, 1],
                          "RT": [1., 1.], "Abund": [100., 100.], "下限": [99.998, 99.998], "上限": [100.002, 100.002]})
    core.process_sheet("synthetic", SimpleNamespace(parse=lambda _: observed.copy()), SimpleNamespace(parse=lambda _: theory.copy()))
    result = _finalize_results(tmp_path / "out-query.csv", tmp_path / "result.xlsx", 1, 0., 1)
    assert set(result.target) == {500., 501.}
    assert result.matched_fragment_count.tolist() == [1, 1]


def test_theoretical_label_merging_is_retained_from_latest_algorithm(tmp_path):
    rows = [[label, 100., "same", "total", target, 1, 1., 100.] for target in (500., 501.) for label in ("diagnostic", "support")]
    csv = tmp_path / "theory.csv"
    pd.DataFrame(rows).to_csv(csv, index=False, header=False)
    assert _write_db_sheets(csv, tmp_path / "theory.xlsx", .00002)
    with pd.ExcelFile(tmp_path / "theory.xlsx") as book:
        for sheet in book.sheet_names:
            frame = book.parse(sheet)
            assert frame.idf.tolist() == ["diagnostic | support"]
            assert len(frame) == 1


@pytest.mark.parametrize("scores,intensities,expected", [([.9, .8], [10, 20], 0), ([.8, .9], [20, 10], 1), ([.9, .9], [30, 20], 0), ([.9, .9], [20, 30], 1), ([.9, .9], [20, 20], 0)])
def test_dedup_keeps_best_row_regardless_of_rt_order(scores, intensities, expected):
    data = pd.DataFrame({"series": ["s", "s"], "x": [16, 16], "y": [1., 1.05], "匹配度分数": scores, "强度总和": intensities, "id": [0, 1]})
    result = deduplicate_by_rt_score_intensity(data)
    assert result.id.tolist() == [expected]


def test_series_frequency_tie_uses_total_score():
    data = pd.DataFrame({"target": [500, 500], "Abund": [100, 100], "index": [1, 2], "总分数": [10, 90]})
    assert filter_target_abundance_duplicates(data)["index"].tolist() == [2]


def test_text_parser_reports_misaligned_arrays():
    block = "index:1\ntarget m/z,500, m/z\nbinaryDataArray:\nbinary: 0 100 200\nbinaryDataArray:\nbinary: 0 50\n"
    with pytest.raises(InputValidationError, match="different lengths"):
        parse_ms2_spectrum_blocks([block])


@pytest.mark.parametrize("block", [
    "index:1\ntarget m/z,500, m/z\n",
    "index:1\ntarget m/z,nan, m/z\nbinaryDataArray:\nbinary: 0 100\nbinaryDataArray:\nbinary: 0 50\n",
    "index:1\ntarget m/z,500, m/z\nbinaryDataArray:\nbinary: 0 bad 200\nbinaryDataArray:\nbinary: 0 50 60\n",
])
def test_malformed_ms2_data_is_reported_instead_of_silent_success(block):
    with pytest.raises(InputValidationError):
        parse_ms2_spectrum_blocks([block])


def test_write_probe_preserves_existing_file(tmp_path):
    probe = tmp_path / ".sphingolipid_write_test"
    probe.write_text("user content")
    ensure_output_dir(tmp_path)
    assert probe.read_text() == "user content"


@pytest.mark.parametrize("settings", [{"top_n": 0}, {"fragment_ppm": float("nan")}, {"fragment_ppm": -1}, {"min_match_score": 1.1}, {"file_indices": ()}, {"file_indices": (1, 2), "result_pattern": "result.xlsx"}])
def test_invalid_settings_are_rejected(tmp_path, settings):
    config = dict(raw_ms2_dir=tmp_path, precursor_dir=tmp_path, ms1_library_path=tmp_path / "library.xlsx", output_dir=tmp_path / "out", file_indices=(1,))
    config.update(settings)
    with pytest.raises(ValueError):
        SphinGOlipIDConfig(**config)


def _write_inputs(tmp_path, target=500., library_mz=500.):
    txt = tmp_path / "hilic-msms-1.txt"
    txt.write_text("spectrum:\nindex:1\nscan start time,1.25,min\ntarget m/z,500, m/z\nbinaryDataArray:\nbinary: 0 264.2686 282.2792 60.044\nbinaryDataArray:\nbinary: 0 25 30 40\n", encoding="utf-8")
    pd.DataFrame({"target": [target], "下限": [target-.1], "上限": [target+.1], "RT": [1.25], "Abund": [1000]}).to_excel(tmp_path / "hilic-msms-1.xlsx", index=False)
    library = tmp_path / "library.xlsx"
    pd.DataFrame({"理论值": [library_mz], "classy": ["Cer"], "name": ["Cer(d34:1)"], "structure": [""]}).to_excel(library, index=False)
    return SphinGOlipIDConfig.from_single_input_dir(tmp_path, tmp_path / "out", library, file_indices=(1,))


def test_no_matching_precursor_produces_empty_result(tmp_path):
    config = _write_inputs(tmp_path, target=900., library_mz=900.)
    result = run_batch(config)
    assert pd.read_excel(result[0]).empty
    assert "lipid_name" in pd.read_excel(config.output_dir / config.final_result_name)


def test_numeric_ms1_strings_are_used_after_validation(tmp_path):
    config = _write_inputs(tmp_path, library_mz="500.0")
    result = run_batch(config)
    assert not pd.read_excel(result[0]).empty


def test_input_output_collision_does_not_modify_input(tmp_path):
    config = _write_inputs(tmp_path)
    original = config.ms1_library_path.read_bytes()
    config = SphinGOlipIDConfig.from_single_input_dir(tmp_path, tmp_path, config.ms1_library_path, file_indices=(1,), result_pattern="library.xlsx")
    with pytest.raises(InputValidationError, match="overwrite"):
        run_batch(config)
    assert config.ms1_library_path.read_bytes() == original


def test_run_log_handler_is_closed_on_success_and_error(tmp_path):
    config = _write_inputs(tmp_path)
    logger = logging.getLogger(f"test.run_handler.{tmp_path.name}")
    run_batch(config, logger=logger)
    assert not any(isinstance(handler, logging.FileHandler) for handler in logger.handlers)
    config.ms1_library_path.unlink()
    with pytest.raises(FileNotFoundError):
        run_batch(config, logger=logger)
    assert not any(isinstance(handler, logging.FileHandler) for handler in logger.handlers)
    (config.output_dir / config.log_file_name).unlink()


def test_rt_empty_and_degenerate_series_keep_schema():
    data = pd.DataFrame({"series": ["s", "s"], "x": [16, 16], "RT": [1., 3.]})
    inliers, models = fit_ransac_models(data, y_col="RT")
    assert inliers.empty and models.empty
    assert "series" in models and "slope" in models
    assert "RT" in inliers and "pred_y" in inliers
    assert predict_rt(data, models).pred_y.isna().all()


def test_mass_tolerance_does_not_include_extra_relative_tolerance():
    predicted = pd.DataFrame({"pred_y": [1.], "MH_T": [1000.]})
    targets = pd.DataFrame({"RT": [1.], "Precur": [1000.105]})
    assert match_predicted_to_targets(predicted, targets, mz_tolerance=.1).empty


def test_mzml_scan_time_seconds_are_normalized():
    spectrum = {"scanList": {"scan": [{"scan start time": unitfloat(120, "second")}]}}
    assert _extract_scan_rt(spectrum) == 2.


def test_failed_rt_series_are_not_removed_with_successful_series():
    data = pd.DataFrame({"identification_id": [f"ID{i}" for i in range(6)], "lipid_name": ["Cer(d18:1/20:0)", "Cer(d18:1/21:0)", "Cer(d18:1/22:0)"] + ["SM(d18:1/20:0)"] * 3, "ms1_rt": [10., 10.3, 10.6, 20., 21., 22.]})
    config = TargetedMzMLConfig(data_dir=Path("."), ms1_library_path=Path("library.xlsx"), output_dir=Path("out"))
    filtered = apply_rt_fit_filter(data, config)
    assert len(filtered) == 6
    assert set(filtered.loc[filtered.lipid_name.str.startswith("SM"), "rt_validation_status"]) == {"not_evaluated"}
