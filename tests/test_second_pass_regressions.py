"""Additional edge cases on the merged GUI/algorithm version."""

import base64
import logging
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.io_utils import InputValidationError, write_excel_atomic
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig, run_batch
from sphingolipid_toolkit.rt_validation import match_predicted_to_targets
from sphingolipid_toolkit.result_cleaning import filter_target_abundance_duplicates
from sphingolipid_toolkit.scoring import score_fragment_matches
from sphingolipid_toolkit.targeted_mzml_pipeline import (
    MS2Spectrum, TargetedMzMLConfig, _read_legacy_or_normalized_ms1_library,
    read_targetlist_csv, run_targeted_mzml_batch,
)
from test_bug_regressions import _write_inputs


@pytest.mark.parametrize("indices", [[1.9], [True], [float("inf")]])
@pytest.mark.parametrize("config_class", [SphinGOlipIDConfig, MS2PipelineConfig])
def test_fractional_or_boolean_file_indices_are_rejected(tmp_path, indices, config_class):
    with pytest.raises(ValueError, match="integer"):
        if config_class is SphinGOlipIDConfig:
            config_class.from_single_input_dir(tmp_path, tmp_path / "out", tmp_path / "library.xlsx", file_indices=indices)
        else:
            config_class(tmp_path, tmp_path / "out", tmp_path / "library.xlsx", file_indices=indices)


def test_result_and_log_collision_is_rejected_before_modifying_files(tmp_path):
    config = _write_inputs(tmp_path)
    conflicting = SphinGOlipIDConfig.from_single_input_dir(tmp_path, config.output_dir, config.ms1_library_path,
                                                         file_indices=(1,), log_file_name="ms2_annotation_results.xlsx")
    config.output_dir.mkdir()
    original = config.output_dir / conflicting.log_file_name
    original.write_bytes(b"previous result")
    with pytest.raises(InputValidationError, match="distinct"):
        run_batch(conflicting)
    assert original.read_bytes() == b"previous result"


def test_failed_analysis_keeps_previous_result_workbook(tmp_path):
    config = _write_inputs(tmp_path)
    result = run_batch(config)[0]
    original = result.read_bytes()
    (tmp_path / "hilic-msms-1.txt").write_text("damaged export", encoding="utf-8")
    with pytest.raises(InputValidationError, match="no spectrum"):
        run_batch(config)
    assert result.read_bytes() == original


def test_missing_later_input_is_detected_before_first_result_is_overwritten(tmp_path):
    config = _write_inputs(tmp_path)
    result = run_batch(config)[0]
    original = result.read_bytes()
    batch = SphinGOlipIDConfig.from_single_input_dir(tmp_path, config.output_dir, config.ms1_library_path, file_indices=(1, 2))
    with pytest.raises(FileNotFoundError, match="index 2"):
        run_batch(batch)
    assert result.read_bytes() == original


def test_atomic_excel_write_preserves_previous_file_on_writer_failure(tmp_path, monkeypatch):
    path = tmp_path / "results.xlsx"
    pd.DataFrame({"value": [1]}).to_excel(path, index=False)
    original = path.read_bytes()
    def fail_during_write(self, temporary, **kwargs):
        Path(temporary).write_bytes(b"partial workbook")
        raise OSError("write interrupted")
    monkeypatch.setattr(pd.DataFrame, "to_excel", fail_during_write)
    with pytest.raises(OSError, match="interrupted"):
        write_excel_atomic(pd.DataFrame({"value": [2]}), path, index=False)
    assert path.read_bytes() == original
    assert not list(tmp_path.glob(".sphingolipid_result_*"))


def test_rt_matching_preserves_both_tables_when_columns_overlap():
    predicted = pd.DataFrame({"pred_y": [1., np.nan], "MH_T": [500., 500.], "annotation": ["candidate", "missing RT"], "annotation_target": ["existing", "existing"]})
    targets = pd.DataFrame({"RT": [1., np.nan], "Precur": [500., 500.], "annotation": ["feature", "missing RT"]})
    result = match_predicted_to_targets(predicted, targets)
    assert len(result) == 1
    assert result.columns.is_unique
    assert result.annotation.tolist() == ["candidate"]
    assert result.annotation_target.tolist() == ["existing"]
    assert result.annotation_target_target.tolist() == ["feature"]


def test_score_summary_does_not_sum_unmatched_intensity():
    matches = pd.DataFrame({"matched": [True, False, pd.NA], "observed_intensity": [10., 100., 1000.]})
    result = score_fragment_matches(matches, total_fragment_intensity=1110.)
    assert result["matched_fragment_count"] == 1
    assert result["matched_intensity_sum"] == 10.
    assert result["matched_intensity_ratio"] == pytest.approx(10 / 1110)


@pytest.mark.parametrize("scores,expected", [(["9", "10"], "best"), ([None, None], "first")])
def test_duplicate_filter_handles_repeated_indices_and_missing_scores(scores, expected):
    data = pd.DataFrame({"target": [500., 500.], "Abund": [100., 100.], "index": [1, 1],
                         "总分数": scores, "name": ["first", "best"]}, index=[5, 5])
    result = filter_target_abundance_duplicates(data)
    assert result.name.tolist() == [expected]


def test_bom_targetlist_and_numeric_enabled_values(tmp_path):
    path = tmp_path / "targets.csv"
    path.write_text("\nTargetedMSMSTable,,,,\nOn,Prec. m/z,Ret. Time (min)\n1,500,1.25\n,600,2\n", encoding="utf-8-sig")
    result = read_targetlist_csv(path)
    assert result.target_mz.tolist() == [500.]


def test_empty_targetlist_reports_a_clear_input_error(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("\n\n", encoding="utf-8-sig")
    with pytest.raises(InputValidationError, match="Empty targetlist"):
        read_targetlist_csv(path)


@pytest.mark.parametrize("masses,intensities", [([100, 200], [50]), ([[100, 200]], [[50, 60]])])
def test_invalid_mzml_arrays_are_rejected(masses, intensities):
    with pytest.raises(InputValidationError, match="equal lengths"):
        MS2Spectrum(Path("sample.mzML"), "scan=1", 500., 1., None, np.array(masses), np.array(intensities))


def test_invalid_native_ms1_library_is_not_silently_reinterpreted(tmp_path):
    path = tmp_path / "library.xlsx"
    pd.DataFrame({"理论值": ["invalid"], "classy": ["Cer"], "name": ["Cer(d34:1)"], "structure": [""]}).to_excel(path, index=False)
    with pytest.raises(InputValidationError, match="non-numeric"):
        _read_legacy_or_normalized_ms1_library(path)


def _write_mzml(path):
    """One standards-shaped MS2 spectrum with actual base64 binary arrays."""
    namespace = "http://psi.hupo.org/ms/mzml"
    ET.register_namespace("", namespace)
    def node(parent, tag, **attributes):
        return ET.SubElement(parent, f"{{{namespace}}}{tag}", attributes)
    def cv(parent, accession, name, value="", **attributes):
        return node(parent, "cvParam", cvRef="MS", accession=accession, name=name, value=value, **attributes)
    root = ET.Element(f"{{{namespace}}}mzML", version="1.1.0", id="synthetic")
    cvs = node(root, "cvList", count="1")
    node(cvs, "cv", id="MS", fullName="PSI-MS", version="4.1.0", URI="https://purl.obolibrary.org/obo/ms.obo")
    run = node(root, "run", id="run1")
    spectra = node(run, "spectrumList", count="1", defaultDataProcessingRef="DP1")
    spectrum = node(spectra, "spectrum", id="scan=1", index="0", defaultArrayLength="3")
    cv(spectrum, "MS:1000511", "ms level", "2")
    scans = node(spectrum, "scanList", count="1")
    scan = node(scans, "scan")
    cv(scan, "MS:1000016", "scan start time", "75", unitCvRef="UO", unitAccession="UO:0000010", unitName="second")
    precursors = node(spectrum, "precursorList", count="1")
    precursor = node(precursors, "precursor")
    ions = node(precursor, "selectedIonList", count="1")
    ion = node(ions, "selectedIon")
    cv(ion, "MS:1000744", "selected ion m/z", "500")
    activation = node(precursor, "activation")
    cv(activation, "MS:1000045", "collision energy", "25")
    arrays = node(spectrum, "binaryDataArrayList", count="2")
    for values, accession, name in (([264.2686, 282.2792, 60.044], "MS:1000514", "m/z array"), ([25., 30., 40.], "MS:1000515", "intensity array")):
        encoded = base64.b64encode(struct.pack("<3d", *values)).decode("ascii")
        array = node(arrays, "binaryDataArray", encodedLength=str(len(encoded)))
        cv(array, "MS:1000523", "64-bit float")
        cv(array, "MS:1000576", "no compression")
        cv(array, accession, name)
        node(array, "binary").text = encoded
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def _targeted_inputs(tmp_path):
    targets = tmp_path / "targetlist"
    targets.mkdir()
    (targets / "targets.CSV").write_text("TargetedMSMSTable,,,\nOn,Prec. m/z,Ret. Time (min)\nTRUE,500,1.25\n", encoding="utf-8-sig")
    _write_mzml(tmp_path / "sample.MZML")
    library = tmp_path / "library.xlsx"
    pd.DataFrame({"理论值": [500.], "classy": ["Cer"], "name": ["Cer(d34:1)"], "structure": [""]}).to_excel(library, index=False)
    return TargetedMzMLConfig(tmp_path, library, tmp_path / "out")


def test_targeted_batch_reads_real_mzml_binary_data_and_bom_targetlist(tmp_path):
    config = _targeted_inputs(tmp_path)
    logger = logging.getLogger(f"test.real_mzml.{tmp_path.name}")
    outputs = run_targeted_mzml_batch(config, logger=logger)
    assert all(path.is_file() for path in outputs.values())
    raw = pd.read_excel(outputs["raw_identifications"])
    assert not raw.empty
    assert raw.ms2_rt.tolist() == [1.25] * len(raw)
    assert raw.matched_fragment_count.max() >= 2
    assert raw.ms1_library_mz_error_ppm.notna().all()
    assert raw.ms1_library_mz_error_ppm.max() <= config.ms1_candidate_ppm
    assert not pd.read_excel(outputs["rt_fitted_identifications"]).empty
    assert not any(isinstance(handler, logging.FileHandler) for handler in logger.handlers)


def test_targeted_output_cannot_overwrite_library(tmp_path):
    config = _targeted_inputs(tmp_path)
    original = config.ms1_library_path.read_bytes()
    config = TargetedMzMLConfig(tmp_path, config.ms1_library_path, tmp_path, feature_table_name="library.xlsx")
    with pytest.raises(InputValidationError, match="overwrite"):
        run_targeted_mzml_batch(config)
    assert config.ms1_library_path.read_bytes() == original
    assert not (tmp_path / config.log_file_name).exists()


def test_targeted_missing_mzml_is_reported_before_export(tmp_path):
    config = _targeted_inputs(tmp_path)
    (tmp_path / "sample.MZML").unlink()
    with pytest.raises(InputValidationError, match="requires"):
        run_targeted_mzml_batch(config)
    assert not config.output_dir.exists()


@pytest.mark.parametrize("changes", [{"feature_mz_ppm": -1}, {"feature_rt_tolerance_min": float("nan")}, {"rt_fit_residual_threshold": 0}])
def test_targeted_matching_settings_are_validated(tmp_path, changes):
    with pytest.raises(ValueError):
        TargetedMzMLConfig(tmp_path, tmp_path / "library.xlsx", tmp_path / "out", **changes)
