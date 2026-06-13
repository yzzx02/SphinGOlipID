import pytest
import pandas as pd

from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig, run_batch, run_one_file


def test_run_one_file_reports_missing_ms2_text(tmp_path):
    config = SphinGOlipIDConfig.from_single_input_dir(
        input_dir=tmp_path / "missing",
        output_dir=tmp_path / "out",
        ms1_library_path=tmp_path / "library.xlsx",
        file_indices=(1,),
    )

    with pytest.raises(FileNotFoundError, match="Missing required MS/MS txt file"):
        run_one_file(config, 1)


def test_legacy_ms2_pipeline_config_exposes_unified_names(tmp_path):
    config = MS2PipelineConfig(
        input_dir=tmp_path / "raw",
        precursor_dir=tmp_path / "targets",
        output_dir=tmp_path / "out",
        ms1_db_path=tmp_path / "library.xlsx",
        file_indices=[1],
    )

    assert config.raw_ms2_dir == tmp_path / "raw"
    assert config.precursor_dir == tmp_path / "targets"
    assert config.ms1_library_path == tmp_path / "library.xlsx"
    assert config.file_indices == (1,)


def test_run_batch_writes_v03_acceptance_outputs(tmp_path):
    raw_dir = tmp_path / "raw"
    precursor_dir = tmp_path / "precursors"
    output_dir = tmp_path / "out"
    raw_dir.mkdir()
    precursor_dir.mkdir()
    _write_mock_ms2_text(raw_dir / "hilic-msms-1.txt")
    _write_mock_precursor_table(precursor_dir / "hilic-msms-1.xlsx")
    library = tmp_path / "library.xlsx"
    _write_mock_ms1_library(library)

    config = SphinGOlipIDConfig(
        raw_ms2_dir=raw_dir,
        precursor_dir=precursor_dir,
        ms1_library_path=library,
        output_dir=output_dir,
        file_indices=(1,),
    )
    results = run_batch(config)

    assert results == [output_dir / "result_msms-1.xlsx"]
    assert (output_dir / "ms2_annotation_results.xlsx").is_file()
    assert (output_dir / "run_log.txt").is_file()
    assert (output_dir / "intermediate").is_dir()


def _write_mock_ms2_text(path):
    path.write_text(
        "\n".join(
            [
                "spectrum:",
                "index:1",
                "scan start time,1.25,min",
                "target m/z,760.5, m/z",
                "binaryDataArray:",
                "binary: 0 100.0 150.0 ",
                "binaryDataArray:",
                "binary: 0 25.0 30.0 ",
            ]
        ),
        encoding="utf-8",
    )


def _write_mock_precursor_table(path):
    pd.DataFrame(
        {
            "target": [760.5],
            "下限": [760.4],
            "上限": [760.6],
            "RT": [1.25],
            "Abund": [1000],
        }
    ).to_excel(path, index=False)


def _write_mock_ms1_library(path):
    pd.DataFrame(
        {
            "理论值": [999.0],
            "classy": ["Cer"],
            "name": ["Cer(d18:1/24:0)"],
            "structure": [""],
        }
    ).to_excel(path, index=False)
