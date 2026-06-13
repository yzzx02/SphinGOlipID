from pathlib import Path

import pandas as pd
import pytest

from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.io_utils import (
    InputValidationError,
    parse_ms2_text_to_dataframe,
    read_ms2_txt_standardized,
    read_excel_table,
    validate_ms2_file_inputs,
    validate_required_columns,
)
from sphingolipid_toolkit.scoring import match_fragments, score_fragment_matches


def test_unified_config_normalizes_paths_and_indices(tmp_path):
    config = SphinGOlipIDConfig.from_single_input_dir(
        input_dir=tmp_path / "input",
        output_dir=tmp_path / "output",
        ms1_library_path=tmp_path / "library.xlsx",
        file_indices=[1, "2"],
    )

    assert config.raw_ms2_dir == tmp_path / "input"
    assert config.precursor_dir == tmp_path / "input"
    assert config.ms1_library_path == tmp_path / "library.xlsx"
    assert config.file_indices == (1, 2)
    assert config.fragment_tolerance_fraction == 20 / 1_000_000


def test_read_excel_table_reports_missing_column(tmp_path):
    path = tmp_path / "targets.xlsx"
    pd.DataFrame({"target": [760.5]}).to_excel(path, index=False)

    with pytest.raises(InputValidationError, match="Missing required column: RT"):
        read_excel_table(path, required_columns=("target", "RT"), context="target table")


def test_read_excel_table_reports_non_numeric_value(tmp_path):
    path = tmp_path / "targets.xlsx"
    pd.DataFrame({"target": ["not-a-number"], "RT": [1.2]}).to_excel(path, index=False)

    with pytest.raises(InputValidationError, match="Column target contains non-numeric value"):
        read_excel_table(path, required_columns=("target", "RT"), numeric_columns=("target",), context="target table")


def test_parse_ms2_text_to_dataframe(tmp_path):
    text_file = tmp_path / "hilic-msms-1.txt"
    text_file.write_text(
        "\n".join(
            [
                "spectrum:",
                "index:42",
                "scan start time,1.25,min",
                "target m/z,760.5, m/z",
                "binaryDataArray:",
                "binary: 0 100.0 150.0",
                "binaryDataArray:",
                "binary: 0 20.0 30.0",
            ]
        ),
        encoding="utf-8",
    )

    parsed = parse_ms2_text_to_dataframe(text_file, encoding="utf-8")
    parsed_from_alias = read_ms2_txt_standardized(text_file, encoding="utf-8")

    assert list(parsed.columns) == [
        "file_name",
        "scan_id",
        "precursor_mz",
        "rt",
        "fragment_mz",
        "fragment_intensity",
    ]
    assert len(parsed) == 2
    assert parsed_from_alias.equals(parsed)
    assert parsed["precursor_mz"].iloc[0] == 760.5
    assert parsed["fragment_intensity"].tolist() == [20.0, 30.0]


def test_validate_required_columns_alias_reports_missing_column():
    with pytest.raises(InputValidationError, match="Missing required column: RT"):
        validate_required_columns(pd.DataFrame({"target": [760.5]}), ("target", "RT"), "target table")


def test_validate_ms2_file_inputs_checks_all_required_tables(tmp_path):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "hilic-msms-1.txt").write_text("spectrum:\n", encoding="utf-8")
    pd.DataFrame(
        {
            "target": [760.5],
            "下限": [760.4],
            "上限": [760.6],
            "RT": [1.2],
            "Abund": [1000],
        }
    ).to_excel(raw_dir / "hilic-msms-1.xlsx", index=False)
    library = tmp_path / "library.xlsx"
    pd.DataFrame(
        {
            "理论值": [760.5],
            "classy": ["Cer"],
            "name": ["Cer(d18:1/24:0)"],
            "structure": [""],
        }
    ).to_excel(library, index=False)

    config = SphinGOlipIDConfig.from_single_input_dir(
        input_dir=raw_dir,
        output_dir=tmp_path / "out",
        ms1_library_path=library,
        file_indices=(1,),
    )
    paths, target_df = validate_ms2_file_inputs(config, 1)

    assert paths.text_file == raw_dir / "hilic-msms-1.txt"
    assert len(target_df) == 1
    assert (tmp_path / "out").is_dir()


def test_match_fragments_and_score_summary():
    observed = pd.DataFrame(
        {
            "fragment_mz": [100.001, 150.0, 200.0],
            "fragment_intensity": [50.0, 5.0, 25.0],
        }
    )
    theoretical = pd.DataFrame(
        {
            "theoretical_mz": [100.0, 200.002],
            "fragment_name": ["diagnostic", "lcb"],
            "fragment_type": ["diagnostic", "lcb"],
        }
    )

    matches = match_fragments(observed, theoretical, ppm_tolerance=20, min_intensity=20)
    summary = score_fragment_matches(matches, total_fragment_intensity=100)

    assert len(matches) == 2
    assert summary["matched_fragment_count"] == 2
    assert summary["matched_intensity_sum"] == 75.0
    assert summary["matched_intensity_ratio"] == 0.75
