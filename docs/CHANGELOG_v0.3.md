# Changelog v0.3

## Summary

v0.3 stabilizes the backend API before GUI work. The focus is configuration, input validation, logging, output layout, tests, and documentation around the converted MS2 workflow.

## Added Files

- `src/sphingolipid_toolkit/config.py`
- `src/sphingolipid_toolkit/io_utils.py`
- `src/sphingolipid_toolkit/logging_utils.py`
- `src/sphingolipid_toolkit/scoring.py`
- `docs/USER_MANUAL.md`
- `docs/INPUT_FORMAT.md`
- `docs/OUTPUT_FORMAT.md`
- `docs/GUI_GUIDE.md`
- `docs/DEVELOPER_NOTES.md`
- `docs/CHANGELOG_v0.3.md`
- `tests/test_config_io_scoring.py`
- `tests/test_ms2_pipeline.py`
- `tests/data/...` mock CLI input files

## Modified Files

- `src/sphingolipid_toolkit/__init__.py`
- `src/sphingolipid_toolkit/cli.py`
- `src/sphingolipid_toolkit/ms2_pipeline.py`
- `src/sphingolipid_toolkit/ms2_legacy_core.py`
- `configs/ms2_match.example.yaml`
- `docs/ARTICLE_ALIGNMENT.md`
- `docs/GUI_PLAN.md`
- `examples/README.md`
- `README.md`
- `pyproject.toml`
- `tests/test_import.py`

## MS2 Legacy Logic Not Changed

The following algorithmic pieces remain owned by `ms2_legacy_core.py` and were not rewritten:

- class-specific theoretical MS/MS fragment generation;
- LCB/headgroup/glycan neutral-loss rule tables;
- measured-theoretical fragment matching using ppm windows;
- original matched-fragment coverage score;
- original candidate ranking by total matched intensity and match score.

## Fixes and Stabilization

- Added `SphinGOlipIDConfig` so CLI, future GUI, and tests share one config shape.
- Added explicit validation for missing files, missing required columns, and non-numeric input values.
- Added standardized MS2 txt parsing output columns for downstream GUI use.
- Added CLI `--save-intermediate` as an alias for `--keep-intermediate`.
- Added `run_log.txt` output and console/file logging support.
- Added final aggregate result file `ms2_annotation_results.xlsx`.
- Separated intermediate files under `intermediate/`.
- Switched one intermediate Excel writer call from `xlsxwriter` to `openpyxl` so the minimal backend run works with the existing required Excel dependency.
- Added pytest coverage for config, IO validation, standardized MS2 parsing, scoring helpers, and the v0.3 output layout.

## Current Limitations

- The included mock dataset exercises the no-candidate path; it validates the workflow and output layout, but not full chemical annotation quality.
- Some legacy functions still print status messages directly; these should be routed through logging in a later cleanup.
- `run_batch()` still returns per-file result paths for backward compatibility, while also writing the aggregate `ms2_annotation_results.xlsx`.
- RT validation, result cleaning, and formula library generation are not yet wired into the GUI.
- No Streamlit GUI is included in v0.3.

## v0.4 Plan

v0.4 will add a Streamlit prototype in `gui/streamlit_app.py` with only:

- Project setup;
- MS2 parameters;
- Run MS2 matching;
- Results review;
- Basic visualization.

The GUI must call `SphinGOlipIDConfig` and `run_batch(config)` only. It must not implement MS2 txt parsing, fragment matching, or scoring logic directly.

