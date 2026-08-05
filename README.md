# SphinGOlipID Python-only Project

This is the Python-only version of the SphinGOlipID sphingolipid annotation toolkit. The uploaded notebooks have been converted into `.py` modules; no original `.ipynb` files are kept in this package.

## Current implemented functions

### 0. v0.3 backend framework

The current backend framework adds the project-level API requested in `plan.txt`:

- unified `SphinGOlipIDConfig` in `src/sphingolipid_toolkit/config.py`;
- input validation helpers in `io_utils.py` for missing files, missing Excel columns, and non-numeric m/z/RT/intensity values;
- package logging helpers in `logging_utils.py`, including a memory handler for future GUI display;
- normalized scoring/result helpers in `scoring.py`;
- `run_batch(config)` accepts the unified config while preserving the older `MS2PipelineConfig` compatibility layer;
- `pytest` and `python -m compileall src tests` are expected to pass from the project root.

### 1. MS2 matching workflow

Implemented in:

- `src/sphingolipid_toolkit/config.py`
- `src/sphingolipid_toolkit/io_utils.py`
- `src/sphingolipid_toolkit/logging_utils.py`
- `src/sphingolipid_toolkit/scoring.py`
- `src/sphingolipid_toolkit/ms2_legacy_core.py`
- `src/sphingolipid_toolkit/ms2_pipeline.py`
- `src/sphingolipid_toolkit/cli.py`

Supported workflow:

1. Read Agilent-style MS/MS txt exports.
2. Read precursor target Excel files.
3. Extract target precursor spectra.
4. Search MS1 theoretical database for candidate isomers.
5. Generate rule-based theoretical MS/MS fragments.
6. Match measured fragments to theoretical fragments by ppm tolerance.
7. Score candidates by matched-fragment coverage and fragment intensity.
8. Export top-N candidates to Excel.

Default parameters match the original notebook:

- fragment tolerance: `20 ppm`
- minimum measured fragment intensity: `20`
- minimum matched fragments: `2`
- minimum MS2 match score: `0.35`
- top candidates per precursor: `3`

Non-algorithmic fixes already applied:

- removed hard-coded Windows paths;
- replaced top-level notebook execution with reusable functions;
- fixed `out_DB.csv` and `out-query.csv` readback so the first generated row is not lost as a header;
- added numeric conversion for measured m/z and intensity columns before filtering;
- added graceful empty-result output when no MS1 candidates, theoretical fragments, or fragment matches are found;
- added pre-run validation for required target/MS1 columns;
- added stable English alias columns such as `lipid_name`, `observed_mz`, `matched_fragment_count`, `ms2_match_score`, and `rank`;
- added final aggregate output `ms2_annotation_results.xlsx`, run log `run_log.txt`, and separated `intermediate/` directory;
- closed Excel handles after matching.

### 2. RT validation workflow

Implemented in:

- `src/sphingolipid_toolkit/rt_validation.py`

Functions include:

- merge multiple Excel result tables;
- parse lipid annotations into homologous `series` and variable carbon number `x`;
- fit RANSAC RT models for homologous series;
- calculate predicted RT and RT deviation;
- find missing homologous-series members;
- match predicted RT/mass candidates back to target tables.

This module replaces the original `rentetion time.ipynb` and `RT_liner.ipynb` notebooks.

### 3. Duplicate removal / result cleaning

Implemented in:

- `src/sphingolipid_toolkit/result_cleaning.py`

Functions include:

- remove close RT duplicates within the same `series`/`x` group;
- keep higher MS2 match score, then higher total fragment intensity;
- reset homologous-series indices;
- filter target/abundance duplicate candidates;
- clean an Excel file directly.

This module replaces the original `去除重复值.ipynb` notebook.

### 4. Notebook-name-compatible Python modules

For traceability, the old notebook names have Python replacements under:

```text
src/sphingolipid_toolkit/legacy_converted/
```

These are `.py` modules only. They re-export the cleaned production functions and avoid invalid notebook syntax such as `%matplotlib inline`, hard-coded file paths, and incomplete cells.

## Folder layout

```text
SphinGOlipID_project_pyonly_v0.2/
├── src/
│   └── sphingolipid_toolkit/
│       ├── __init__.py
│       ├── config.py
│       ├── io_utils.py
│       ├── logging_utils.py
│       ├── scoring.py
│       ├── ms2_legacy_core.py
│       ├── ms2_pipeline.py
│       ├── cli.py
│       ├── formula_library.py
│       ├── rt_validation.py
│       ├── result_cleaning.py
│       └── legacy_converted/
├── configs/
│   └── ms2_match.example.yaml
├── docs/
│   ├── GUI_PLAN.md
│   ├── INPUT_FORMAT.md
│   ├── OUTPUT_FORMAT.md
│   ├── USER_MANUAL.md
│   ├── GUI_GUIDE.md
│   ├── DEVELOPER_NOTES.md
│   └── ARTICLE_ALIGNMENT.md
├── examples/
├── tests/
├── requirements.txt
├── pyproject.toml
└── README.md
```

## Installation

```bash
cd SphinGOlipID_project_pyonly_v0.2
python -m pip install -e .
```

## Run MS2 matching

The default file naming follows the original MS2 notebook:

- `hilic-msms-1.txt` ... `hilic-msms-6.txt`
- `hilic-msms-1.xlsx` ... `hilic-msms-6.xlsx`

```bash
sphingolipid-ms2 \
  --input-dir ./data/raw_ms2 \
  --precursor-dir ./data/raw_ms2 \
  --output-dir ./outputs/ms2_results \
  --ms1-db ./data/library/MS1_DB_new_3.0.xlsx \
  --files 1-6 \
  --fragment-ppm 20 \
  --min-intensity 20 \
  --min-fragments 2 \
  --min-score 0.35 \
  --top-n 3 \
  --save-intermediate
```

Default final outputs in `--output-dir`:

```text
ms2_annotation_results.xlsx
run_log.txt
intermediate/
result_msms-{i}.xlsx
```

## Run Desktop GUI

The Tkinter desktop GUI keeps the earlier lab software layout while calling the
current `run_batch(config)` backend:

```bash
python -m sphingolipid_toolkit.gui
```

After editable installation, the console script is also available:

```bash
sphingolipid-gui
```

The GUI accepts one raw TXT file, one precursor XLSX file, one MS1 database XLSX
file, and an output folder. The selected single TXT/XLSX pair is mapped to the
backend's configurable file-name patterns and processed as file index `1`.

## Retention Time IUP Validation

The cleaned RT/IUP workflow is available from `sphingolipid_toolkit.rt_iup`.
It fits linear or quadratic RANSAC models for each lipid subclass and
unsaturation, removes obvious IUP ordering violations, rescues strict candidates
that sit between adjacent valid IUP curves, and draws final RT plots with solid
fits and capped 95% confidence bands.

```python
import pandas as pd
from sphingolipid_toolkit.rt_iup import fit_rt_iup, prepare_rank_table, write_rt_iup_plots

rank_table = pd.read_excel("SphinGOlipID_rank5_normalized_raw_2d_rt.xlsx")
prepared = prepare_rank_table(rank_table, dataset_name="top5")
result = fit_rt_iup(prepared)
write_rt_iup_plots(result.plot_rows, result.lines, "rt_plots")
```

The current six-fraction workflow treats ECN and IUP as complementary RT
validation modes:

- ECN fits each unsaturation series independently after collapsing repeated
  RT values at the same carbon number to their median. Outliers are removed
  iteratively and the final model is refit on the stable inlier set.
- IUP searches the original candidates independently within the locked
  `0.20 min` RT window and `0.05 min` IUP allowance. Parallelism and ordering
  guide candidate selection but do not silently relax either threshold.
- Linear models are preferred unless a quadratic model has at least four
  distinct carbon numbers and improves R² by at least `0.002`.
- Plot legends report only the unsaturation and fitted R².

## Six-fraction supplementary information

The compact SI release is under
`supplementary_data/SI_6fractions_20260805/`. It contains the six-fraction
initial TOP3 identifications, the three ECN/IUP classification summaries,
the retained RT results with explicit high/low score tiers, and the fitted-line table.
The wide internal QA tables and raw Agilent `.d`/converted `.mzML` files are
not included.

The flat SI tables can be regenerated with:

```bash
python scripts/export_si_core_tables.py \
  --input path/to/TOP3_RT_input_6_fractions_corrected.csv \
  --results-root path/to/final_RT_results \
  --output-dir supplementary_data/SI_6fractions_20260805 \
  --include-low-score
```

The six-fraction reconstruction and figure workflow is implemented in:

- `scripts/prepare_missing_fraction_workbooks.py`
- `scripts/rebuild_six_fraction_initial_tables.py`
- `scripts/build_top3_rt_input_from_fixed_csv.py`
- `scripts/regenerate_top3_rt_ecn_iup_preview.py`
- `scripts/build_final_filter_summary.py`
- `scripts/build_ms1_rt_rescue.py`
- `scripts/export_si_core_tables.py`

## Use from Python

```python
from sphingolipid_toolkit import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import run_batch

config = SphinGOlipIDConfig(
    raw_ms2_dir="./data/raw_ms2",
    precursor_dir="./data/raw_ms2",
    output_dir="./outputs/ms2_results",
    ms1_library_path="./data/library/MS1_DB_new_3.0.xlsx",
    file_indices=(1, 2, 3, 4, 5, 6),
)

result_files = run_batch(config)
```

The older `MS2PipelineConfig(input_dir=..., ms1_db_path=...)` form is still supported for compatibility.

## Verify the backend

```bash
python -m pytest -q
python -m compileall src tests
```

## RT validation example

```python
import pandas as pd
from sphingolipid_toolkit.rt_validation import add_series_columns, fit_ransac_models

df = pd.read_excel("annotation_results.xlsx")
df = add_series_columns(df, annotation_col="注释", drop_unparsed=True)
# Make sure the observed RT column is named y, or pass y_col="RT".
inliers, models = fit_ransac_models(df, x_col="x", y_col="RT")
```

## Result cleaning example

```python
import pandas as pd
from sphingolipid_toolkit.result_cleaning import remove_duplicate_annotations

df = pd.read_excel("annotation_results.xlsx")
cleaned = remove_duplicate_annotations(df, annotation_col="注释", rt_col="RT")
cleaned.to_excel("annotation_results_cleaned.xlsx", index=False)
```

## What is still reserved for later

- Full theoretical formula-library construction from class templates.
- GUI implementation.
- More end-to-end tests once real example input data are provided.

## Targeted mzML workflow

For Agilent targeted MS2 mzML data with multiple targetlist CSV files:

```bash
sphingolipid-targeted-mzml \
  --data-dir "E:/李奕晓数据备份/20240422-Agilent/2D-MSMS" \
  --ms1-db "path/to/MS1_DB_new_3.0.xlsx" \
  --output-dir ./outputs/targeted_mzml \
  --top-n 3
```

See `docs/TARGETED_MZML_WORKFLOW.md`.
