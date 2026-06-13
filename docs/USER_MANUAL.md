# User Manual

## Install

From the project root:

```bash
python -m pip install -e .
```

For local development without installation, the project test config already adds `src` to the Python path for `pytest`.

## Run MS2 matching from CLI

```bash
sphingolipid-ms2 \
  --input-dir ./data/raw_ms2 \
  --precursor-dir ./data/raw_ms2 \
  --output-dir ./outputs/ms2_results \
  --ms1-db ./data/library/MS1_DB_new_3.0.xlsx \
  --files 1-6
```

Useful parameters:

```text
--fragment-ppm 20
--min-intensity 20
--min-fragments 2
--min-score 0.35
--top-n 3
--save-intermediate
--log-level INFO
```

## Run MS2 matching from Python

```python
from sphingolipid_toolkit import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import run_batch

config = SphinGOlipIDConfig(
    raw_ms2_dir="./data/raw_ms2",
    precursor_dir="./data/raw_ms2",
    ms1_library_path="./data/library/MS1_DB_new_3.0.xlsx",
    output_dir="./outputs/ms2_results",
    file_indices=(1, 2, 3, 4, 5, 6),
)

result_files = run_batch(config)
```

## Validate the Project

```bash
python -m pytest -q
python -m compileall src tests
```

## Run the Included Mock Dataset

```bash
sphingolipid-ms2 \
  --input-dir tests/data/raw_ms2 \
  --precursor-dir tests/data/precursors \
  --ms1-db tests/data/library/mock_ms1_library.xlsx \
  --output-dir tests/output \
  --fragment-ppm 20 \
  --min-intensity 20 \
  --min-fragments 2 \
  --min-score 0.35 \
  --top-n 3
```

Expected output files:

```text
tests/output/ms2_annotation_results.xlsx
tests/output/run_log.txt
tests/output/intermediate/
```

## Current Scope

v0.3 stabilizes the backend API. Streamlit GUI, formula library construction, RT validation integration, and final desktop packaging are planned follow-up stages.
