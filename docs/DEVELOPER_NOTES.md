# Developer Notes

## v0.3 Backend Boundary

The backend is organized around a stable public API:

- `SphinGOlipIDConfig` in `config.py`
- `run_batch(config)` and `run_one_file(config, file_index)` in `ms2_pipeline.py`
- validation helpers in `io_utils.py`
- log helpers in `logging_utils.py`
- standard result/score helpers in `scoring.py`

`ms2_legacy_core.py` still owns the converted rule-based MS/MS fragment generation and matching logic. Do not rewrite that chemistry logic unless a specific bug is identified and covered by a test.

## Input Validation

Run validation before entering the legacy core. Missing files, missing columns, and non-numeric target/MS1 values should raise clear exceptions from `io_utils.py`.

## Logging

New backend code should log through package loggers, usually:

```python
from sphingolipid_toolkit.logging_utils import get_logger

logger = get_logger("sphingolipid_toolkit.ms2")
```

Avoid adding new `print()` calls in backend code. Existing legacy prints can be removed or redirected later as part of cleanup.

## Tests

Required local checks:

```bash
python -m pytest -q
python -m compileall src tests
```

Tests should use small mock data and should not depend on large real instrument exports.

## Planned Next Steps

1. Add Streamlit prototype in `gui/streamlit_app.py`.
2. Connect RT validation outputs to stable result files.
3. Expand result-cleaning reports.
4. Implement formula-library generation after the MS2/GUI API is stable.

