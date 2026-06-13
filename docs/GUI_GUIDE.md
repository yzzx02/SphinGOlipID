# GUI Guide

The GUI is planned for v0.4. The current v0.3 backend exposes the API that the GUI should call.

## Backend Entry Point

GUI code should construct `SphinGOlipIDConfig` from user controls and call:

```python
from sphingolipid_toolkit import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import run_batch

results = run_batch(config)
```

The GUI should not call functions in `ms2_legacy_core.py` directly.

## Planned Streamlit Pages

1. Project setup: paths, project name, ion mode.
2. MS2 matching parameters: ppm, intensity cutoff, min fragments, score cutoff, top-N.
3. Run monitor: progress, logs, warnings, result paths.
4. Results review: table filters and downloads.
5. Visualization: class distribution, score distributions, matched-fragment plot.
6. RT validation: planned after backend integration.
7. Result cleaning: planned after backend integration.

## Logging

Use `logging_utils.attach_memory_handler()` to collect backend log messages for display in the GUI.

## Result Display

Prefer the standard alias columns added in v0.3:

```text
lipid_name
observed_mz
observed_rt
matched_fragment_count
ms2_match_score
rank
```

