# Conversion notes

This package is Python-only. The original notebooks were not preserved.

Converted sources:

| Original source | Python replacement |
|---|---|
| `MS2 Match (1).ipynb` | `ms2_legacy_core.py`, `ms2_pipeline.py`, `cli.py` |
| `version3.2plus.ipynb` | `legacy_converted/ms2_version3_2plus.py` using the configurable MS2 pipeline |
| `rentetion time.ipynb` | `rt_validation.py` |
| `RT_liner.ipynb` | `rt_validation.py` |
| `去除重复值.ipynb` | `result_cleaning.py` |

Notebook-specific problems fixed during conversion:

1. Removed `%matplotlib inline`, which is invalid in normal Python files.
2. Removed top-level hard-coded file paths such as `D:/...`, `C:/...`, and `E:/...`.
3. Replaced deprecated `DataFrame.append` with `pd.concat`.
4. Added missing imports required by the duplicate-removal notebook.
5. Removed or replaced an incomplete duplicate-removal cell that had a missing closing parenthesis.
6. Replaced `ExcelWriter.save()` with context-manager based writing.
7. Added numeric conversion before m/z and intensity comparisons.
8. Added empty-result handling for MS2 batches with no candidate or no matched fragment.
