# Output Format

The MS2 backend writes one result Excel file per processed file index and one aggregate result file for the whole batch.

Default output file names:

```text
ms2_annotation_results.xlsx
run_log.txt
intermediate/
result_msms-1.xlsx
result_msms-2.xlsx
...
result_msms-6.xlsx
```

## Legacy result columns

The converted notebook columns are preserved:

```text
注释
target
匹配度分数
母离子
RT
Abund
实际mz
强度
碎片
强度总和
总分数
相对强度
```

## Standard alias columns

v0.3 adds stable English aliases for CLI, GUI, tests, and downstream modules:

```text
lipid_name
observed_mz
observed_rt
signal_intensity
matched_fragment_count
matched_intensity_sum
matched_intensity_ratio
ms2_match_score
rank
```

The aliases are derived from the legacy output and do not change the original scoring logic.

## Intermediate files

When `save_intermediate` / `--keep-intermediate` is enabled, the pipeline keeps files such as:

```text
actual_ms2_data-{i}.xlsx
DB_sheet-{i}.xlsx
Isomer.xlsx
out_DB.csv
out-query.csv
```

When disabled, temporary `Isomer.xlsx`, `out_DB.csv`, and `out-query.csv` are removed after each file.
