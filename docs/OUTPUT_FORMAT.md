# Output Format

The TXT command writes one result workbook per processed index and one
aggregate workbook. Direct mzML writes one identification workbook. Targeted
mzML writes raw and RT-filtered identifications, the targetlist, and the MS1
feature table. The desktop and non-targeted CLI write a six-sheet workbook
documented in the [desktop guide](GUI_GUIDE.md).

Result workbooks are written to temporary files and replace the previous
workbook after a successful write. A failed rerun preserves the previous result;
use the run status/log to distinguish an old workbook from a new completed run.

Default TXT output file names:

```text
ms2_annotation_results.xlsx
run_log.txt
result_msms-1.xlsx
result_msms-2.xlsx
...
result_msms-6.xlsx
```

## Legacy result columns

Compatibility columns for existing review workflows include:

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
```

## Standard alias columns

The current engine also exports stable English columns:

```text
lipid_name
observed_mz
observed_rt
signal_intensity
matched_fragment_count
matched_intensity_sum
total_score
ms2_match_score
rank
```

`total_score` is the current 0–100 score. `ms2_match_score` stores it on a 0–1
scale for compatibility. `result_rank` and `resolution_level` record the
engine's ranking and reporting scope. `ms2_evidence_json` retains fragment
evidence. When it exceeds the Excel cell limit, the row records its digest
and a reference to an adjacent `ms2-evidence-*.json.gz` file.

## Intermediate files

For TXT, `--keep-intermediate` retains `intermediate/parsed_ms2-{i}.csv` with
the parsed scan IDs, precursor m/z, RT and measured fragments. For non-targeted
mzML, `--save-intermediate` or the desktop **保留 XCMS 峰表** option retains
per-file centWave outputs beside the workbook.
