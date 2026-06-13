# Targeted mzML Annotation Workflow

This workflow is for targeted Agilent mzML MS2 data where one MS1 feature may have multiple targeted MS2 spectra.

## Goal

1. Merge all `targetlist/*.csv` files into one MS1-like feature table.
2. Read every targeted mzML MS2 spectrum.
3. Match each MS2 spectrum back to the MS1 feature by precursor m/z and MS2 RT.
4. Run the legacy rule-based MS2 annotation and keep top 3 candidates per MS2 spectrum.
5. Annotate all top 1-3 results back to the MS1 feature row.
6. Preserve both RT values:
   - `ms1_rt`: unified feature RT from targetlist.
   - `ms2_rt`: actual MS2 spectrum RT from mzML.
7. Use all top 1-3 candidates in RT fitting, then remove only RANSAC-excluded outliers.
8. Keep rows from series with too few points for fitting as `not_evaluated`.

## Command

```bash
sphingolipid-targeted-mzml \
  --data-dir "E:/path/to/20240422-Agilent/2D-MSMS" \
  --ms1-db "E:/path/to/MS1_library.csv" \
  --output-dir "./outputs/targeted_mzml" \
  --ms1-ppm 10 \
  --feature-ppm 10 \
  --fragment-ppm 10 \
  --min-intensity 20 \
  --min-fragments 2 \
  --min-score 0.35 \
  --top-n 3
```

Optional:

```text
--targetlist-dir  Override data-dir/targetlist.
--mzml-dir        Override where mzML files are searched.
--ms1-ppm         MS1 library candidate tolerance in ppm.
--feature-ppm     mzML precursor to targetlist feature tolerance in ppm.
```

## Outputs

```text
combined_targetlist_raw.xlsx
combined_ms1_feature_table.xlsx
targeted_mzml_raw_identifications.xlsx
targeted_mzml_rt_fitted_identifications.xlsx
targeted_mzml_run_log.txt
```

## MS1 Library

The MS1 library may be either the normalized workbook used by the legacy pipeline or an Agilent MassHunter compound database CSV. For the MassHunter CSV format, neutral monoisotopic mass is computed from `Formula` when the `Mass` field is empty, then converted to positive-mode `[M+H]+` theoretical m/z for candidate lookup.

The MS1 library is required because the legacy MS2 matcher is a candidate-based rule system: precursor m/z is first used to find possible lipid structures, those structures generate theoretical fragments, and the observed MS2 spectrum is scored against those fragments. Without the library there are no candidate structures or theoretical fragments to score.

For strict 10 ppm matching, use all three ppm parameters together:

```text
--ms1-ppm 10 --feature-ppm 10 --fragment-ppm 10
```

`--ms1-ppm` controls MS1 library candidate lookup, `--feature-ppm` controls targeted mzML precursor assignment back to the merged targetlist feature, and `--fragment-ppm` controls MS2 fragment matching.

## RT Fitting

RT fitting is a post-processing filter. Series with too few points, unparsed lipid series, or RANSAC fitting failures are kept as `not_evaluated`; only rows classified as fitted outliers are removed from `targeted_mzml_rt_fitted_identifications.xlsx`.

## Notes

The MS2 fragment-generation and scoring logic still delegates to `ms2_legacy_core.py`. The targeted mzML workflow only changes data orchestration: feature merging, mzML reading, spectrum-to-feature mapping, and final long-form output.
