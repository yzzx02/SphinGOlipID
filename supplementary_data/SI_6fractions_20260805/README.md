# Six-fraction SphinGOlipID supplementary information

This directory contains the compact supplementary-information release for the
six HILIC fractions processed on 2026-08-04/05.

## Scope

- Six fractions: `0-5`, `5-15`, `15-17.5`, `17.5-20.5`, `20.5-22`, and
  `22-30` min.
- Initial TOP3 universe: 9,693 unique candidate-pair rows and 1,546 unique MS1
  feature IDs.
- Classification modes: detailed chain, long-chain-base series, and total
  carbon/total unsaturation.
- RT rules: ECN and IUP.
- Locked thresholds: fitted R² >= 0.99, RT residual <= 0.20 min, and IUP
  allowance 0.05 min.

## Files

- `SphinGOlipID_SI_results_6fractions_20260805.xlsx`: compact five-sheet SI
  workbook.
- `SI01_raw_identifications.csv`: all 9,693 initial candidate-pair rows, reduced
  to 16 core identification columns.
- `SI02_filter_summary.csv`: the six formal filter definitions for each of the
  three classification modes.
- `SI03_rt_filtered_all_scores.csv`: retained ECN/IUP rows for all score tiers,
  reduced to 20 core columns.
- `SI04_fit_lines.csv`: fitted-line parameters and RT/IUP audit fields for the
  fixed formal `match_score >= 0.50` analysis.

## Score tiers in `SI03_rt_filtered_all_scores.csv`

- `high_score_ge_0.50`: the fixed formal results fitted from candidates with
  `match_score >= 0.50`.
- `low_score_lt_0.50`: additional RT-retained rows with `match_score < 0.50`,
  obtained by fitting the complete 9,693-row universe.

The `analysis_scope` column keeps these two calculations explicit. Low-score
rows supplement the SI table and do not overwrite the fixed formal high-score
results.

## Exclusions

Wide internal QA tables, superseded workbook versions, hundreds of individual
plot images, converted `.mzML` files, and original Agilent `.d` directories are
not included in this repository release.
