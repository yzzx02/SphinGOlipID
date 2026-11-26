# Changelog

## 1.0.0 — 2026-11-26

- Review measured peaks, MS/MS spectra and ECN fits from the desktop results tab. Export an annotation table with a measured peak area for each sample.
- Fit IOLS and RANSAC independently. Select the qualifying linear model with the most recovered inliers before evaluating quadratic fits. The adjustable absolute RT window defaults to 0.20 min.
- Update positive-mode fragment pools and the curated chain references. The bundled indexes contain 753,564 positive and 348,895 negative records.
- Preserve boolean and unjudged IUP results when using pandas 3.

## 0.6.3 — 2026-11-24

- Admit positive GM2 and GM3 when independent peaks match at least half of the combined core HG pool. Chain-level reporting requires the composition-dependent LCB evidence.

## 0.6.2 — 2026-11-24

- Limit each positive precursor to eight core HG channels, retain auxiliary products in Other, and add the sum-composition ammonium entry for GM2(36:1;O2).

## 0.6.1 — 2026-11-15

- Allow positive LCB compositions with up to three double bonds at any carbon
  length. Restore the previously curated SL triene records that meet the current
  acyl limits and include the explicit SHexSPB(d22:3) library candidate.

## 0.6.0 — 2026-11-10

- Retain the existing desktop layout and replace mandatory TXT/target-table
  inputs with multiple mzML files and one output XLSX; default to XCMS centWave.
- Expose peak-width limits, noise, S/N, positive/negative polarity and ppm/Da
  MS2 tolerance. Process files serially without alignment or gap filling.
- Retain native local peak IDs, areas, measured EIC apices, per-spectrum
  evidence and native ranking. Export features, annotations, all spectrum
  matches, association audit, exact settings and actual dependency versions.
- Add the independent `sphingolipid-untargeted` command and packaged R runner;
  original scientific libraries, pool rules, gates and score calculations remain
  unchanged.

- Use the displayed 0–100 total-score scale for GUI result filtering, including
  current total-score-only workbooks and read-only aliases for older results.
- Save complete result bundles with external MS2 evidence and all worksheets;
  retain existing files when evidence is missing, conflicting or a save fails.
- Preserve the existing Tk interface and all scientific rules, libraries,
  gates, candidate ranking and score values.

## 0.5.1 — 2026-11-10

- Add only the native-CID-verified d18:1 DMS composition with a two-HG gate including 110.0964, reusing HG/Other 60/40.
- Correct LM1 reference formulas/masses and trimethylSo intrinsic ion masses; keep seven rare types reference-only.
- Freshly rebuild and verify independent indexes with stable IDs and all previous scientific fields preserved; retain original libraries and migration checkpoint.

## 0.5.0 — 2026-11-10

- All identification entry points use the current migrated sphingolipid engine,
  including the historical MS2PipelineConfig constructor. Retired backend
  selectors and coverage/intensity/count parameters now fail explicitly.
- Retire the notebook generators and experimental evidence classifier/scorer;
  preserve archived result data containers and the existing SphinGOlipID frontend.
- Preserve complete external evidence references when combining text results.
- Document the eight reference-only historical types without enabling new rules.
- Scientific engine, fragment libraries, gates, weights and confirmation unchanged.

## 0.4.0 — 2026-11-10

- Import the attributed independent scientific-engine snapshot, sphingolipid
  libraries, type/name crosswalk and production adapters; retain the original UI.
