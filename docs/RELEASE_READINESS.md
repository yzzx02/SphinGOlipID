# Release readiness for manuscript-linked GitHub repository

This checklist records the items that should be complete before the repository URL is cited in the manuscript or Supporting Information.

## Already present

- Installable `src/`-layout Python package.
- CLI entry points for legacy TXT-based MS2 matching, targeted mzML processing, and the Tkinter GUI.
- Unified configuration, input validation, logging, scoring helpers, result cleaning, RT validation, and RT/IUP modeling modules.
- Unit/integration-style tests with compact mock data.
- User, input/output, GUI, developer, targeted-mzML, and manuscript-alignment documentation.
- Six-fraction RT/SI result tables under `supplementary_data/SI_6fractions_20260805/`.
- Cross-version GitHub Actions CI for Python 3.10–3.12 on the release-readiness branch.

## Must be resolved before the manuscript links to this repository

### 1. Publish the actual theoretical MS1 library used in the manuscript

The manuscript describes a theoretical library containing 84,240 sphingolipid entries. The repository currently contains only a mock MS1 library for tests, while `formula_library.py` is still a placeholder. Before release, do **one** of the following:

1. add the exact XLSX/CSV theoretical library used for the manuscript, or
2. implement and document a deterministic library-generation workflow that reproduces the deposited library.

The deposited file should have a stable, descriptive name and should be referenced directly from the README.

### 2. Publish the manuscript-locked final annotation table

Add the final annotation table corresponding to the numbers and figures in the submitted manuscript. Keep an immutable copy with a version/date in the file name, and document the total number of annotations and the relationship to the RT-validated subset.

### 3. Reconcile MS/MS mass-tolerance documentation

The stable TXT workflow currently defaults to a fragment tolerance of 20 ppm, whereas the targeted mzML workflow defaults to 10 ppm. The manuscript must report the parameter actually used for the final analysis, and the README/config examples should make that manuscript setting explicit rather than implying that all workflows share the same default.

### 4. Add a software license

No license should be inferred automatically. Choose the intended license (for example MIT, BSD-3-Clause, GPL-3.0, or another institution-approved license) and add a root-level `LICENSE` file before public citation.

### 5. Add citation metadata after author order is final

Add a root-level `CITATION.cff` once the final software author/contributor order and manuscript citation are known. Include the repository URL and, after archiving a release, the DOI if one is minted (for example through Zenodo).

### 6. Create a versioned manuscript release

Before submission or revision, create a stable release/tag (recommended pattern: `v1.0.0-manuscript` or `v1.0.0`) from the exact commit that matches the manuscript. The manuscript/SI should link to the repository and preferably cite the archived release DOI rather than only the moving `main` branch.

## Strongly recommended

- Replace internal-development wording such as “project scaffold” and “reserved for later” where it no longer reflects the current code.
- Keep `plan.txt` as development history, but do not present it as current user documentation.
- Add one compact end-to-end example that produces a non-empty annotation result from distributable mock/example data.
- Add tests for result cleaning and the manuscript-used RT-validation path if those modules are part of the claimed software workflow.
- Record the exact environment used for manuscript analysis (Python version and dependency versions) in a release/environment file.
- Keep raw proprietary/institutional data out of the public repository unless sharing is explicitly permitted.

## Manuscript-to-code consistency checks

Before tagging the release, verify that the repository and manuscript agree on:

- MS1 mass tolerance;
- MS/MS fragment tolerance;
- minimum fragment intensity and minimum matched-fragment requirement;
- candidate ranking/scoring definitions;
- Top-N candidate retention;
- RT normalization and fraction correction;
- iterative OLS/RANSAC selection rules and thresholds;
- final annotation counts shown in the manuscript figures/tables.
