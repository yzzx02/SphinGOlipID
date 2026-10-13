# Alignment with the SphinGOlipID manuscript logic

The project is organized around the manuscript's main technical storyline: two-dimensional LC separation plus software-assisted sphingolipid annotation.

## Module mapping

1. **Theoretical formula library** (`formula_library.py`)
   - Enumerates sphingolipid classes, LCB/FA chain ranges, hydroxylation states, unsaturation, headgroup, and glycan composition.
   - Feeds MS1-level candidate retrieval.

2. **Rule-based MS/MS fragment library** (`ms2_legacy_core.py`)
   - Encodes diagnostic headgroup ions, neutral losses, LCB-related fragments, and glycan-sequence fragments.
   - Retains class-specific residue masses and rules, with the scientific corrections documented below.
   - Final-paper parentheses notation is parsed into main chain and attached branches. Legacy space-separated tokens and integer `classy` indices remain supported. CSV preserves explicit metadata and has documented GM1/type I B topology fallbacks. Unverified `#` syntax is rejected.

3. **MS2 matching and scoring** (`ms2_pipeline.py`)
   - Parses measured MS/MS spectra.
   - Matches measured fragments to theoretical fragments using ppm windows.
   - Scores candidates by matched-fragment coverage and relative intensity ranking.
   - Uses `SphinGOlipIDConfig`, `io_utils.py`, `logging_utils.py`, and `scoring.py` as the stable v0.3 API layer around the converted legacy algorithm.
   - Uses shared one-to-one assignment: maximum pair count, minimum total absolute ppm error, then intensity. Each centroid and each theoretical m/z contributes once per candidate.

4. **RT correction and validation** (`rt_iup.py`, `rt_validation.py`)
   - `rt_iup.py`: carbon-wise median RT → RANSAC inlier seed → iterative residual pruning → final OLS refit → Linear-first evaluation with Quadratic fallback → ECN/IUP validation. It does not compare independent iterative-OLS and RANSAC model families.
   - `rt_validation.py` is a separate older RT implementation used by the targeted batch RT filter; calling that filter is not equivalent to running the `rt_iup.py` ECN/IUP workflow.

5. **Result cleaning** (`result_cleaning.py`)
   - Provides deterministic score/intensity-prioritized RT near-duplicate removal. This utility is not automatically invoked by either MS2 batch runner.

## Recommended manuscript wording for the algorithm section

Describe the specific executed workflow: candidate retrieval, legacy class-specific fragment generation, configured ppm matching, matched-theoretical-fragment coverage and intensity ranking, followed by the RT workflow actually run. For `rt_iup.py`, use “RANSAC-seeded OLS refit”; describe the supplied glycan topology and independent-centroid matching. Do not claim support for `#` syntax or universal branch reconstruction from lipid names.

See [the submission algorithm audit](ALGORITHM_AUDIT_20260914.md) for behavior changes, regression coverage, result-impact limits, and unresolved manuscript/code gaps. This audit did not review or change the theoretical formula library.

Current corrections and synthetic comparisons: [scientific logic report](SCIENTIFIC_LOGIC_UPDATE_REPORT.md). The preceding audit is a historical snapshot.
