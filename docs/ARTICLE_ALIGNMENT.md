# Alignment with the SphinGOlipID manuscript logic

The project is organized around the manuscript's main technical storyline: two-dimensional LC separation plus software-assisted sphingolipid annotation.

## Module mapping

1. **Theoretical formula library** (`formula_library.py`)
   - Enumerates sphingolipid classes, LCB/FA chain ranges, hydroxylation states, unsaturation, headgroup, and glycan composition.
   - Feeds MS1-level candidate retrieval.

2. **Rule-based MS/MS fragment library** (`ms2_legacy_core.py`)
   - Encodes diagnostic headgroup ions, neutral losses, LCB-related fragments, and glycan-sequence fragments.
   - Preserves the uploaded MS2 matching algorithm as the current stable core.
   - Legacy glycan input uses space-separated residue-loss tokens plus comma-separated position indices in `classy`. Leading spaces are significant. This is not a `#`-encoded glycan-tree parser; the name-only CSV fallback does not supply branch indices.

3. **MS2 matching and scoring** (`ms2_pipeline.py`)
   - Parses measured MS/MS spectra.
   - Matches measured fragments to theoretical fragments using ppm windows.
   - Scores candidates by matched-fragment coverage and relative intensity ranking.
   - Uses `SphinGOlipIDConfig`, `io_utils.py`, `logging_utils.py`, and `scoring.py` as the stable v0.3 API layer around the converted legacy algorithm.
   - Counts matched theoretical m/z values after choosing the strongest centroid per theoretical m/z. One observed peak may support multiple nearby theoretical m/z values; the count is not necessarily the number of independent observed peaks.

4. **RT correction and validation** (`rt_iup.py`, `rt_validation.py`)
   - `rt_iup.py`: carbon-wise median RT → RANSAC inlier seed → iterative residual pruning → final OLS refit → linear/quadratic evaluation → ECN/IUP validation. It does not compare independent iterative-OLS and RANSAC model families.
   - `rt_validation.py` is a separate older RT implementation used by the targeted batch RT filter; calling that filter is not equivalent to running the `rt_iup.py` ECN/IUP workflow.

5. **Result cleaning** (`result_cleaning.py`)
   - Provides deterministic score/intensity-prioritized RT near-duplicate removal. This utility is not automatically invoked by either MS2 batch runner.

## Recommended manuscript wording for the algorithm section

Describe the specific executed workflow: candidate retrieval, legacy class-specific fragment generation, configured ppm matching, matched-theoretical-fragment coverage and intensity ranking, followed by the RT workflow actually run. For `rt_iup.py`, use “RANSAC-seeded OLS refit”; report branch-index metadata only when supplied. Do not claim support for `#` syntax, universal branch reconstruction from lipid names, or independent-peak counting.

See [the submission algorithm audit](ALGORITHM_AUDIT_20260914.md) for behavior changes, regression coverage, result-impact limits, and unresolved manuscript/code gaps. This audit did not review or change the theoretical formula library.
