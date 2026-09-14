# Alignment with the SphinGOlipID manuscript logic

The project is organized around the manuscript's main technical storyline: two-dimensional LC separation plus software-assisted sphingolipid annotation.

## Module mapping

1. **Theoretical formula library / MS1 candidate retrieval** (`formula_library.py`, `ms2_pipeline.py`)
   - The production MS2 workflow currently consumes a prebuilt theoretical MS1 library for accurate-mass candidate retrieval.
   - `formula_library.py` is presently an interface placeholder; the exact manuscript library (84,240 entries) or a deterministic generator for that library must be deposited before the manuscript release is tagged.

2. **Rule-based MS/MS fragment library** (`ms2_legacy_core.py`)
   - Encodes diagnostic headgroup ions, neutral losses, LCB-related fragments, and glycan-sequence fragments.
   - Preserves the uploaded MS2 matching algorithm as the current stable core.

3. **MS2 matching and scoring** (`ms2_pipeline.py`, `targeted_mzml_pipeline.py`)
   - Parses measured MS/MS spectra and/or targeted mzML data.
   - Matches measured fragments to candidate-specific theoretical fragments using ppm windows.
   - Scores/ranks candidates using matched-fragment evidence and fragment-intensity information.
   - Uses `SphinGOlipIDConfig`, `io_utils.py`, `logging_utils.py`, and `scoring.py` as the stable project API around the converted legacy algorithm.

4. **RT correction and validation** (`rt_validation.py`, `rt_iup.py`)
   - Supports homologous-series extraction, robust RT modeling, candidate validation, and RT/IUP evaluation used in the six-fraction workflow.

5. **Result cleaning** (`result_cleaning.py`)
   - Supports duplicate annotation removal and final result cleaning/export.

6. **User interfaces** (`gui.py`, `cli.py`)
   - Provides a Tkinter desktop GUI and command-line entry points that call the backend modules rather than reimplementing the annotation logic.

## Recommended manuscript wording for the algorithm section

The code structure supports describing SphinGOlipID as a rule-based annotation workflow consisting of MS1 candidate retrieval from a theoretical sphingolipid formula library, candidate-specific in silico MS/MS fragment generation using sphingolipid class-specific fragmentation rules and glycan-sequence encoding, measured-theoretical fragment matching within a defined ppm tolerance, evidence-based candidate scoring/ranking, and retention-time-based orthogonal validation using homologous-series models.

## Release note

For publication reproducibility, manuscript claims should be tied to a versioned repository release that contains (or directly links to) the exact theoretical MS1 library and the manuscript-locked final annotation table. See `docs/RELEASE_READINESS.md`.
