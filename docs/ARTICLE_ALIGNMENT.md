# Alignment with the SphinGOlipID manuscript logic

The project is organized around the manuscript's main technical storyline: two-dimensional LC separation plus software-assisted sphingolipid annotation.

## Module mapping

1. **Theoretical formula library** (`formula_library.py`)
   - Enumerates sphingolipid classes, LCB/FA chain ranges, hydroxylation states, unsaturation, headgroup, and glycan composition.
   - Feeds MS1-level candidate retrieval.

2. **Rule-based MS/MS fragment library** (`ms2_legacy_core.py`)
   - Encodes diagnostic headgroup ions, neutral losses, LCB-related fragments, and glycan-sequence fragments.
   - Preserves the uploaded MS2 matching algorithm as the current stable core.

3. **MS2 matching and scoring** (`ms2_pipeline.py`)
   - Parses measured MS/MS spectra.
   - Matches measured fragments to theoretical fragments using ppm windows.
   - Scores candidates by matched-fragment coverage and relative intensity ranking.
   - Uses `SphinGOlipIDConfig`, `io_utils.py`, `logging_utils.py`, and `scoring.py` as the stable v0.3 API layer around the converted legacy algorithm.

4. **RT correction and validation** (`rt_validation.py`)
   - Reserved for homologous series extraction and RANSAC-based RT model validation.

5. **Result cleaning** (`result_cleaning.py`)
   - Reserved for duplicate annotation removal and final export formatting.

## Recommended manuscript wording for the algorithm section

The code structure supports describing the algorithm as a rule-based annotation workflow consisting of: MS1 candidate retrieval from a theoretical formula library, in silico MS/MS fragment generation using sphingolipid class-specific fragmentation rules and glycan-sequence encoding, measured-theoretical fragment matching within a ppm tolerance window, evidence scoring based on matched-fragment coverage and fragment intensity, and RT-based validation using homologous-series correction models.
