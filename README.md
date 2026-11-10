# SphinGOlipID

Sphingolipid identification with the current fragment-pool engine, an independent
curated library and the existing SphinGOlipID desktop frontend. LipidGate's
frontend is not imported. The attributed scientific-engine snapshot runs without
a LipidGate installation, pyOpenMS or a Qt application.

## Install

```powershell
git lfs pull
python -m pip install -e .
```

Python 3.10 or later is required. Compressed MSP files, SQLite indexes and the
historical reference-name catalog are tracked by Git LFS and bundled in the
package. An incomplete or fingerprint-mismatched bundled index is rejected.

The non-targeted desktop workflow also requires R with `xcms`, `MSnbase`,
`mzR`, `jsonlite` and `BiocParallel`. Add `Rscript` to PATH; Windows installations
under `Program Files/R` are also detected. There is no silent detector fallback.
XCMS 4.8.0 / R 4.5.2 is the locally validated combination.

## Identify

```powershell
sphingolipid-gui
sphingolipid-untargeted sample1.mzML sample2.mzML --output results.xlsx
sphingolipid-identify sample.mzML --output identifications.xlsx --fragment-ppm 20
sphingolipid-ms2 --input-dir RAW --precursor-dir TARGETS --output-dir RESULTS --files 1-6
sphingolipid-targeted-mzml --data-dir DATA --output-dir RESULTS --fragment-da 0.3
```

The original Tkinter setup and review layout is retained. Select one or more
centroided mzML files containing MS1/MS2 and an output XLSX path; no target table
is required. XCMS centWave extracts each file's actual local MS1 chromatographic
peaks serially. Peak-width limits (seconds), noise and S/N are editable, along
with MS1 tolerance, MS2 tolerance (ppm or Da), minimum score and Top-N. The MS1
tolerance also supplies centWave's consecutive-scan ROI ppm setting.
The GUI uses the same current MS2 engine as the text and mzML commands.
The separate historical text command still uses paired
`hilic-msms-{i}.txt` and `hilic-msms-{i}.xlsx` target tables; the two folders can
be supplied separately. Targeted mzML input uses mzML files and targetlist CSVs.
`sphingolipid-identify` can identify native mzML directly without those targetlists.
Historical MS1 workbooks do not provide identification candidates; no such
workbook is required for current matching.

| Setting | Default |
| --- | --- |
| Precursor tolerance | 10 ppm |
| Fragment tolerance, text / direct mzML | 20 ppm |
| Fragment tolerance, targeted mzML | 10 ppm |
| Absolute fragment tolerance | `--fragment-da` replaces ppm |
| Minimum total score | 50 / 100 |
| Candidates requested | 3; existing score/count ranking and ties retained |

For the desktop / `sphingolipid-untargeted` workflow, Top-N defaults to **1**.
Its centWave settings default to peakwidth 5–60 s, noise 1000, S/N 3,
prefilter `(3, noise)`, MS1/ROI tolerance 10 ppm and mzdiff −0.001.
These are application defaults; they do not claim to be the upstream XCMS
defaults. Positive mode is selected initially; negative mode can be selected.
See [centWave parameter documentation](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.html).

No cross-file alignment or gap filling is performed. Local XCMS peak IDs,
measured boundaries and areas are retained. Apex RT, intensity and centroid m/z
are measured from the MS1 EIC inside those boundaries. An acquired MS2 needs a
confirmed measured preceding MS1 and one unambiguous local peak of the same
file and polarity. Missing detection never creates an MS2-only feature.
Each feature retains the native candidate order from one best acquired
spectrum, chosen by its first-ranked score, matched-fragment count, then
acquisition order. Fragments from different scans are never combined.

The defaults above preserve existing workflow choices. Actual tolerances are
recorded in results. Use `--top-n 1` to request Top1 under the current ranking
policy; equal scores do not override the matched-fragment-count split.

## One current backend

Version 0.5.0 retired the notebook/coverage backend and the separate experimental
multi-evidence scorer. Every formal entry point, including the compatibility
constructor `MS2PipelineConfig`, now invokes the same current engine.
`backend="legacy"` and old absolute-intensity, global matched-count and 0–1
coverage-score parameters raise an explicit migration error instead of being
accepted or silently ignored. Their CLI flags have been removed.
The internal backend name `lipidgate` identifies the attributed scientific
snapshot, not a dependency on the source application or its GUI.

The archived evidence vocabulary remains readable for historical tables. Its
old classifier, generators, gates and scorers cannot annotate new spectra.
Historical research outputs are not automatically rescored.

Scientific fragments, pools, admission, weights and chain confirmation are
unchanged by this interface cleanup. Other never imposes admission and counts
raw experimental support strictly above 0.5% of the base peak. Shared d/t chain
confirmation uses composition-specific LCB-H2O / LCB-2H2O at inclusive 1%.
The source's candidate ordering, scores, matched-count Top1 split and full
evidence provenance are retained. See
[migration and scientific coverage](docs/LIPIDGATE_ENGINE_MIGRATION.md).

## Python API

```python
from sphingolipid_toolkit import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import run_batch
from sphingolipid_toolkit.engine import IdentificationConfig, SphingolipidIdentifier

config = SphinGOlipIDConfig.from_single_input_dir(
    "./data/raw_ms2", "./outputs/ms2_results", file_indices=(1,),
    fragment_ppm=20, min_total_score=50, top_n=1,
)
result_files = run_batch(config)

with SphingolipidIdentifier(IdentificationConfig(fragment_ppm=20, top_n=1)) as identifier:
    rows = identifier.identify(538.5199, [(264.2686, 100.0)], polarity="+")
```

Positive and negative modes use their separate current libraries. New library
indexes must be built through `scripts/build_sphingolipid_library.py` with fresh
normalization and verification; changing a stored rules hash is insufficient.

## Results and reference names

Outputs contain reference-engine pool scores, matched-fragment counts, rank,
reporting resolution and `ms2_evidence_json`. Compatibility columns support the
existing review and RT tools without rescoring candidates. Excel cells exceeding
the format limit retain complete evidence in an adjacent content-addressed JSON
file; aggregate outputs preserve relative evidence-file references.

The non-targeted workbook contains `annotations` (one selected acquisition per
local feature, with its Top-N candidates), `features` (all detected local peaks,
including unannotated peaks), `spectrum_matches` (every passing acquired-spectrum
result), `spectrum_audit` (local-association decisions), `parameters` and
`input_files` (actual R/XCMS versions and source paths). Scores, gates, native
scan IDs, ranking and full MS2 evidence remain unchanged. Keep adjacent
`ms2-evidence-*.json.gz` files with the workbook; GUI Save As copies the bundle.

The historical 85,710-row MS1 list is a type/name crosswalk, not a library-size
target or an identification candidate grid. Current eligible chain combinations
remain authoritative. Of 112 retained sphingolipid type names, 105 have current
MS2 class/family coverage (DMS only at d18:1); seven remain reference-only. This is type availability,
not measured recall, exact-isomer resolution or coverage of every historical
composition. Myriocin is archived separately as non-sphingolipid.
See [the complete type audit](docs/sphingolipid_type_coverage.md) and
[the eight-type literature review](docs/sphingolipid_literature_review_20261110.md).

## Retained RT and review tools

`rt_validation.py`, `rt_iup.py` and `result_cleaning.py` retain the project's
existing RT modeling, review and cleaning workflows. RT validation is orthogonal
to MS2 evidence and does not replace pool admission or chain confirmation.
Archived notebook identification wrappers fail explicitly; RT wrappers remain
available. Historical documents describing notebook or experimental scoring
are development history rather than current identification settings.

## Development

```powershell
python -m pytest -q
python -m compileall src tests
```

Continue on `codex/lipidgate-sphingolipid-migration`; do not merge into main
without an explicit user request. Preserve the original HG2/B2 reference and
published checkpoints. Scientific code and library revisions are independently
traceable. Raw samples, research outputs, caches and builds stay untracked.


## Literature revision 2026-11-10

Version 0.5.1 adds one native-CID-supported N,N-dimethylSo(d18:1) composition,
requiring two N-related HG products including 110.0964 with HG/Other 60/40.
The independent reference catalog corrects LM1's missing NeuAc residue and
trimethylSo's intrinsic [M]+ ion; these remain reference-only, without a direct
GM1 alias or a permissive single-amine gate. No other rare-type grid is added.
Existing classes retain their fragments, evidence rules, scores and ranking.
The local scientific revision is separate from the immutable migration manifest.
