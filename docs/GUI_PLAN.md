# GUI construction plan for SphinGOlipID

## 1. Recommended technology route

For the first usable version, build the GUI with **Streamlit** or **PySide6**:

- **Streamlit** is fastest for lab-internal use: upload files, set parameters, run, preview tables, download results.
- **PySide6** is better for a polished standalone desktop program: local file selectors, progress bars, tabs, and packaged `.exe`.

Recommended development order: Streamlit prototype first, then PySide6 desktop packaging after the workflow is stable.

## 2. GUI page structure

### Page A. Project setup

Inputs:

- raw MS/MS text folder;
- target precursor Excel folder;
- theoretical MS1 database file;
- output folder;
- processing range, e.g. 1-6.

Validation:

- check whether `hilic-msms-{i}.txt` and `hilic-msms-{i}.xlsx` exist;
- check required columns: `target`, `下限`, `上限`, `RT`, `Abund`;
- check MS1 database columns: `理论值`, `classy`, `name`, `structure`.

### Page B. MS2 matching parameters

Controls:

- fragment tolerance, ppm;
- minimum measured fragment intensity;
- minimum matched fragment number;
- minimum MS2 match score;
- top-N candidates per precursor;
- keep/delete intermediate files.

Preset buttons:

- manuscript/default mode: 20 ppm, intensity ≥ 20, matched fragments ≥ 2, score ≥ 0.35, top 3;
- strict mode: lower ppm and higher score;
- exploratory mode: relaxed score and keep intermediates.

### Page C. Run monitor

Functions:

- show current file index;
- show current stage: spectrum parsing → candidate retrieval → theoretical fragment generation → fragment matching → scoring → export;
- show warning/error messages without stopping the whole batch unless a required file is missing;
- provide a run log.

### Page D. Results review

Tables:

- final candidates per precursor;
- matched fragments, measured m/z, intensity, theoretical fragment annotation;
- scores: match coverage, total intensity, total score, relative intensity.

Filters:

- lipid class;
- score threshold;
- RT range;
- precursor m/z;
- top candidate only / top 3 / all retained candidates.

Visualizations:

- bar plot of candidate scores per precursor;
- matched fragment mirror/stick plot;
- class distribution summary;
- later: RT predicted-vs-observed scatter plot and deviation histogram.

### Page E. Library/RT modules, later versions

After other modules are added:

- formula-library builder tab;
- glycan-sequence editor/checker;
- RANSAC RT model fitting tab;
- RT correction report export;
- final annotation report generator.

## 3. Backend API design

The GUI should not directly call notebook-style functions. It should call stable backend functions:

```python
from sphingolipid_toolkit import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import run_batch

config = SphinGOlipIDConfig(
    raw_ms2_dir=raw_dir,
    precursor_dir=target_dir,
    output_dir=out_dir,
    ms1_library_path=ms1_db,
    file_indices=(1, 2, 3, 4, 5, 6),
    fragment_ppm=20,
    min_fragment_intensity=20,
    min_matched_fragments=2,
    min_match_score=0.35,
    top_n=3,
)
result_files = run_batch(config)
```

This keeps the GUI thin and makes the same algorithm available from CLI, GUI, and future manuscript reproducibility scripts.

## 4. Development milestones

### Milestone 1: runnable MS2 prototype

- file selection;
- parameter panel;
- run button;
- output result Excel;
- basic log display.

### Milestone 2: result inspection

- interactive final result table;
- matched-fragment detail view;
- export filtered results;
- simple score plots.

### Milestone 3: full SphinGOlipID workflow

- formula-library generation;
- MS1 candidate retrieval;
- MS2 matching;
- RT correction/validation;
- duplicate removal;
- final annotation export.

### Milestone 4: desktop packaging

- PyInstaller packaging;
- bundled default config;
- user manual;
- example dataset and tutorial;
- versioned release notes.

## 5. Suggested GUI folder to add later

```text
gui/
├── streamlit_app.py
├── pages/
│   ├── 1_MS2_Matching.py
│   ├── 2_Result_Review.py
│   ├── 3_RT_Validation.py
│   └── 4_Library_Builder.py
└── assets/
    └── sphingolipid_logo.png
```
