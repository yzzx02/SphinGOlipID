# Scientific logic correction

## Scope and source decisions

This update continues `codex/algorithm-submission-lock` from audit baseline
`4ceeae448bee9851c69ead63d323d7f5a8eb50a9`. It changes future production behavior.
It does not rerun samples, overwrite historical outputs, or regenerate the
1,723 identifications or 1,173 RT-related historical results. The theoretical
formula library and the disabled structural annotation resolution framework
are outside this change.

The author clarified that the final manuscript's displayed glycan notation
should be used and that historical impact counts should **not** be computed.
The final `中文初稿改1(1).docx`, Figure 4 and adjacent paragraph, explicitly gives
`Gal-Gal(-Fuc)-GlcNAc-Gal-Glc`: consecutive residues form the main chain and
parenthesized Fuc is a branch substituent. No verified `#` grammar was found.
Consequently `#` raises an explicit error instead of being interpreted speculatively.

## Glycan production path

`glycan_encoding.py` provides `parse_glycan_encoding()` and immutable normalized
main-chain/branch records. In final-paper notation the parenthesized residue
attaches to the immediately preceding main-chain residue (zero-based position
in the normalized record); sequence order is non-reducing end toward ceramide.
Supported syntax is `Residue-Residue(-BranchResidue)-Residue`, with one
single-residue branch per attachment. Nested branches, multi-residue branches,
multiple substituents at the same attachment, unknown residues, malformed
separators and `#` are rejected. These encodings do not specify linkage carbon
numbers or stereochemistry; those are not inferred.

`GSL_fragments()` passes its existing residue mass dictionary to the new loss
generator. Main-chain prefixes retain branches until their attachment residue
leaves; branch-first losses and their subsequent main-chain losses are also
generated. Residues are subtracted only once per path. The existing water
loss value and NeuAc diagnostic ions (292.10 and 274.09) are preserved.
Production tables deduplicate by exact theoretical m/z within each annotation
and retain alternative labels as explanations, not additional evidence counts.

Historical space-separated `-Residue` tokens and comma-separated integer
positions remain supported by the original GSL enumeration. Leading empty
tokens remain significant: the historical GM1 template is
` -Gal -GalNAc -NeuAc -Gal -Glc`, positions `0,3`. It normalizes to
`Gal-GalNAc-Gal(-NeuAc)-Glc`. The new representation preserves the legacy core
loss set, excluding the empty sentinel's unfragmented precursor/water pair,
and can add topology-consistent branch-retaining paths. It is not claimed that
the old incomplete enumeration and the new enumeration have identical full sets.

Targeted MassHunter CSV accepts optional `structure`/`glycan_encoding` and
`classy`/`branch_positions` columns, validates them and passes them into production.
Names exactly identifying GM1 and type I B antigen have documented topology
fallbacks. Other existing name fallbacks remain; they do not establish a
universal branched-isomer reconstruction capability. Supply explicit metadata
for other branched structures. Both Excel and CSV reach the same core parser.

## Independent centroid matching

`fragment_assignment.py` is shared by legacy text `process_sheet`, targeted
mzML matching/finalization and the normalized `scoring.match_fragments` API.
Within each candidate, admissible edges use the existing inclusive ppm window.
The deterministic assignment optimizes, in order:

1. Maximum number of observed–theoretical pairs.
2. Minimum sum of absolute ppm errors.
3. Maximum sum of observed intensities.

Residual augmenting paths can revise earlier pairings. Exact rational costs
avoid epsilon weights or new scientific thresholds. Equal theoretical masses
are one theory node even when several fragment labels explain them. Each
observed centroid and each theory mass can appear at most once per candidate.
Different candidates remain separate competing structural hypotheses.

The existing coverage and intensity-rank score equation is preserved, but its
inputs now count independent assignments; therefore scores and ranks can change.
Top-N ties are resolved consistently by score, coverage, intensity, then lexical
annotation. The text runner's low-level `query()` still constructs all edges;
one-to-one selection is applied before its output is finalized.

## Linear-first RT selection

`rt_iup.py` retains carbon-wise median RT, RANSAC initial inliers, iterative
residual filtering and final OLS refit. A valid positive-trend Linear fit with
sufficient distinct carbon numbers and R² ≥ 0.99 is accepted immediately.
Quadratic is evaluated by `fit_line()` only after Linear fails; it still needs
at least four distinct carbon numbers, R² ≥ 0.99 and the existing physical
monotonicity constraints. `quadratic_min_r2_gain` is removed.

Bracket rescue likewise tries Linear first and considers Quadratic only if
Linear fails fitting or the physical/bracket constraints. Among equally
complete feasible global IUP candidate sets, Linear count precedes support and
R² preferences; the existing bounded ordering/shift feasibility constraints
still apply. Existing ECN/IUP windows are unchanged. The separate older
`rt_validation.py` workflow is unchanged and should not be described as this
`rt_iup.py` implementation. RT remains an orthogonal validation step and does
not upgrade structural annotation resolution.

## Old-vs-new comparisons and impact limits

Reproduce the committed synthetic comparisons from the repository root with:

```sh
python scripts/scientific_logic_shadow.py
```

The script loads old functions directly from the fixed Git baseline. All
spectra, precursor masses and RT point sets are constructed in the script;
it does not open experimental data. The Figure 4 legacy integer counterpart
is a constructed encoding of the displayed structure, not an experimental row.
Every comparison row is labeled `synthetic_only`.

| Comparison | Output |
|---|---|
| Candidate count/intensity/score/rank/Top1 | `audit/fragment_matching_shadow_comparison.csv` |
| GM1 and Figure 4 branch fragment sets | `audit/glycan_branch_shadow_comparison.csv` |
| RT model choice and R² | `audit/rt_linear_first_shadow_comparison.csv` |

The synthetic duplicate-peak example demonstrates a possible Top1 change;
the curved RT example demonstrates Quadratic → Linear when Linear is already
valid. These are regression examples, **not historical impact statistics**.

| Requested historical impact | Status |
|---|---|
| Candidates changed by branch fragments | Not computed, per author instruction |
| Candidate score/rank changes | Not computed, per author instruction |
| Precursor Top1 changes | Not computed, per author instruction |
| RT series model-type changes | Not computed, per author instruction |
| Pass/fail changes among 1,173 RT-valid annotations | Not computed, per author instruction |

The historical files remain untouched; this does not establish that rerunning
the corrected code would reproduce them. A branch-only reassessment could be
limited to affected glycan candidates, but matching and RT selection apply
across classes. If a future release claims the final 1,723 annotations were
generated by this corrected implementation, validation should cover the full
original candidate search and RT inputs in a separate output location; checking
only the previously selected rows would miss changed competitors. No such rerun
was performed or scheduled here.

## Verification

Tests use artificial spectra and RT points. Coverage includes the final-paper
structure, historical GM1 core fragments, CSV metadata transmission, explicit
syntax errors, all possible three-by-three assignment edge graphs against an
exhaustive optimality oracle, inclusive ppm boundaries, intensity tie-breaking,
duplicate labels, deterministic Top-N, real Linear/Quadratic fits, lazy fallback,
and Linear-first rescue/IUP preference.

Final full-suite command: `python -m pytest -q` — **138 passed**, 30 pandas
`GroupBy.apply` deprecation warnings. The warning category already existed;
additional production-path tests exercise it more often. `git diff --check`
passed. No experimental dataset was used for this verification.
