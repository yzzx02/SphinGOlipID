# Independent Figure 4 annotation framework

This is an architectural extension only. Existing sample data and historical
identifications are not reprocessed, recounted, filtered, ranked or rewritten.
Neither existing pipeline imports this framework. No legacy output columns,
fragment generation, matching, scoring, Top-N selection or RT behavior change.

## Workflow and terminology

The immutable `FIGURE4_WORKFLOW` metadata is exactly:

`MS1 Candidates → MS/MS Targets → Species Level → Molecular Species → RT Validation`

The formal annotation levels are `Species level` and `Molecular species level`.
An unresolved annotation has `level=None`; `UNSPECIFIED` describes missing rules,
not an additional structural level. MS1 candidates and matching are upstream
stages; RT validation is a downstream, independent orthogonal check.

`annotation_evidence.py` separates chemical `EvidenceType` from potential
`EvidenceRole`. HG/GLYCAN suggest `SPECIES_EVIDENCE`; LCB/FA suggest
`MOLECULAR_SPECIES_EVIDENCE`. These roles do not assert diagnostic sufficiency.
PRECURSOR, SUPPORT, unclassified fragments and ordinary neutral losses default
to `SUPPORTING_EVIDENCE`. Supporting evidence remains available for future
confidence assessment; this framework does not implement or change scoring.

## Rule contract

`RuleRegistry` uses exact subclass keys. No scientific subclass rules or numeric
thresholds ship with this extension. An absent rule returns `UNSPECIFIED` with
a TODO explanation, even if an HG or LCB fragment matched. Different subclasses
can register different callbacks and record their scientific provenance.
Synthetic rules in tests are not proposed scientific rules.

Each `EvidenceSummary` belongs to one candidate, subclass and spectrum. Callers
must provide actual matching provenance; theoretical presence alone is not a
match (`MatchedEvidence.matched` defaults to false). Species rules receive only
matched Species evidence. Future subclass-specific necessary evidence can be
explicitly assigned that role with documented provenance; the adapter does not
infer such assignments. Molecular rules receive only matched molecular-species
evidence with a recorded chain composition, and run only after the same summary
satisfies Species level. A single matched LCB never guarantees promotion.
Isobaric ambiguity and necessary evidence combinations remain TODO/UNSPECIFIED.

Rules return `EvidenceEvaluation` with SATISFIED, NOT_SATISFIED or UNSPECIFIED.
The framework enforces prerequisites; subclass rules remain responsible for
scientifically justified combinations, applicability and chain specificity.
No RT result is passed to a rule. `attach_rt_validation()` attaches an externally
supplied RT status without recomputing RT or changing any annotation level.

## Legacy metadata inventory

All locations below refer to `ms2_legacy_core.py` unless stated otherwise.

| Source | Metadata and limitation |
| --- | --- |
| `LCB_fragment`, `LCB_mz`, `almost_fuc()` | Paired LCB ID/mass entries, including OH variants, are potential LCB evidence. `lcb_metadata()` reads supplied mappings without modifying them. Chain ambiguity/sufficiency is unspecified. |
| `SM()`, `Lyso_SM()` | Named phosphocholine ions are HG. `P-O断裂` remains unclassified. |
| `PE_cer()`, `PG_cer()`, `CAEP()`, `N_CAEP()` | Named headgroup-derived ions are HG candidates, not proof of subclass specificity. Original malformed N_CAEP closing bracket is preserved. |
| `CerP()`, `S1P()`, `Lyso_sulfo()` | Phosphate/sulfate-related ions are HG candidates; subclass sufficiency is unspecified. |
| `PI_cer()` | `[phosphate+H]+` naming/mass interpretation requires verification; remains unclassified. |
| `GSL_fragments()` | NeuAc ions are GLYCAN. Sequential/branch glycan losses can use NEUTRAL_LOSS plus GLYCAN loss origin. Empty glycan tokens/zero losses must not be treated as diagnostic glycans. |
| `FMC_1/3/5()`, `EO_Glccer()`, `Glu_So()`, `Gb3_So()` | Glycan/acetyl loss provenance hooks only; no automatic subclass rule. Preserve Glu_So ID spelling. |
| `learn_fuc()` and subclass headgroup-loss branches | Ordinary water and headgroup losses: NEUTRAL_LOSS, supporting by default. |
| `EO_cer()`, `EO_Glccer()`, `OCer()` | Acyl-related neutral losses, not established independent FA ions. Chain specificity and sufficiency remain TODO. |
| `almost_fuc()` FA_C/FA_U | Composition subtraction is not matched FA evidence. No FA classifier is implemented. |
| `So()`, `O_FA_cer()` | Empty bodies do not establish diagnostic requirements; OCer is dispatched separately. Cer/lyso requirements remain unspecified. |
| Pipeline initialization `[C2H5NO+H]+` | Common SUPPORT fragment, never automatic chain evidence. |
| MS1 target/precursor | PRECURSOR context, not structural MS/MS evidence. |

`lookup_legacy_metadata()` uses exact source/ID pairs for the conservative named
ion subset. Unknown IDs use `evidence_type=None`, supporting role and a TODO
note. It does not guess chemistry from m/z or lipid-name substrings.
`neutral_loss_metadata()` is an explicit provenance hook, not a loss parser.
Original IDs are never corrected in legacy tables. FA remains a vocabulary and
future metadata extension point only; no independent legacy FA evidence is claimed.

## Isolation and future integration points

`AnnotationResolutionConfig.annotation_resolution_enabled` defaults to false.
`annotate_sidecar()` returns the identical legacy object plus `None` when disabled,
without inspecting evidence. When explicitly enabled by a caller it returns a
separate resolution object and still preserves the legacy object. The pure
resolver can also be called directly for tests or future integrations.

Potential future boundaries are theoretical metadata lookup after generation,
and independent evidence interpretation after `process_sheet()` or
`match_observed_to_theoretical()`. No hooks are installed in either pipeline.
`_finalize_results()`, `finalize_spectrum_matches()`, `process_group()`,
`rt_iup.py` and all existing result writers remain untouched. Enabling the
standalone config does not activate a production pipeline feature.

## Safe verification

Run only the new synthetic tests:

```text
python -m pytest tests/test_annotation_evidence.py tests/test_annotation_resolution.py tests/test_annotation_legacy_compatibility.py
```

Tests use in-memory artificial evidence/tables and source inspection. They do
not invoke experimental files, sample pipelines, fragment generation or RT
fitting. Compatibility coverage verifies identical output object, values,
columns, ordering and scores at the sidecar boundary, and absent imports from
existing execution modules. This is not a claim of historical-data replay.
