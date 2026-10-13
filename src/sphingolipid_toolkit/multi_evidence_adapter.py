"""Opt-in adapters for existing production theory/observed tables.

No historical files are read or rewritten. Batch runners retain their existing
output contract. Call these adapters explicitly to obtain a separate evidence
table and evaluation; never reconstruct full spectra from final matched rows.
"""
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from .evidence import annotate_theoretical_fragments, capture_fragment_origins
from .multi_evidence_scoring import evaluate_candidate_evidence


def generate_evidence_fragments_for_candidates(candidates, observed_mz, ms1_rt, abundance, *, eligibility_overrides=None):
    """Use existing generation, capture source, append metadata in a new table.

    Temporary generation uses an isolated new directory, not the pipeline's
    out_DB.csv. Chemistry, fragment counts and masses are not changed.
    """
    from .targeted_mzml_pipeline import generate_legacy_fragments_for_candidates
    results = []
    for _, candidate in candidates.iterrows():
        with TemporaryDirectory(prefix="sphingolipid_evidence_") as temporary:
            with capture_fragment_origins() as origins:
                theory = generate_legacy_fragments_for_candidates(
                    pd.DataFrame([candidate]),observed_mz,ms1_rt,abundance,Path(temporary))
            if theory.empty:
                continue
            lipid_class = str(candidate["name"]).split("(",1)[0]
            results.append(annotate_theoretical_fragments(theory,lipid_class,origins,eligibility_overrides))
    return pd.concat(results,ignore_index=True) if results else pd.DataFrame()


def evaluate_legacy_candidate(lipid_class, annotation, theoretical, full_observed_spectrum, precursor,
    *, legacy_row=None, enabled=True, **kwargs):
    """Text A/B, targeted mz/intensity or normalized data use the same evaluator.

    ``enabled=False`` returns a copy of legacy_row with no appended fields.
    With True, return CandidateEvidenceResult (summary plus fragment evidence).
    The caller must supply the actual MS1 decision and complete single spectrum.
    """
    if not enabled:
        return dict(legacy_row) if legacy_row is not None else None
    observed = full_observed_spectrum.copy()
    for old_mz, old_intensity in (("A","B"),("relmz","intense"),("mz","intensity")):
        if old_mz in observed and old_intensity in observed:
            observed = observed.rename(columns={old_mz:"fragment_mz",old_intensity:"fragment_intensity"})
            break
    for scan_column in ("idx","scan_id"):
        if scan_column in observed and observed[scan_column].nunique() > 1:
            raise ValueError("Supply one complete MS/MS spectrum, not pooled scans")
    old = dict(legacy_row) if legacy_row is not None else {}
    result = evaluate_candidate_evidence(lipid_class,annotation,theoretical,observed,precursor,
        legacy_match_score=old.get("匹配度分数",old.get("legacy_match_score")),
        legacy_total_score=old.get("总分数",old.get("legacy_total_score")),**kwargs)
    result.summary = {**old, **result.summary}
    return result
