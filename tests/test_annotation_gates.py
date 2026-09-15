import json
import pandas as pd
import pytest

from sphingolipid_toolkit.evidence import annotate_theoretical_fragments
from sphingolipid_toolkit.evidence_config import EvidenceScoringConfig
from sphingolipid_toolkit.multi_evidence_scoring import PrecursorEvidence, evaluate_candidate_evidence
from sphingolipid_toolkit.multi_evidence_adapter import evaluate_legacy_candidate, generate_evidence_fragments_for_candidates


def evaluate(lipid, labels, matched, *, config=EvidenceScoringConfig(), annotation=None, origins=None):
    masses = [100.+10*i for i in range(len(labels))]
    theory = annotate_theoretical_fragments(pd.DataFrame({"idf":labels,"mz":masses}),lipid,origins)
    observed = pd.DataFrame({"fragment_mz":[masses[i] for i in matched],"fragment_intensity":[100.]*len(matched)})
    return evaluate_candidate_evidence(lipid,annotation or lipid+"(d18:1/16:0)",theory,observed,
        PrecursorEvidence(True,0.,"[M+H]+"),config=config,total_composition=lipid+"(d34:1)")


@pytest.mark.parametrize("matched,level", [([0],"species"),([0,1],"species"),([0,1,2],"molecular_species"),([1,2,3,4],"candidate")])
def test_hg_and_lcb_independently_require_half(matched,level):
    labels = ["[phosphocholine+H]+","d18:1","d18:1+OH","d18:1+2OH","SM(d34:1)-H2O"]
    result = evaluate("SM",labels,matched)
    assert result.summary["annotation_level"] == level
    assert json.loads(result.summary["gate_summary"])["LCB"]["required"] == 2


@pytest.mark.parametrize("matched,level", [([0],"species"),([0,2],"species"),([0,2,3],"molecular_species"),([2,3,4,5],"candidate")])
def test_cer_author_dehydration_nl_is_required_in_addition_to_lcb(matched,level):
    labels = ["Cer(d34:1)-1H2O","Cer(d34:1)-2H2O","d18:1","d18:1+OH","d18:1+2OH","[C2H5NO+H]+"]
    result = evaluate("Cer",labels,matched)
    assert result.summary["annotation_level"] == level
    assert json.loads(result.summary["gate_summary"])["NL"]["required"] == 1
    if level == "molecular_species":
        assert result.summary["lcb_directly_supported"]
        assert result.summary["fa_assignment_mode"] == "inferred_from_total_composition_minus_lcb"
        assert result.summary["fa_directly_fragment_confirmed"] is False


def test_strong_common_cannot_rescue_hg_or_nl_gate():
    for lipid in ("SM","Cer"):
        result = evaluate(lipid,["d18:1","d18:1+OH","[C2H5NO+H]+"],[0,1,2])
        assert result.summary["annotation_level"] == "candidate"
        assert result.summary["S_supp"] > .9
        assert result.summary["confidence_score"] is None


def test_neutral_gsl_nl_only_species_and_nl_plus_lcb_molecular():
    candidates = pd.DataFrame({"name":["LacCer(d34:1)"],"classy":[""],"structure":["Gal-Glc"]})
    all_theory = generate_evidence_fragments_for_candidates(candidates,1000.,5.,1.)
    theory = all_theory[all_theory.anno.eq("LacCer(d18:1/16:0)")]
    nl_mass = theory[theory.fragment_type.eq("NL") & theory.scoring_eligible].mz.iloc[0]
    lcb = theory[theory.fragment_type.eq("LCB")].mz.tolist()[:2]
    for masses,level in [([nl_mass],"species"),([nl_mass]+lcb,"molecular_species")]:
        result = evaluate_candidate_evidence("LacCer","LacCer(d18:1/16:0)",theory,
            pd.DataFrame({"fragment_mz":masses,"fragment_intensity":[100.]*len(masses)}),PrecursorEvidence(True,0.,"[M+H]+"))
        assert result.summary["annotation_level"] == level


def test_failed_precursor_skips_msms_even_if_inputs_are_invalid():
    result = evaluate_candidate_evidence("SM","SM(d18:1/16:0)",None,None,PrecursorEvidence(False,None,None))
    assert result.summary["score_status"] == "precursor_failed"
    assert result.summary["annotation_level"] == "candidate"
    assert result.fragments.empty
    with pytest.raises(ValueError):
        PrecursorEvidence(True,None,None)


def test_minimum_normalized_intensity_is_inclusive_and_centralized():
    labels = ["[phosphocholine+H]+","[C2H5NO+H]+"]
    theory = annotate_theoretical_fragments(pd.DataFrame({"idf":labels,"mz":[100.,200.]}),"SM")
    observed = pd.DataFrame({"fragment_mz":[100.,200.],"fragment_intensity":[5.,100.]})
    for minimum,expected in [(.05,True),(.05001,False)]:
        result = evaluate_candidate_evidence("SM","SM(d18:1/16:0)",theory,observed,PrecursorEvidence(True,0.,"[M+H]+"),
            config=EvidenceScoringConfig(minimum_normalized_intensity=minimum))
        assert result.summary["hg_pass"] is expected


def test_single_chain_never_claims_inferred_fa():
    result = evaluate("So",["d18:1","d18:1+OH","d18:1+2OH"],[0,1],annotation="So(d18:1)")
    assert result.summary["annotation_level"] == "molecular_species"
    assert result.summary["fa_assignment_mode"] == "not_applicable_single_chain"
    assert not result.summary["fa_directly_fragment_confirmed"]
    assert result.summary["score_status"] == "UNSPECIFIED_policy"


def test_actual_so_generator_does_not_invent_lcb_evidence():
    table = generate_evidence_fragments_for_candidates(pd.DataFrame({"name":["So(d18:1)"],"classy":["So"],"structure":[""]}),300.2898,5.,1.)
    assert not table.fragment_type.eq("LCB").any()


def test_wrong_lcb_identity_cannot_upgrade_candidate():
    result = evaluate("SM",["[phosphocholine+H]+","d16:1","d16:1+OH","d18:1"],[0,1,2,3])
    assert result.summary["annotation_level"] == "species"
    assert not result.summary["lcb_pass"]


def test_disabled_adapter_returns_exact_legacy_fields_without_evaluation():
    old = {"匹配度分数":.8,"总分数":80.,"相对强度":.5,"注释":"test"}
    result = evaluate_legacy_candidate(None,None,None,None,None,legacy_row=old,enabled=False)
    assert result == old and result is not old


def test_both_legacy_formats_preserve_old_score_columns_and_share_semantics():
    old = {"匹配度分数":.8,"总分数":80.,"相对强度":.5}
    theory = pd.DataFrame({"idf":["[phosphocholine+H]+"],"mz":[184.0733]})
    results = []
    for columns in (("A","B"),("relmz","intense")):
        observed = pd.DataFrame({columns[0]:[184.0733],columns[1]:[100.]})
        result = evaluate_legacy_candidate("SM","SM(d18:1/16:0)",theory,observed,PrecursorEvidence(True,0.,"[M+H]+"),legacy_row=old)
        assert all(result.summary[key] == value for key,value in old.items())
        assert result.summary["legacy_match_score"] == .8
        assert result.summary["legacy_total_score"] == 80.
        results.append(result.summary)
    assert results[0] == results[1]
