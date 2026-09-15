import pandas as pd
import pytest

from sphingolipid_toolkit.evidence import annotate_theoretical_fragments
from sphingolipid_toolkit.evidence_config import EvidenceScoringConfig
from sphingolipid_toolkit.multi_evidence_scoring import PrecursorEvidence, evaluate_candidate_evidence, add_rt_validation


def synthetic_sm():
    # Artificial masses distinguish scoring mathematics from chemistry tests.
    return annotate_theoretical_fragments(pd.DataFrame({
        "idf":["[phosphocholine+H]+","d18:1","d18:1+OH","d18:1+2OH","SM(d34:1)-C5H14NO4P","SM(d34:1)-H2O"],
        "mz":[100.,200.,210.,220.,300.,400.]}),"SM")


def evaluate(theory=None, masses=(100.,200.,210.,300.,400.), intensities=(100.,50.,50.,25.,10.), **kwargs):
    return evaluate_candidate_evidence("SM","SM(d18:1/16:0)",synthetic_sm() if theory is None else theory,
        pd.DataFrame({"fragment_mz":masses,"fragment_intensity":intensities}),
        PrecursorEvidence(True,0.,"[M+H]+"),total_composition="SM(d34:1)",**kwargs)


def test_structural_quality_is_max_while_gate_keeps_unmatched_theory():
    result = evaluate()
    assert result.summary["secondary_theoretical_count"] == 3
    assert result.summary["S_LCB"] == pytest.approx(1.05*.5/(.5+.05))
    assert result.summary["S_HG"] == pytest.approx(1.)
    expected = 60 + 20*(1.05*.5/.55) + 20
    assert result.summary["multi_evidence_score"] == pytest.approx(expected)
    assert result.summary["annotation_level"] == "molecular_species"


def test_full_spectrum_normalization_is_scale_invariant():
    a = evaluate(masses=(100.,200.,210.,300.,400.,999.),intensities=(100.,50.,50.,25.,10.,200.))
    b = evaluate(masses=(100.,200.,210.,300.,400.,999.),intensities=(1000.,500.,500.,250.,100.,2000.))
    assert a.summary["multi_evidence_score"] == pytest.approx(b.summary["multi_evidence_score"])
    assert a.summary["S_HG"] == pytest.approx(1.1*.5/.6)


def test_k_behavior():
    a = evaluate(masses=(100.,999.),intensities=(50.,100.),config=EvidenceScoringConfig(primary_half_intensity=.10))
    b = evaluate(masses=(100.,999.),intensities=(50.,100.),config=EvidenceScoringConfig(primary_half_intensity=.05))
    assert a.summary["S_HG"] < b.summary["S_HG"]


def test_explanation_fragments_do_not_dilute_denominator():
    theory = synthetic_sm()
    extras = annotate_theoretical_fragments(pd.DataFrame({"idf":["unreviewed"+str(i) for i in range(20)],"mz":[500.+i for i in range(20)]}),"SM")
    assert evaluate(pd.concat([theory,extras])).summary["multi_evidence_score"] == evaluate(theory).summary["multi_evidence_score"]


def test_duplicate_mass_priority_scores_once_and_no_centroid_reuse():
    theory = synthetic_sm()
    # Force HG and a LCB isobar; only HG may consume this experimental evidence.
    theory.loc[1,"theoretical_mz"] = 100.
    result = evaluate(theory,masses=(100.,210.,300.,400.),intensities=(100.,50.,25.,10.))
    assert result.fragments.theoretical_mz.is_unique
    assert result.fragments[result.fragments.theoretical_mz.eq(100.)].fragment_type.iloc[0] == "HG"
    assert result.summary["primary_theoretical_count"] == 1 and result.summary["secondary_theoretical_count"] == 2
    nearby = synthetic_sm()
    nearby.loc[1,"theoretical_mz"] = 100.0001
    result = evaluate(nearby,masses=(100.,),intensities=(100.,))
    assert result.fragments.matched.sum() == 1
    assert result.summary["S_LCB"] == 0.


def test_diagnostic_nl_unmapped_for_sm_does_not_affect_score():
    theory = synthetic_sm()
    result = evaluate(theory[theory.fragment_type.ne("NL")])
    assert result.summary["multi_evidence_score"] == evaluate().summary["multi_evidence_score"]
    assert result.summary["S_NL"] is None


def test_precursor_fragment_cannot_enter_msms_score():
    from sphingolipid_toolkit.evidence import FragmentOrigin
    precursor = annotate_theoretical_fragments(pd.DataFrame({"idf":["[M+H]+"],"mz":[500.],"fragment_origin":[FragmentOrigin("MS1")]}),"SM")
    a = evaluate()
    b = evaluate(pd.concat([synthetic_sm(),precursor]))
    assert a.summary["multi_evidence_score"] == b.summary["multi_evidence_score"]
    assert "Precursor" not in set(b.fragments.fragment_type)


def test_rt_cannot_change_level_gates_or_score():
    result = evaluate()
    rt = add_rt_validation(result,"fail",error=1.)
    assert all(rt.summary[key] == value for key,value in result.summary.items())


@pytest.mark.parametrize("values", [(float("nan"),),(-1.,),(float("inf"),)])
def test_invalid_intensity_rejected(values):
    with pytest.raises(ValueError):
        evaluate(masses=(100.,),intensities=values)


def test_empty_or_zero_spectrum_has_zero_scores_and_fails_gates():
    for masses,intensities in [((),()),((100.,),(0.,))]:
        result = evaluate(masses=masses,intensities=intensities)
        assert result.summary["multi_evidence_score"] == 0
        assert result.summary["confidence_score"] is None
        assert result.summary["annotation_level"] == "candidate"
