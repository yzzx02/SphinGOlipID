"""Author-final 60/20/20 policy; artificial spectra and existing rule labels."""
import pandas as pd
import pytest

from sphingolipid_toolkit.evidence import annotate_theoretical_fragments
from sphingolipid_toolkit.evidence_config import EvidenceScoringConfig
from sphingolipid_toolkit.multi_evidence_scoring import (
    PrecursorEvidence, evaluate_candidate_evidence, structural_quality, support_count_score,
)


def run(lipid, labels, masses, observed_masses, intensities, *, precursor_mz=1000., overrides=None, config=EvidenceScoringConfig()):
    theory = annotate_theoretical_fragments(pd.DataFrame({"idf":labels,"mz":masses}),lipid,
        eligibility_overrides=overrides)
    return evaluate_candidate_evidence(lipid,lipid+"(d18:1/16:0)",theory,
        pd.DataFrame({"fragment_mz":observed_masses,"fragment_intensity":intensities}),
        PrecursorEvidence(True,0.,"[M+H]+",observed_mz=precursor_mz),config=config,
        total_composition=lipid+"(d34:1)")


@pytest.mark.parametrize("intensity",[0.,.1,.5,1.])
def test_quality_endpoints_bounds_and_lipidgate_formula(intensity):
    for k in (.10,.05):
        q = structural_quality(intensity,k)
        assert q == pytest.approx((1+k)*intensity/(intensity+k))
        assert 0 <= q <= 1
    if intensity in (0.,1.):
        assert structural_quality(intensity,.10) == intensity
    else:
        assert structural_quality(intensity,.10) < structural_quality(intensity,.05)


def test_structural_pool_uses_max_not_average_or_coverage():
    # Inverse q gives exactly q=.3 and q=.8 at k=.10. Unmatched base peak
    # intensity=1 sets normalization without contributing a theory match.
    intensities = [q*.10/(1.10-q) for q in (.3,.8)] + [1.]
    result = run("GM3",["NeuAc+H","NeuAc+H-H2O"],[292.10,274.09],[292.10,274.09,600.],intensities)
    assert result.summary["S_primary"] == pytest.approx(.8)
    assert result.summary["S_HG"] == pytest.approx(.8)
    assert result.summary["primary_matched_count"] == 2
    assert result.summary["primary_gate_pass"]


@pytest.mark.parametrize("lipid,labels,primary,secondary", [
    ("GM3",["NeuAc+H","d18:1","GM3(d34:1)-H2O"],"HG","LCB"),
    ("SM",["[phosphocholine+H]+","d18:1","SM(d34:1)-H2O"],"HG","LCB"),
    ("HexCer",["HexCer(d34:1)-Hex","d18:1","HexCer(d34:1)-H2O"],"diagnostic_nl","LCB"),
    ("LacCer",["LacCer(d34:1)-Gal-Glc","d18:1","LacCer(d34:1)-H2O"],"diagnostic_nl","LCB"),
    ("Cer",["d18:1","Cer(d34:1)-1H2O","[C2H5NO+H]+"],"LCB","diagnostic_nl"),
])
def test_fixed_score_mappings_and_aliases(lipid,labels,primary,secondary):
    result = run(lipid,labels,[100.,200.,300.],[100.,200.,300.],[100.,50.,1.],overrides={labels[0]:True})
    s = result.summary
    assert (s["primary_evidence_type"],s["secondary_evidence_type"]) == (primary,secondary)
    assert s["multi_evidence_score"] == pytest.approx(60+20*(1.05*.5/.55)+20)
    assert (s["primary_weight"],s["secondary_weight"],s["support_weight"]) == (60,20,20)
    for evidence,pool in ((primary,"primary"),(secondary,"secondary")):
        assert s["S_"+("NL" if evidence == "diagnostic_nl" else evidence)] == s["S_"+pool]
    assert s["S_Supp"] == s["S_supp"] == 1
    assert s["S_HG"] is None if primary != "HG" else s["S_NL"] is None


@pytest.mark.parametrize("total,matched,expected",[(8,1,1/3),(8,2,2/3),(8,3,1),(8,5,1),(2,1,.5),(2,2,1),(1,1,1)])
def test_support_count_saturation(total,matched,expected):
    labels = [f"SM(d34:1)-{i+1}H2O" for i in range(total)]
    masses = [100.+i for i in range(total)]
    result = run("SM",labels,masses,masses[:matched],[1.]*matched)
    assert result.summary["S_supp"] == pytest.approx(expected)
    assert result.summary["N_support_total"] == total
    assert result.summary["N_support_matched"] == matched
    assert result.fragments.q.isna().all()
    assert result.summary["annotation_level"] == "candidate"


def test_support_intensity_independence_and_supporting_nl_membership():
    labels = ["SM(d34:1)-H2O","SM(d34:1)-C3H9N","SM(d34:1)-C5H14NO4P"]
    a = run("SM",labels,[100.,200.,300.],[100.,200.,300.],[1.,1.,1.])
    b = run("SM",labels,[100.,200.,300.],[100.,200.,300.],[10000.,.01,1.])
    assert a.summary["S_supp"] == b.summary["S_supp"] == 1
    assert a.summary["N_support_total"] == 2  # diagnostic NL excluded
    roles = a.fragments.set_index("theoretical_mz")
    assert roles.loc[200.,"evidence_role"] == "supporting_nl"
    assert roles.loc[200.,"scoring_pool"] == "support"
    assert not roles.loc[200.,"gate_eligible"]
    assert roles.loc[300.,"evidence_role"] == "diagnostic_nl"
    assert pd.isna(roles.loc[300.,"scoring_pool"])


def test_missing_support_is_na_and_does_not_redistribute_twenty_points():
    result = run("SM",["[phosphocholine+H]+","d18:1"],[100.,200.],[100.,200.],[1.,1.])
    assert result.summary["support_status"] == "no_support_theory"
    assert result.summary["S_supp"] is None
    assert result.summary["S_Supp"] is None
    assert result.summary["multi_evidence_score"] == 80
    assert support_count_score(0,0) == (None,0)


def test_highest_quality_does_not_replace_half_gate():
    result = run("SM",["[phosphocholine+H]+"]*4,[100.,110.,120.,130.],[100.],[1.])
    assert result.summary["S_primary"] == 1
    assert result.summary["primary_theoretical_count"] == 4
    assert result.summary["primary_matched_count"] == 1
    assert result.summary["primary_gate_pass"] is False
    assert result.summary["annotation_level"] == "candidate"


@pytest.mark.parametrize("lipid,labels,overrides,selected",[
    ("Cer",["Cer(d34:1)-1H2O","d18:1","[C2H5NO+H]+"],None,"LCB"),
    ("HexCer",["HexCer(d34:1)-Hex","d18:1","[C2H5NO+H]+"],{"HexCer(d34:1)-Hex":True},"NL"),
])
def test_identical_mass_primary_priority_is_class_specific(lipid,labels,overrides,selected):
    for seed in range(4):
        table = pd.DataFrame({"idf":labels,"mz":[100.]*3}).sample(frac=1,random_state=seed)
        enriched = annotate_theoretical_fragments(table,lipid,eligibility_overrides=overrides)
        result = evaluate_candidate_evidence(lipid,lipid+"(d18:1/16:0)",enriched,
            pd.DataFrame({"fragment_mz":[100.],"fragment_intensity":[100.]}),PrecursorEvidence(True,0.,"[M+H]+"))
        assert len(result.fragments) == 1
        assert result.fragments.fragment_type.iloc[0] == selected
        assert result.fragments.scoring_pool.iloc[0] == "primary"
        assert result.summary["primary_matched_count"] == 1
        assert result.summary["secondary_matched_count"] == 0
        assert result.summary["N_support_matched"] == 0
        assert len(result.fragments.fragment_name.iloc[0].split(" | ")) == 3


def test_strong_precursor_does_not_depress_structural_quality():
    labels = ["[phosphocholine+H]+","d18:1","SM(d34:1)-H2O"]
    results = [run("SM",labels,[100.,200.,300.],[100.,200.,300.,1000.,1001.],[50.,100.,10.,p,p/2]) for p in (100.,1e10)]
    a,b = [r.summary for r in results]
    assert a["multi_evidence_score"] == b["multi_evidence_score"]
    assert b["non_precursor_base_peak_intensity"] == 100
    assert b["structural_normalization_status"] == "non_precursor_base_peak"


def test_window_is_inclusive_configurable_and_falls_back():
    # Entire spectrum lies inside exclusion window; non-precursor signal absent.
    result = run("SM",["[phosphocholine+H]+"],[1002.],[1000.,1002.],[100.,50.])
    assert result.summary["structural_normalization_status"] == "spectrum_relative_fallback"
    assert result.summary["S_primary"] == pytest.approx(structural_quality(.5,.1))
    tighter = run("SM",["[phosphocholine+H]+"],[1002.],[1000.,1002.],[100.,50.],
        config=EvidenceScoringConfig(precursor_cluster_exclusion_da=1.))
    assert tighter.summary["S_primary"] == 1
    for offset in (-4.1,4.1):
        # Use a precursor for which the floating-point difference is exactly
        # the same as the passed window, to lock inclusive comparison.
        mass = 1000.+offset
        width = abs(mass-1000.)
        boundary = run("SM",["[phosphocholine+H]+"],[mass],[1000.,mass],[100.,50.],
            config=EvidenceScoringConfig(precursor_cluster_exclusion_da=width))
        assert boundary.summary["non_precursor_base_peak_intensity"] == 0


def test_actual_precursor_centroid_cannot_be_structural_evidence():
    result = run("SM",["[phosphocholine+H]+"],[1000.],[1000.],[100.])
    assert result.summary["S_primary"] == 0
    assert result.summary["primary_matched_count"] == 0
    assert not result.summary["primary_gate_pass"]
