from dataclasses import replace
import json
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.evidence_export import evidence_json
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.rules import DEFAULT_RULES
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _calculate_pool_scores, _pool_for_fragment
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from test_indexed_library import identity
from test_sphingolipid_isomer_ranking import pair


@pytest.mark.parametrize('lcb_strengths,top_names',[
    ((19.5726,11.5532),['HexCer(d18:1/h24:0)']),
    ((42.2477,10.1385),['HexCer(d18:1/h24:0)']),
    ((25.7131,33.3267),['HexCer(d18:1/h24:0)','HexCer(t18:0/24:1)']),
])
def test_hexcer_two_lcb_peaks_receive_independent_support_and_preserve_ranking(tmp_path,lcb_strengths,top_names):
    h,t=pair()
    path=tmp_path/'hexcer-support.sqlite'
    write_library_index(path,[h,t],identity())
    spectrum=ExperimentalSpectrum('two-lcb',h.precursor_mz,1.,'+',normalize_peaks([
        (264.2686,lcb_strengths[0]),(282.2791,lcb_strengths[1]),
        (630.6183,48.),(648.6289,65.),(810.6817,100.)]))
    fast=LipidMS2Searcher(path,native_min_candidates=0,fragment_tolerance_da=.001)
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False,fragment_tolerance_da=.001)
    try:
        rows=fast.score_spectrum(spectrum,top_n=2)
        assert rows==slow.score_spectrum(spectrum,top_n=2)
        with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
            assert rows==slow.score_spectrum(spectrum,top_n=2)
        assert sorted(r['matched_name'] for r in rows if r['result_rank']==1)==sorted(top_names)
        assert all(r['matched_fragment_count']==5 and r['resolution_level']=='chain_level' for r in rows)
        assert all(r['other_score']==pytest.approx(20*2/3) for r in rows)
        scores=[fast._score_sphingo_candidate(spectrum,r) for r in (h,t)]
        assert scores[0].pool_scores['hg']==scores[1].pool_scores['hg']
        assert scores[0].pool_scores['lcb'].matched_count==2
        assert scores[1].pool_scores['lcb'].matched_count==1
        audit=json.loads(evidence_json(spectrum,scores[0]))['fragments']
        auxiliary=[f for f in audit if f.get('auxiliary_lcb_support')]
        assert len(auxiliary)==1 and auxiliary[0]['library_pool']=='lcb'
        assert auxiliary[0]['scoring_pool']=='other' and auxiliary[0]['fragment_type']=='LCB碎片'
        credited=[f for f in audit if f.get('contributes_to_score')]
        assert len(credited)==len({f['observed_mz'] for f in credited})==3
        hg=[f for f in audit if f['library_pool']=='hg' and f['observed_mz'] is not None]
        assert len(hg)==2 and all(f.get('contributes_to_score',True) for f in hg)
    finally:
        fast.close();slow.close()


@pytest.mark.parametrize('strength,auxiliary',[(.4999,False),(.5,False),(.5001,True),(1.,True)])
def test_hexcer_auxiliary_threshold_keeps_hg_gate_and_shared_one_percent_confirmation(strength,auxiliary):
    h,_=pair();searcher=object.__new__(LipidMS2Searcher);searcher.fragment_tolerance_da=.001
    before=evidence_program(h)
    def score(hg,primary):
        peaks=[(264.2686,primary),(282.2791,strength),(810.6817,100.)]
        if hg:peaks.append((648.6289,100.))
        spectrum=ExperimentalSpectrum('boundary',h.precursor_mz,1.,'+',normalize_peaks(peaks))
        return searcher._score_sphingo_candidate(spectrum,h)
    result=score(True,10.)
    assert result.passed_required_gates and result.resolution_level=='chain_level'
    assert bool(result.pool_scores['other'].auxiliary_fragment_indexes)==auxiliary
    assert result.pool_scores['other'].pool_score==pytest.approx(20*(1+auxiliary)/3)
    assert result.pool_scores['hg'].pool_score==60.
    assert result.pool_scores['lcb'].pool_score<=20.
    assert not score(False,10.).passed_required_gates
    assert score(True,0.).resolution_level==('chain_level' if strength>=1. else 'species_level')
    assert evidence_program(h)==before
    assert all(_pool_for_fragment(h,f)=='lcb' for f in h.fragments[:2])


def test_hexcer_centroid_aliases_do_not_repeat_auxiliary_or_primary_credit():
    h,_=pair()
    primary=ExperimentalPeak(282.2791,100.,1.)
    secondary=ExperimentalPeak(264.2686,20.,.2)
    matches=[FragmentMatch(h.fragments[1],primary,0.),FragmentMatch(h.fragments[0],secondary,0.),
             FragmentMatch(h.fragments[0],secondary,0.),FragmentMatch(h.fragments[-1],primary,0.)]
    native=_calculate_pool_scores(matches,h,DEFAULT_RULES.get('HexCer'))
    with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
        assert native==_calculate_pool_scores(matches,h,DEFAULT_RULES.get('HexCer'))
    assert native['lcb'].scoring_fragment_indexes==(1,)
    assert native['other'].auxiliary_fragment_indexes==(0,)
    assert native['other'].matched_count==1
    assert native['other'].pool_score==pytest.approx(20/3)


@pytest.mark.parametrize('adduct',['[M+Na]+','[M-H]-'])
def test_hexcer_other_adducts_do_not_receive_auxiliary_support(adduct):
    h,_=pair();h=replace(h,adduct=adduct)
    matches=[FragmentMatch(f,ExperimentalPeak(f.mz,100.,1.),0.) for f in h.fragments]
    pools=_calculate_pool_scores(matches,h,DEFAULT_RULES.get('HexCer'))
    assert all(not p.auxiliary_fragment_indexes and not p.scoring_fragment_indexes for p in pools.values())
