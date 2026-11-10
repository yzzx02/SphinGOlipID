from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import CandidateScore,ExperimentalSpectrum,FragmentRecord,LibraryRecord,normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import CHAIN_CONFIRMATION_CLASSES,has_chain_evidence
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher


def pair(extra_h=False):
    p=828.6923
    common=[FragmentRecord(630.6183,'M+H-C6H10O5-2H2O','Diagnostic_HG'),
            FragmentRecord(648.6289,'M+H-C6H10O5-H2O','Diagnostic_HG'),
            FragmentRecord(666.6395,'M+H-C6H10O5','Common'),
            FragmentRecord(810.6817,'M+H-H2O','Common'),
            FragmentRecord(p,'[M+H]+','Precursor Ion')]
    h=LibraryRecord(0,'HexCer','HexCer(t42:1)','HexCer(d18:1/h24:0)',p,'[M+H]+',fragments=[
        FragmentRecord(264.2686,'LCB-2H2O','LCB碎片'),
        FragmentRecord(282.2791,'LCB-H2O','LCB碎片'),*common])
    if extra_h:h.fragments.insert(0,FragmentRecord(252.2686,'LCB-CH2O-H2O','Common'))
    t=replace(h,record_id=1,lipid_chain_name='HexCer(t18:0/24:1)',fragments=[
        FragmentRecord(264.2686,'LCB-3H2O','Common'),
        FragmentRecord(282.2791,'LCB-2H2O','LCB碎片'),
        FragmentRecord(300.2897,'LCB-H2O','LCB碎片'),
        FragmentRecord(318.3003,'LCB','Common'),*common])
    return h,t


def spectrum(h,kind):
    peaks=[(630.6183,1000.),(648.6289,1000.),(h.precursor_mz,30.)]
    if kind=='t_higher':peaks += [(264.2686,200.),(282.2791,100.),(318.3003,100.)]
    elif kind=='h_higher':peaks += [(252.2686,200.),(282.2791,200.)]
    elif kind=='weak_t':peaks += [(264.2686,200.),(282.2791,5.)]
    else:peaks=[(f.mz,1000.) for f in h.fragments]
    return ExperimentalSpectrum(kind,h.precursor_mz,1.,'+',normalize_peaks(peaks))


@pytest.mark.parametrize('kind',['t_higher','h_higher','equal','equal_extra_fragment'])
@pytest.mark.parametrize('reverse',[False,True])
def test_indexed_positive_isomers_reach_normal_topn_without_t_predeletion(tmp_path,kind,reverse):
    h,t=pair(kind in {'h_higher','equal_extra_fragment'})
    spec=spectrum(h,kind)
    path=tmp_path/'isomers.sqlite'
    write_library_index(path,[t,h] if reverse else [h,t],dict(version=LIBRARY_CACHE_VERSION,
        rules_sha256=code_fingerprint(Path(__file__).resolve().parents[2]/'src/sphingolipid_toolkit/_vendor/lipidgate/ms2')))
    fast=LipidMS2Searcher(path,native_min_candidates=0,fragment_tolerance_da=.001)
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False,fragment_tolerance_da=.001)
    try:
        scores=[fast._score_sphingo_candidate(spec,r) for r in (h,t)]
        snapshot=deepcopy(scores)
        assert all(r.passed_required_gates and r.resolution_level=='chain_level' for r in scores)
        assert all(has_chain_evidence(r.record,r.matched_fragments) for r in scores)
        assert not any(m.fragment.name=='LCB-H2O' for m in scores[1].matched_fragments)
        assert fast._resolve_trihydroxy_lcb_isomers(scores)==scores==snapshot
        if kind=='h_higher':assert scores[0].total_score>scores[1].total_score
        elif kind=='t_higher':assert scores[1].total_score>scores[0].total_score
        else:assert scores[0].total_score==scores[1].total_score==100.
        actual=fast.score_spectrum(spec,top_n=2)
        assert actual==slow.score_spectrum(spec,top_n=2)
        with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
            assert actual==slow.score_spectrum(spec,top_n=2)
        assert {r['matched_name'] for r in actual}=={h.lipid_chain_name,t.lipid_chain_name}
        if kind=='equal':assert [r['result_rank'] for r in actual]==[1,1]
        elif kind=='equal_extra_fragment':
            assert [(r['matched_name'],r['result_rank']) for r in actual]==[(h.lipid_chain_name,1),(t.lipid_chain_name,2)]
        else:
            expected=h if kind=='h_higher' else t
            assert actual[0]['matched_name']==expected.lipid_chain_name
            assert [r['result_rank'] for r in actual]==[1,2]
            assert len(fast.score_spectrum(spec,top_n=1))==1
    finally:fast.close();slow.close()


def test_weak_t_confirmation_is_not_promoted_or_given_secondary_exception():
    h,t=pair();s=object.__new__(LipidMS2Searcher);s.fragment_tolerance_da=.001
    spec=spectrum(h,'weak_t')
    scores=[s._score_sphingo_candidate(spec,r) for r in (h,t)]
    assert scores[0].resolution_level=='chain_level' and scores[1].resolution_level=='species_level'
    snapshot=deepcopy(scores)
    assert s._resolve_trihydroxy_lcb_isomers(scores)==scores==snapshot
    assert s._confirmed_isomer_candidate_ids(scores)==set()
    collapsed=s._collapse_report_equivalent_results(scores)
    assert collapsed==[scores[0]]


@pytest.mark.parametrize('cls',sorted(CHAIN_CONFIRMATION_CLASSES))
def test_all_shared_positive_classes_keep_both_confirmed_isomers_despite_unequal_scores(cls):
    h,t=pair()
    h=replace(h,compound_class=cls,lipid_name=f'{cls}(t42:1)',lipid_chain_name=f'{cls}(d18:1/h24:0)')
    t=replace(t,compound_class=cls,lipid_name=f'{cls}(t42:1)',lipid_chain_name=f'{cls}(t18:0/24:1)')
    from sphingolipid_toolkit._vendor.lipidgate.ms2.models import FragmentMatch,ExperimentalPeak
    a=CandidateScore(h,90.,True,[],0.,resolution_level='chain_level',matched_fragments=[
        FragmentMatch(h.fragments[1],ExperimentalPeak(282.2791,100.,.1),0.)])
    b=CandidateScore(t,80.,True,[],0.,resolution_level='chain_level',matched_fragments=[
        FragmentMatch(t.fragments[1],ExperimentalPeak(282.2791,100.,.1),0.)])
    assert LipidMS2Searcher._resolve_trihydroxy_lcb_isomers([a,b])==[a,b]
    assert LipidMS2Searcher._confirmed_isomer_candidate_ids([a,b])=={id(a),id(b)}
    b.passed_required_gates=False
    assert LipidMS2Searcher._confirmed_isomer_candidate_ids([a,b])==set()


@pytest.mark.parametrize('marker',[False,True])
def test_negative_isomer_selection_is_unchanged(marker):
    h,t=pair();s=object.__new__(LipidMS2Searcher);s.fragment_tolerance_da=.001
    spec=spectrum(h,'t_higher')
    if marker:spec=ExperimentalSpectrum('negative-control',h.precursor_mz,1.,'+',normalize_peaks([
        *[(p.mz,p.intensity) for p in spec.peaks],(300.2897,1000.)]))
    scores=[s._score_sphingo_candidate(spec,r) for r in (h,t)]
    for result in scores:result.record=replace(result.record,adduct='[M-H]-')
    assert s._resolve_trihydroxy_lcb_isomers(scores)==[scores[1] if marker else scores[0]]
    assert s._confirmed_isomer_candidate_ids(scores)==set()
