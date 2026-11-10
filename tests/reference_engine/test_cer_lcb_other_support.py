from dataclasses import replace
import json
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.evidence_export import evidence_json
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalPeak,ExperimentalSpectrum,FragmentMatch,FragmentRecord,LibraryRecord,normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.rules import DEFAULT_RULES
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _calculate_pool_scores,_pool_for_fragment
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from test_indexed_library import identity


def records():
    other=[FragmentRecord(mz,name,'Common') for mz,name in (
        (250.2529,'LCB-CH2O-H2O'),(544.5452,'M+H-CH2O-H2O'),
        (556.5452,'M+H-2H2O'),(574.5558,'M+H-H2O'))]
    other.append(FragmentRecord(592.5663,'[M+H]+','Precursor Ion'))
    a=LibraryRecord(0,'Cer','Cer(d38:2)','Cer(d18:2/20:0)',592.5663,'[M+H]+',fragments=[
        FragmentRecord(262.2529,'LCB-2H2O','LCB碎片'),FragmentRecord(280.2635,'LCB-H2O','LCB碎片'),*other])
    b=replace(a,record_id=1,lipid_chain_name='Cer(d19:0/19:2)',fragments=[
        FragmentRecord(280.2999,'LCB-2H2O','LCB碎片'),FragmentRecord(298.3104,'LCB-H2O','LCB碎片'),
        replace(other[0],mz=268.2999),*other[1:]])
    return a,b


def searcher():
    s=object.__new__(LipidMS2Searcher);s.fragment_tolerance_da=.3
    return s


def test_cer_example_two_lcb_vs_one_becomes_equal_score_top1_without_changing_gates_or_matches(tmp_path):
    a,b=records();path=tmp_path/'cer-support.sqlite'
    write_library_index(path,[a,b],identity())
    spectrum=ExperimentalSpectrum('case',592.5663,1.,'+',normalize_peaks([
        (262.3401,165.2),(280.3414,190.9),(268.3037,32.3),(592.5663,20.),(9000.,1000.)]))
    fast=LipidMS2Searcher(path,native_min_candidates=0,fragment_tolerance_da=.3)
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False,fragment_tolerance_da=.3)
    try:
        rows=fast.score_spectrum(spectrum,top_n=1)
        assert rows==slow.score_spectrum(spectrum,top_n=1)
        assert {r['matched_name'] for r in rows}=={a.lipid_chain_name,b.lipid_chain_name}
        assert all(r['result_rank']==1 and r['total_score']==69.9785 and r['matched_fragment_count']==3 for r in rows)
        for r in rows:assert r['lcb_score']==pytest.approx(43.31179099346854)
        evidence=json.loads(next(r['ms2_evidence_json'] for r in rows if r['matched_name']==a.lipid_chain_name))
        aux=next(f for f in evidence['fragments'] if f['name']=='LCB-2H2O')
        assert aux['auxiliary_lcb_support'] and aux['scoring_pool']=='other' and aux['library_pool']=='lcb'
        assert aux['fragment_type']=='LCB碎片' and aux['contributes_to_score']
    finally:fast.close();slow.close()


@pytest.mark.parametrize('strength,aux',[(4.999,0),(5.,0),(5.001,1),(165.2,1)])
@pytest.mark.parametrize('overrides',[None,{0:1.,1:.2}])
def test_other_aux_uses_raw_threshold_and_the_effective_strongest_lcb(strength,aux,overrides):
    a,_=records();s=searcher()
    spectrum=ExperimentalSpectrum('boundary',a.precursor_mz,1.,'+',normalize_peaks([
        (262.3401,strength),(280.3414,190.9),(592.5663,20.),(9000.,1000.)]))
    matches=s._match_fragments_for_record(spectrum,a)
    adjustment=None if overrides is None else {id(m):overrides[i] for i,m in enumerate(matches[:2])}
    kwargs=dict(spectrum=spectrum,quality_relative_intensity_overrides=adjustment)
    x=_calculate_pool_scores(matches,a,DEFAULT_RULES.get('Cer'),**kwargs)
    with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
        y=_calculate_pool_scores(matches,a,DEFAULT_RULES.get('Cer'),**kwargs)
    assert x==y
    assert x['lcb'].matched_count==2 and x['lcb'].total_count==2
    if overrides is None:
        assert len(x['other'].auxiliary_fragment_indexes)==aux
        assert x['lcb'].scoring_fragment_indexes==(1,)
        assert x['other'].matched_count==1+aux
    else:
        assert x['lcb'].scoring_fragment_indexes==(0,)
        assert x['other'].auxiliary_fragment_indexes==(1,)


def test_duplicate_centroid_aliases_cannot_credit_both_pools_or_repeat_other():
    a,_=records()
    p=ExperimentalPeak(280.2635,100.,1.);q=ExperimentalPeak(262.2529,20.,.2)
    matches=[FragmentMatch(a.fragments[1],p,0.),FragmentMatch(a.fragments[0],q,0.),
             FragmentMatch(a.fragments[0],q,0.),FragmentMatch(a.fragments[-1],p,0.)]
    x=_calculate_pool_scores(matches,a,DEFAULT_RULES.get('Cer'))
    assert x['lcb'].scoring_fragment_indexes==(1,)
    assert x['other'].matched_count==1 and x['other'].auxiliary_fragment_indexes==(0,)
    credited=x['lcb'].scoring_fragment_indexes+x['other'].scoring_fragment_indexes
    assert len(credited)==len(set(credited))==2
    assert x['other'].pool_score==pytest.approx(40/3)


@pytest.mark.parametrize('adduct,name',[('[M-H]-','Cer(d18:2/20:0)'),('[M+H]+','Cerd18:2/20:0(O-16:0)')])
def test_negative_and_esterified_cer_do_not_receive_auxiliary_routing(adduct,name):
    a,_=records();a=replace(a,adduct=adduct,lipid_chain_name=name)
    matches=[FragmentMatch(f,ExperimentalPeak(f.mz,100.,1.),0.) for f in a.fragments]
    x=_calculate_pool_scores(matches,a,DEFAULT_RULES.get('Cer'))
    assert all(not p.auxiliary_fragment_indexes and not p.scoring_fragment_indexes for p in x.values())


def test_aux_does_not_create_chain_confirmation_or_raise_other_past_40():
    a,_=records();s=searcher();before=evidence_program(a)
    spectrum=ExperimentalSpectrum('weak',a.precursor_mz,1.,'+',normalize_peaks([
        (262.2529,6.),(280.2635,8.),*[(f.mz,1000.) for f in a.fragments[2:]]]))
    x=s._score_sphingo_candidate(spectrum,a)
    assert x.passed_required_gates and x.resolution_level=='species_level'
    assert x.pool_scores['other'].pool_score==40.
    assert evidence_program(a)==before
    assert all(_pool_for_fragment(a,f)=='lcb' for f in a.fragments[:2])


@pytest.mark.parametrize('cls',['SM','Hex2Cer','LacCer','Hex3Cer'])
def test_other_lipid_classes_keep_their_existing_scoring_routes(cls):
    a,_=records();a=replace(a,compound_class=cls,lipid_name=f'{cls}(d38:2)',lipid_chain_name=f'{cls}(d18:2/20:0)')
    matches=[FragmentMatch(f,ExperimentalPeak(f.mz,100.,1.),0.) for f in a.fragments]
    x=_calculate_pool_scores(matches,a,DEFAULT_RULES.get(cls))
    assert all(not p.auxiliary_fragment_indexes and not p.scoring_fragment_indexes for p in x.values())
