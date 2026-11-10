from dataclasses import replace
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import CandidateScore, FragmentMatch, ExperimentalPeak
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import (
    CHAIN_CONFIRMATION_CLASSES, HEXCER_LCB_POOLS, C, H, N, O, PROTON, WATER,
    has_chain_evidence, normalize_fragments,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.matching import match_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _pool_for_fragment
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import _expand_hydroxy_fa_sphingolipids, _without_d18_1_300_lcb_marker


@pytest.mark.parametrize('cls',sorted(CHAIN_CONFIRMATION_CLASSES))
@pytest.mark.parametrize('carbon,db',[(18,0),(20,1),(16,0)])
def test_t_removed_support_is_not_restored_in_any_pool(cls,carbon,db):
    name=f'{cls}(t{carbon}:{db}/16:0)'
    old=[FragmentRecord(252.2686,'LCB-CH6O3','Common')]
    fragments=normalize_fragments(cls,name,800.,old)
    assert 'LCB-CH6O3' not in {f.name for f in fragments}
    assert normalize_fragments(cls,name,800.,fragments)==fragments
    if cls in {'Cer','HexCer','Hex2Cer','LacCer'}:
        assert {f.name for f in fragments if f.fragment_type=='LCB碎片'}=={'LCB-H2O','LCB-2H2O'}


def record(cls,series,db=None):
    db=(1 if series=='d' else 0) if db is None else db
    intact=18*C+(39-2*db)*H+N+(3 if series=='t' else 2)*O+PROTON
    fragments=[FragmentRecord(round(intact-i*WATER,4),name,'Common')
               for i,name in enumerate(('LCB','LCB-H2O','LCB-2H2O','LCB-3H2O'))]
    fragments += [FragmentRecord(538.5194,'[M+H]+','Precursor Ion'),FragmentRecord(252.2686,'LCB-CH2O-H2O','Common')]
    return LibraryRecord(0,cls,f'{cls}({series}{34}:{db})',f'{cls}({series}18:{db}/16:0)',538.5194,'[M+H]+',fragments=fragments)


@pytest.mark.parametrize('cls',sorted(CHAIN_CONFIRMATION_CLASSES))
@pytest.mark.parametrize('series,db',[('d',1),('d',2),('t',0)])
def test_shared_confirmation_uses_two_series_specific_channels_independent_of_pool(cls,series,db):
    r=record(cls,series,db)
    for f in r.fragments:
        for relative in (.009999,.01,.010001):
            peaks=normalize_peaks([(f.mz,1000*relative),(1500.,1000.)])
            matches=match_fragments(peaks,r.fragments,lambda mz:.00001)
            assert has_chain_evidence(r,matches)==(f.name in HEXCER_LCB_POOLS[series] and relative>=.01)


@pytest.mark.parametrize('series',['d','t'])
def test_cer_half_lcb_gate_other_cannot_admit_and_raw_one_percent_only_changes_resolution(series):
    r=record('Cer',series)
    r=replace(r,fragments=normalize_fragments('Cer',r.lipid_chain_name,r.precursor_mz,r.fragments))
    expected=HEXCER_LCB_POOLS[series]
    assert {f.name for f in r.fragments if f.fragment_type=='LCB碎片'}==expected
    mask=sum(1<<i for i,f in enumerate(r.fragments) if f.name in expected)
    assert evidence_program(r)==(((mask,1),),)
    search=object.__new__(LipidMS2Searcher);search.fragment_tolerance_da=.00001
    for f in r.fragments:
        if f.name=='[M+H]+':continue
        assert _pool_for_fragment(r,f)==('lcb' if f.name in expected else 'other')
        for strength in (5.,10.):
            spec=ExperimentalSpectrum('cer',r.precursor_mz,1.,'+',normalize_peaks([(r.precursor_mz,1000.),(f.mz,strength)]))
            actual=search._score_sphingo_candidate(spec,r)
            with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
                reference=search._score_sphingo_candidate(spec,r)
            assert actual==reference
            assert actual.passed_required_gates==(f.name in expected)
            if f.name in expected:
                assert actual.resolution_level==('chain_level' if strength>=10. else 'species_level')
                assert actual.pool_scores['lcb'].pool_score==60.
            else:assert actual.resolution_level=='class_level'
    for f in r.fragments:
        if f.name in expected:
            spec=ExperimentalSpectrum('no-other',r.precursor_mz,1.,'+',normalize_peaks([(f.mz,1000.)]))
            score=search._score_sphingo_candidate(spec,r)
            assert score.passed_required_gates and score.total_score==60.


def test_cer_m_retains_existing_lcb_definitions_and_ceil_half_gate():
    fragments=[FragmentRecord(284.2948,'LCB','LCB碎片'),FragmentRecord(266.2842,'LCB-H2O','LCB碎片'),
               FragmentRecord(256.2635,'Ceramide fragment U','LCB碎片'),FragmentRecord(522.5245,'[M+H]+','Precursor Ion')]
    r=LibraryRecord(0,'Cer','Cer(m34:1)','Cer(m18:1/16:0)',522.5245,'[M+H]+',fragments=fragments)
    assert normalize_fragments('Cer',r.lipid_chain_name,r.precursor_mz,fragments)==fragments
    assert evidence_program(r)==(((7,2),),)
    search=object.__new__(LipidMS2Searcher);search.fragment_tolerance_da=.00001
    for n in (1,2):
        spec=ExperimentalSpectrum('m',r.precursor_mz,1.,'+',normalize_peaks([(f.mz,1000.) for f in fragments[:n]]))
        score=search._score_sphingo_candidate(spec,r)
        assert score.passed_required_gates==(n==2)
        if n==2:
            assert score.pool_scores['lcb'].pool_score==60. and score.pool_scores['other'].pool_score==pytest.approx(40./3)
            assert score.resolution_level=='chain_level'


def test_cer_hydroxy_expansion_shifts_precursor_losses_in_other_but_not_lcb_support():
    r=record('Cer','d')
    r=replace(r,fragments=[*r.fragments,FragmentRecord(520.5088,'M+H-H2O','Common')])
    r=replace(r,fragments=normalize_fragments('Cer',r.lipid_chain_name,r.precursor_mz,r.fragments))
    [original,derived]=_expand_hydroxy_fa_sphingolipids([r])
    masses={f.name:f.mz for f in original.fragments}
    for f in derived.fragments:
        delta=O if f.name in {'[M+H]+','M+H-H2O'} else 0.
        assert f.mz==pytest.approx(masses[f.name]+delta)


def test_sm_d18_1_common_annotation_does_not_restore_excluded_t_marker():
    r=record('SM','d')
    revised=_without_d18_1_300_lcb_marker(r)
    assert 'LCB' not in {f.name for f in revised.fragments}
    assert {'LCB-H2O','LCB-2H2O'} <= {f.name for f in revised.fragments}
    assert _without_d18_1_300_lcb_marker(replace(r,lipid_chain_name='SM(t18:0/16:0)'))==replace(r,lipid_chain_name='SM(t18:0/16:0)')


@pytest.mark.parametrize('cls',['Cer','SM','HexCer'])
@pytest.mark.parametrize('marker,mz',[('LCB-H2O',300.2897),('LCB-2H2O',282.2791)])
@pytest.mark.parametrize('relative',[.005,.01])
def test_hydroxy_vs_t_resolution_reuses_shared_confirmation_and_raw_threshold(cls,marker,mz,relative):
    h=LibraryRecord(0,cls,f'{cls}(t34:1)',f'{cls}(d18:1/h16:0)',554.5143,'[M+H]+',fragments=[])
    t=replace(h,record_id=1,lipid_chain_name=f'{cls}(t18:0/16:1)')
    fragment=FragmentRecord(mz,marker,'Common' if cls=='SM' else 'LCB碎片')
    match=FragmentMatch(fragment,ExperimentalPeak(mz,1000*relative,relative),0.)
    h_score=CandidateScore(h,80.,True,[],0.,resolution_level='chain_level')
    t_score=CandidateScore(t,80.,True,[],0.,resolution_level='chain_level' if relative>=.01 else 'species_level',matched_fragments=[match])
    assert LipidMS2Searcher._resolve_trihydroxy_lcb_isomers([h_score,t_score])==[h_score,t_score]
    assert t_score.resolution_level==('chain_level' if relative>=.01 else 'species_level')
