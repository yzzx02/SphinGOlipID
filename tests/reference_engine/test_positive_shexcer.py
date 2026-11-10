from dataclasses import replace
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import normalize_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _pool_for_fragment
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from test_indexed_library import identity


def raw_record(hydroxy):
    cls = 'SHexCer+O' if hydroxy else 'SHexCer'
    precursor = 824.5552 if hydroxy else 808.5603
    lcb = [FragmentRecord(mz, name, 'LCB碎片', intensity=80.+i) for i,(mz,name) in enumerate([
        (252.2686,'LCB-CH2O-H2O'),(264.2686,'LCB-2H2O'),(282.2791,'LCB-H2O')])]
    other = ([FragmentRecord(564.5350,'Cer+H-H2O','Common')] if hydroxy else [
        FragmentRecord(530.5295,'M+H-H2SO4-C6H10O5-H2O','Common'),
        FragmentRecord(548.5401,'M+H-H2SO4-C6H10O5','Diagnostic_HG')])
    return LibraryRecord(0,cls,f'{cls}({"t" if hydroxy else "d"}36:1)',
        f'{cls}(d18:1/{"h" if hydroxy else ""}18:0)',precursor,'[M+H]+',polarity='+',fragments=[
            *lcb,*other,FragmentRecord(726.5878 if hydroxy else 710.5929,'M+H-H2SO4','Diagnostic_HG'),
            FragmentRecord(precursor,'[M+H]+','Precursor Ion')])


def normalized_record(hydroxy):
    r = raw_record(hydroxy)
    return replace(r,fragments=normalize_fragments(r.compound_class,r.lipid_chain_name,r.precursor_mz,r.fragments))


@pytest.mark.parametrize('hydroxy',[False,True])
def test_shexcer_family_preserves_source_peaks_and_normalizes_only_roles(hydroxy):
    before=raw_record(hydroxy);after=normalized_record(hydroxy)
    assert [(f.name,f.mz,f.intensity) for f in before.fragments]==[(f.name,f.mz,f.intensity) for f in after.fragments]
    assert normalize_fragments(after.compound_class,after.lipid_chain_name,after.precursor_mz,after.fragments)==after.fragments
    assert [f.name for f in after.fragments if f.fragment_type=='LCB碎片']==['LCB-2H2O','LCB-H2O']
    assert len([f for f in after.fragments if f.fragment_type=='Diagnostic_HG'])==2
    if hydroxy:
        assert 'M+H-H2SO4-C6H10O5-H2O' not in {f.name for f in after.fragments}
    assert all(_pool_for_fragment(after,f)=='other' for f in after.fragments if f.name in {'LCB-CH2O-H2O','[M+H]+'})


@pytest.fixture(params=[False,True],ids=['SHexCer','SHexCer+O'])
def family_search(request,tmp_path):
    record=normalized_record(request.param)
    bank=tmp_path/'shexcer.sqlite'
    write_library_index(bank,[record],identity())
    fast=LipidMS2Searcher(bank,native_min_candidates=0,fragment_tolerance_da=.001,min_total_score=0.)
    slow=LipidMS2Searcher(bank,use_native_engine=False,use_score_prefilter=False,fragment_tolerance_da=.001,min_total_score=0.)
    yield record,fast,slow
    fast.close();slow.close()


def compare(record,fast,slow,peaks):
    spectrum=ExperimentalSpectrum('shexcer',record.precursor_mz,1.,'+',normalize_peaks(peaks))
    rows=fast.score_spectrum(spectrum,top_n=3)
    assert rows==slow.score_spectrum(spectrum,top_n=3)
    with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
        assert rows==slow.score_spectrum(spectrum,top_n=3)
    return spectrum,rows


def test_shexcer_family_either_hg_admits_without_lcb_or_other_and_no_fallback(family_search):
    record,fast,slow=family_search
    mask=sum(1<<i for i,f in enumerate(record.fragments) if f.fragment_type=='Diagnostic_HG')
    assert evidence_program(record)==(((mask,1),),)
    for f in record.fragments:
        _,rows=compare(record,fast,slow,[(f.mz,1000.)])
        if f.fragment_type=='Diagnostic_HG':
            assert len(rows)==1 and rows[0]['matched_name']==record.lipid_name
            assert rows[0]['resolution_level']=='species_level'
            assert (rows[0]['hg_score'],rows[0]['lcb_score'],rows[0]['other_score'])==(60.,0.,0.)
        else:
            assert rows==[]
    no_hg=[(f.mz,1000.) for f in record.fragments if f.fragment_type!='Diagnostic_HG']
    _,rows=compare(record,fast,slow,no_hg)
    assert rows==[]
    assert not fast._score_sphingo_candidate(ExperimentalSpectrum('no-hg',record.precursor_mz,1.,'+',normalize_peaks(no_hg)),record).passed_required_gates


@pytest.mark.parametrize('stronger_index',[0,1])
@pytest.mark.parametrize('stronger_intensity',[100.,1000.])
def test_shexcer_family_default_hg_score_uses_strongest_peak_without_weak_peak_dilution(
    family_search,stronger_index,stronger_intensity,
):
    record,fast,slow=family_search
    hg=[f for f in record.fragments if f.fragment_type=='Diagnostic_HG']
    strongest,weaker=hg[stronger_index],hg[1-stronger_index]
    single_peaks=[(strongest.mz,stronger_intensity),(9000.,1000.)]
    single,single_rows=compare(record,fast,slow,single_peaks)
    both,both_rows=compare(record,fast,slow,[*single_peaks,(weaker.mz,10.)])
    assert len(single_rows)==len(both_rows)==1
    assert both_rows[0]['hg_score']==single_rows[0]['hg_score']
    assert both_rows[0]['lcb_score']==both_rows[0]['other_score']==0.
    assert both_rows[0]['matched_fragment_count']==2
    for spectrum in (single,both):
        candidate=fast._score_sphingo_candidate(spectrum,record)
        assert candidate.passed_required_gates and candidate.resolution_level=='species_level'
        assert candidate.pool_scores['hg'].total_count==2
    assert fast._score_sphingo_candidate(both,record).pool_scores['hg'].matched_count==2
    if stronger_intensity==1000.:
        assert both_rows[0]['hg_score']==60.


@pytest.mark.parametrize('marker',['LCB-H2O','LCB-2H2O'])
@pytest.mark.parametrize('relative,level',[(.009999,'species_level'),(.01,'chain_level'),(.010001,'chain_level')])
def test_shexcer_family_shared_one_percent_confirmation_independent_of_pool(family_search,marker,relative,level):
    record,fast,slow=family_search
    f=next(f for f in record.fragments if f.name==marker)
    hg=next(f for f in record.fragments if f.name=='M+H-H2SO4')
    _,rows=compare(record,fast,slow,[(hg.mz,1000.),(f.mz,relative*1000.)])
    assert len(rows)==1 and rows[0]['resolution_level']==level
    assert rows[0]['matched_name']==(record.lipid_chain_name if level=='chain_level' else record.lipid_name)
    assert rows[0]['lcb_score']<=20. and rows[0]['other_score']<=20.
    assert rows[0]['lcb_score']>0. and rows[0]['other_score']==0.


@pytest.mark.parametrize('strength,other_score',[(4.999,0.),(5.,0.),(5.001,20/3)])
def test_shexcer_family_other_strict_half_percent_and_d252_cannot_confirm_chains(family_search,strength,other_score):
    record,fast,slow=family_search
    peaks={f.name:f.mz for f in record.fragments}
    _,rows=compare(record,fast,slow,[(peaks['M+H-H2SO4'],1000.),(peaks['LCB-CH2O-H2O'],strength)])
    assert len(rows)==1 and rows[0]['other_score']==pytest.approx(other_score)
    assert rows[0]['resolution_level']=='species_level'


def test_shexcer_family_raw_confirmation_is_not_raised_by_precursor_normalization(family_search):
    record,fast,slow=family_search
    f={f.name:f.mz for f in record.fragments}
    _,rows=compare(record,fast,slow,[(f['M+H-H2SO4'],100.),(f['LCB-2H2O'],5.),(record.precursor_mz,1000.)])
    assert len(rows)==1 and rows[0]['total_score']>50.
    assert rows[0]['resolution_level']=='species_level'
    assert (rows[0]['hg_score'],rows[0]['lcb_score'],rows[0]['other_score'])==pytest.approx((60.,10.5,20/3))


@pytest.mark.parametrize('primary',['LCB-H2O','LCB-2H2O'])
def test_shexcer_family_two_lcb_supports_are_scored_once_and_keep_library_identities(family_search,primary):
    record,fast,slow=family_search
    f={f.name:f.mz for f in record.fragments}
    secondary='LCB-2H2O' if primary=='LCB-H2O' else 'LCB-H2O'
    spectrum,rows=compare(record,fast,slow,[(f['M+H-H2SO4'],1000.),(f[primary],300.),
                                         (f[secondary],100.),(record.precursor_mz,30.)])
    assert len(rows)==1 and rows[0]['matched_fragment_count']==4
    assert rows[0]['lcb_score']==pytest.approx(18.) and rows[0]['other_score']==pytest.approx(40/3)
    candidate=fast._score_sphingo_candidate(spectrum,record)
    assert candidate.pool_scores['lcb'].matched_count==2 and candidate.pool_scores['lcb'].total_count==2
    evidence=json.loads(rows[0]['ms2_evidence_json'])['fragments']
    auxiliary=[p for p in evidence if p.get('auxiliary_lcb_support')]
    assert len(auxiliary)==1 and auxiliary[0]['name']==secondary
    assert auxiliary[0]['fragment_type']=='LCB碎片' and auxiliary[0]['library_pool']=='lcb'
    assert auxiliary[0]['scoring_pool']=='other'
    credited=[p for p in evidence if p.get('contributes_to_score')]
    assert len(credited)==len({p['observed_mz'] for p in credited})==3


@pytest.mark.parametrize('strength,other_score',[(4.999,0.),(5.,0.),(5.001,20/3)])
def test_shexcer_family_auxiliary_lcb_uses_strict_raw_half_percent(family_search,strength,other_score):
    record,fast,slow=family_search
    f={p.name:p.mz for p in record.fragments}
    _,rows=compare(record,fast,slow,[(f['M+H-H2SO4'],1000.),(f['LCB-2H2O'],100.),(f['LCB-H2O'],strength)])
    assert len(rows)==1 and rows[0]['other_score']==pytest.approx(other_score)
    assert rows[0]['resolution_level']=='chain_level'


def test_hydroxy_headgroup_alias_is_one_physical_peak_in_both_engines(tmp_path):
    record=normalized_record(True)
    alias=next(f for f in record.fragments if f.name=='Cer+H-H2O')
    record=replace(record,fragments=[*record.fragments,replace(alias,name='M+H-H2SO4-C6H10O5')])
    bank=tmp_path/'alias.sqlite';write_library_index(bank,[record],identity())
    fast=LipidMS2Searcher(bank,native_min_candidates=0)
    slow=LipidMS2Searcher(bank,use_native_engine=False,use_score_prefilter=False)
    try:
        _,rows=compare(record,fast,slow,[(alias.mz,1000.)])
        assert len(rows)==1 and rows[0]['hg_score']==60. and rows[0]['matched_fragment_count']==1
        fragments=json.loads(rows[0]['ms2_evidence_json'])['fragments']
        assert sum(f['observed_mz'] is not None for f in fragments)==1
    finally:fast.close();slow.close()


@pytest.mark.parametrize('adduct',['[M-H]-','[M+Na]+'])
def test_shexcer_curation_leaves_other_adducts_untouched(adduct):
    script=Path(__file__).resolve().parents[2]/'scripts/curate_positive_hexosylceramide_library.py'
    spec=importlib.util.spec_from_file_location('shexcer_curator',script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    block=['Name: SHexCer+O(d18:1/h18:0)','PrecursorMZ: 824.5552',f'PrecursorType: {adduct}',
           'CompoundClass: SHexCer+O','Num Peaks: 1','264.2686 80 "LCB-2H2O" "LCB碎片"']
    assert module.curate_block(block)==block
