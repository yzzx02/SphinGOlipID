"""Local measured features, unchanged MS2 decisions and complete XLSX output."""
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from sphingolipid_toolkit import untargeted as workflow
from sphingolipid_toolkit.engine import SphingolipidIdentifier


SOURCE=Path(__file__).parent/'fixtures'/'untargeted_two_local_peaks.mzML'


def config(tmp_path,**values):
    return workflow.UntargetedConfig((SOURCE,),tmp_path/'result.xlsx',**values)


@pytest.mark.parametrize('values',[
    {'peakwidth_min':0},{'peakwidth_min':61},{'peakwidth_max':float('nan')},
    {'noise':-1},{'noise':float('inf')},{'snthresh':0},
    {'polarity':'auto'},{'top_n':0},{'top_n':1.5},
    {'precursor_ppm':0},{'fragment_ppm':None},{'min_total_score':101},
    {'fragment_da':.3},
])
def test_invalid_untargeted_settings_fail_before_running_r(tmp_path,values):
    with pytest.raises(ValueError):config(tmp_path,**values)


def test_da_tolerance_and_no_target_list_are_supported(tmp_path):
    settings=config(tmp_path,fragment_ppm=None,fragment_da=.3)
    assert settings.identification_config().fragment_da==.3
    assert settings.top_n==1 and settings.noise==1000
    assert not hasattr(settings,'target_file')
    with pytest.raises(ValueError):
        workflow.UntargetedConfig((SOURCE,SOURCE),tmp_path/'result.xlsx')


def native_features(path=SOURCE):
    return workflow.normalize_xcms_features(pd.DataFrame({
        'xcms_peak_id':['CP1','CP2'],'mz':[538.5225925995]*2,
        'mzmin':[538.5225925995]*2,'mzmax':[538.5225925995]*2,
        'rt':[40.,110.],'rtmin':[32.,102.],'rtmax':[48.,118.],
        'into':[971137.326171875,700000.],'maxo':[100100.,80100.],
        'sn':[558.,500.]}),path,'+')


def test_native_peak_ids_areas_and_seconds_are_preserved():
    frame=native_features()
    assert frame.xcms_peak_id.tolist()==['CP1','CP2']
    assert frame.peak_area.tolist()==[971137.326171875,700000.]
    assert frame.rt_min.tolist()==pytest.approx([32/60,102/60])
    assert frame.xcms_rt_min.tolist()==pytest.approx([40/60,110/60])


def test_apices_come_from_actual_local_ms1_not_nominal_xcms_time():
    frame=native_features()
    frame.loc[:,'xcms_rt_min']=[.68,1.84]
    measured=workflow.measure_local_apices(frame,SOURCE,'+')
    assert measured.apex_rt_min.tolist()==pytest.approx([40/60,110/60])
    assert measured.apex_intensity.tolist()==pytest.approx([100100.,80100.])
    assert measured.apex_native_scan_id.str.startswith('scan=').all()
    assert measured.peak_area.tolist()==frame.peak_area.tolist()


def spectrum(path=SOURCE,rt=40.2/60):
    return SimpleNamespace(precursor_mz=538.5225925995,rt_minutes=rt,polarity='+',
        metadata={'source':str(path),'precursor_refinement':dict(status='confirmed',
            ms1_mz=538.5225925995,ms1_native_id='scan=41',ms1_rt_seconds=rt*60-.2)})


def measured_features(path=SOURCE):
    frame=native_features(path)
    frame['apex_rt_min']=[40/60,110/60]
    return frame


def test_same_mass_neighboring_peaks_remain_separate():
    frame=measured_features()
    assert workflow.match_local_feature(frame,spectrum(),10)==0
    assert workflow.match_local_feature(frame,spectrum(rt=110.2/60),10)==1
    assert workflow.match_local_feature(frame,spectrum(rt=80/60),10) is None


@pytest.mark.parametrize('case',['missing_ms1','wrong_file','wrong_mode','outside_mass','overlapping_peaks'])
def test_ms2_cannot_create_or_borrow_a_local_feature(case,tmp_path):
    frame=measured_features()
    query=spectrum()
    if case=='missing_ms1':query.metadata['precursor_refinement']['status']='no_ms1'
    if case=='wrong_file':query.metadata['source']=str(tmp_path/'another.mzML')
    if case=='wrong_mode':query.polarity='-'
    if case=='outside_mass':query.precursor_mz+=.1
    if case=='overlapping_peaks':frame.loc[1,['rt_min','rt_max']]=frame.loc[0,['rt_min','rt_max']].to_numpy()
    assert workflow.match_local_feature(frame,query,10) is None


def test_one_best_acquisition_keeps_its_candidate_order_ranks_and_scores():
    matches=pd.DataFrame({
        'source_path':['one']*4+['two'],'polarity':['+']*5,'xcms_peak_id':['CP1']*5,
        'native_scan_id':['scan=1','scan=1','scan=2','scan=2','scan=1'],
        'total_score':[80.,79.,85.,84.,99.], 'matched_fragment_count':[3,5,4,2,6],
        'rank':[1,2,1,2,1],'ms2_evidence_json':['A','B','C','D','E']})
    original=matches.copy(deep=True)
    selected=workflow.select_feature_annotations(matches)
    assert selected.ms2_evidence_json.tolist()==['C','D','E']
    assert selected['rank'].tolist()==[1,2,1]
    assert selected.total_score.tolist()==[85.,84.,99.]
    pd.testing.assert_frame_equal(matches,original)


def test_equal_scores_choose_more_fragments_without_merging_spectra():
    matches=pd.DataFrame({'source_path':['one']*2,'polarity':['+']*2,'xcms_peak_id':['CP1']*2,
        'native_scan_id':['scan=1','scan=2'],'total_score':[80.,80.],
        'matched_fragment_count':[2,3],'rank':[1,1],'ms2_evidence_json':['first','second']})
    selected=workflow.select_feature_annotations(matches)
    assert selected.ms2_evidence_json.tolist()==['second']


@pytest.fixture
def detected_features(monkeypatch):
    monkeypatch.setattr(workflow,'find_rscript',lambda path=None:Path('synthetic-Rscript'))
    def detect(path,settings,work,logger):
        return native_features(path),{'xcms':'test fixture','algorithm':'centWave'}
    monkeypatch.setattr(workflow,'run_xcms',detect)


def test_feature_pipeline_preserves_exact_reference_results_and_per_spectrum_evidence(tmp_path,detected_features):
    settings=config(tmp_path)
    path=workflow.run_untargeted(settings)
    rows=pd.read_excel(path,sheet_name='spectrum_matches')
    annotations=pd.read_excel(path,sheet_name='annotations')
    features=pd.read_excel(path,sheet_name='features')
    with SphingolipidIdentifier(settings.identification_config()) as identifier:
        direct=identifier.identify_mzml(SOURCE)
    assert rows.native_scan_id.tolist()==direct.native_scan_id.tolist()
    assert rows.matched_name.tolist()==direct.matched_name.tolist()
    assert rows.total_score.tolist()==direct.total_score.tolist()
    assert rows.ms2_evidence_json.tolist()==direct.ms2_evidence_json.tolist()
    assert rows.result_rank.tolist()==direct.result_rank.tolist()
    assert len(annotations)==1 and annotations.matched_name.iloc[0]=='Cer(d18:1/16:0)'
    assert features.annotation_status.tolist()==['identified','no_passing_ms2_candidate']
    assert features.qualified_ms2_count.tolist()==[2,1]
    assert features.identified_ms2_count.tolist()==[2,0]
    assert rows.peak_area.tolist()==[971137.326171875]*2


def test_missing_local_detection_never_creates_a_new_peak(tmp_path,monkeypatch,detected_features):
    def detect(path,settings,work,logger):
        return native_features(path).iloc[:1].copy(),{'xcms':'test fixture'}
    monkeypatch.setattr(workflow,'run_xcms',detect)
    path=workflow.run_untargeted(config(tmp_path))
    features=pd.read_excel(path,sheet_name='features')
    audit=pd.read_excel(path,sheet_name='spectrum_audit')
    assert features.xcms_peak_id.tolist()==['CP1']
    assert audit.association_status.tolist()==['linked','linked','no_unique_local_peak']
    assert pd.isna(audit.xcms_peak_id.iloc[-1])


def test_empty_detection_exports_a_valid_workbook_without_ms2_only_fallback(tmp_path,monkeypatch,detected_features):
    monkeypatch.setattr(workflow,'run_xcms',lambda path,*args:(native_features(path).iloc[:0],{'xcms':'test fixture'}))
    path=workflow.run_untargeted(config(tmp_path))
    assert pd.read_excel(path).empty
    assert pd.read_excel(path,sheet_name='features').empty


def test_xcms_failure_does_not_replace_existing_output(tmp_path,monkeypatch,detected_features):
    settings=config(tmp_path)
    settings.output_file.write_bytes(b'previous workbook')
    def fail(*args):raise RuntimeError('actual XCMS failure')
    monkeypatch.setattr(workflow,'run_xcms',fail)
    with pytest.raises(RuntimeError,match='XCMS failure'):workflow.run_untargeted(settings)
    assert settings.output_file.read_bytes()==b'previous workbook'
    assert not list(tmp_path.glob('.sphingolipid_xcms_*'))
