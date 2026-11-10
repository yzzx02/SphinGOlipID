from dataclasses import replace
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sphingolipid_toolkit.backend_policy import RemovedBackendError
from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig,run_batch
from sphingolipid_toolkit.targeted_mzml_pipeline import TargetedMzMLConfig,MS2Spectrum,annotate_one_spectrum
from sphingolipid_toolkit.engine import IdentificationConfig,SphingolipidIdentifier
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import LibraryRecord,FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint


@pytest.fixture
def current_inputs(tmp_path):
    import sphingolipid_toolkit._vendor.lipidgate.ms2 as ms2
    record=LibraryRecord(17,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',538.5199,'[M+H]+',
        'C34H67NO3','+',[
            FragmentRecord(264.2686,'LCB-2H2O','LCB碎片'),
            FragmentRecord(282.2791,'LCB-H2O','LCB碎片'),
            FragmentRecord(520.5093,'M+H-H2O','Common')])
    bank=tmp_path/'current.sqlite'
    write_library_index(bank,[record],dict(version=LIBRARY_CACHE_VERSION,source_sha256='fixture',
        rules_sha256=code_fingerprint(Path(ms2.__file__).parent)))
    raw=tmp_path/'raw';raw.mkdir()
    (raw/'hilic-msms-1.txt').write_text('\n'.join([
        'spectrum:','index:1','scan start time,1.25,min','target m/z,538.5199, m/z',
        'binaryDataArray:','binary: 0 264.2686','binaryDataArray:','binary: 0 100.0']),encoding='utf-8')
    pd.DataFrame(dict(target=[538.5199],下限=[538.51],上限=[538.53],RT=[1.25],Abund=[1000.])).to_excel(raw/'hilic-msms-1.xlsx',index=False)
    return raw,bank


def make_config(kind,tmp_path,**changes):
    if kind=='compat':return MS2PipelineConfig(tmp_path,tmp_path/'out',**changes)
    if kind=='unified':return SphinGOlipIDConfig(tmp_path,tmp_path,None,tmp_path/'out',**changes)
    return TargetedMzMLConfig(tmp_path,None,tmp_path/'out',**changes)


@pytest.mark.parametrize('kind',['compat','unified','targeted'])
@pytest.mark.parametrize('backend',['legacy','unknown'])
def test_every_configuration_rejects_noncurrent_backends(tmp_path,kind,backend):
    with pytest.raises(RemovedBackendError):make_config(kind,tmp_path,backend=backend)


@pytest.mark.parametrize('kind',['compat','unified','targeted'])
@pytest.mark.parametrize('name,value',[
    ('min_fragment_intensity',20.),('min_matched_fragments',2),('min_match_score',.35)])
def test_retired_scoring_options_cannot_silently_affect_new_identification(tmp_path,kind,name,value):
    with pytest.raises(RemovedBackendError,match=name):make_config(kind,tmp_path,**{name:value})


@pytest.mark.parametrize('kind',['compat','unified'])
def test_both_text_entrypoints_use_current_one_lcb_gate_and_scores(current_inputs,tmp_path,kind):
    raw,bank=current_inputs
    config=(MS2PipelineConfig(raw,tmp_path/'out',file_indices=(1,),ms2_library_path=bank)
        if kind=='compat' else SphinGOlipIDConfig(raw,raw,None,tmp_path/'out',file_indices=(1,),ms2_library_path=bank))
    assert config.backend=='lipidgate'
    result=pd.read_excel(run_batch(config)[0])
    assert result.matched_name.tolist()==['Cer(d18:1/16:0)']
    assert result.total_score.tolist()==[60.]
    assert result.matched_fragment_count.tolist()==[1]
    assert result.lcb_score.tolist()==[60.]
    assert result.other_score.tolist()==[0.]
    assert result.annotation_engine.tolist()==['SphinGOlipID']


def test_targeted_single_spectrum_matches_same_current_engine(current_inputs,tmp_path):
    raw,bank=current_inputs
    config=TargetedMzMLConfig(tmp_path,None,tmp_path/'out',ms2_library_path=bank)
    spectrum=MS2Spectrum(raw/'sample.mzML','scan=11',538.5199,1.25,60.,np.array([264.2686]),np.array([100.]))
    feature=pd.Series(dict(feature_id='target-row-1',feature_mz=538.5199,ms1_rt=1.25))
    actual=annotate_one_spectrum(spectrum,feature,pd.DataFrame(),config)
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank,fragment_ppm=config.fragment_ppm)) as identifier:
        expected=identifier.identify(538.5199,[(264.2686,100.)],scan_id='scan=11',rt_minutes=1.25)
    assert actual.matched_name.tolist()==[expected[0]['matched_name']]
    assert actual.total_score.tolist()==[expected[0]['total_score']]
    assert actual.ms2_evidence_json.tolist()==[expected[0]['ms2_evidence_json']]


def test_compat_config_is_revalidated_before_running_mutated_backend(current_inputs,tmp_path):
    raw,bank=current_inputs
    config=MS2PipelineConfig(raw,tmp_path/'out',file_indices=(1,),ms2_library_path=bank)
    config.backend='legacy'
    with pytest.raises(RemovedBackendError):run_batch(config)


@pytest.mark.parametrize('module,name',[
    ('ms2_legacy_core','learn_fuc'),('ms2_pipeline','_finalize_results'),
    ('targeted_mzml_pipeline','generate_legacy_fragments_for_candidates'),
    ('targeted_mzml_pipeline','finalize_spectrum_matches'),
    ('multi_evidence_adapter','evaluate_legacy_candidate'),
    ('multi_evidence_adapter','generate_evidence_fragments_for_candidates'),
    ('multi_evidence_scoring','evaluate_candidate_evidence'),
    ('multi_evidence_scoring','evaluate_annotation_gates')])
def test_old_and_experimental_identification_routes_fail_explicitly(module,name):
    with pytest.raises(RemovedBackendError):
        getattr(importlib.import_module('sphingolipid_toolkit.'+module),name)()


@pytest.mark.parametrize('module',['ms2_match_original','ms2_version3_2plus'])
def test_converted_notebook_entrypoints_cannot_reactivate_old_backend(module):
    with pytest.raises(RemovedBackendError):
        importlib.import_module('sphingolipid_toolkit.legacy_converted.'+module)


@pytest.mark.parametrize('flag',['--backend','--min-intensity','--min-fragments','--min-score','--ms1-db'])
def test_text_cli_rejects_old_backend_flags(tmp_path,flag):
    from sphingolipid_toolkit.cli import build_parser
    with pytest.raises(SystemExit) as error:
        build_parser().parse_args(['--input-dir',str(tmp_path),'--output-dir',str(tmp_path/'out'),flag,'legacy'])
    assert error.value.code==2


def test_aggregate_preserves_external_evidence_paths_in_subfolders(tmp_path):
    import logging
    from sphingolipid_toolkit.ms2_pipeline import _write_batch_result
    from sphingolipid_toolkit.result_writer import write_identifications
    source=tmp_path/'fraction1';source.mkdir()
    (source/'existing.json.gz').write_bytes(b'opaque evidence')
    row=pd.DataFrame(dict(ms2_evidence_json=['{"external_evidence_sha256":"abc"}'],
        evidence_file=['existing.json.gz'],evidence_sha256=['abc'],total_score=[60.],result_rank=[1],rank=[1]))
    part=source/'result.xlsx';write_identifications(row,part)
    out=tmp_path/'combined.xlsx'
    _write_batch_result([part],out,logging.getLogger('test'))
    result=pd.read_excel(out)
    assert (out.parent/result.evidence_file.iloc[0]).resolve()==(source/'existing.json.gz').resolve()
    assert result.total_score.tolist()==[60.] and result.result_rank.tolist()==[1]


@pytest.mark.parametrize('route',['text','spectrum','directory'])
@pytest.mark.parametrize('retired',[
    {'backend':'legacy'}, {'backend':'lipidgate','min_match_score':.35}])
def test_direct_reference_adapters_reject_retired_controls(route,retired):
    from types import SimpleNamespace
    from sphingolipid_toolkit import reference_pipeline
    config=SimpleNamespace(**retired)
    with pytest.raises(RemovedBackendError):
        if route=='text':reference_pipeline.run_text_file(config,1)
        elif route=='spectrum':reference_pipeline.annotate_targeted_spectrum(None,None,config,None)
        else:reference_pipeline.annotate_targeted_directory(None,config,None)
