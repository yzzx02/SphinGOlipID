from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from sphingolipid_toolkit.config import TEXT_FRAGMENT_PPM,TARGETED_FRAGMENT_PPM
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig,run_one_file
from sphingolipid_toolkit.targeted_mzml_pipeline import TargetedMzMLConfig,MS2Spectrum,annotate_one_spectrum
from test_current_backend import current_inputs

@pytest.mark.parametrize('ppm,expected',[(10.,1),(20.,2)])
def test_current_text_tolerance_reaches_fragment_matching(current_inputs,tmp_path,ppm,expected):
    raw,bank=current_inputs
    shifted=282.2791*(1+15e-6)
    path=raw/'hilic-msms-1.txt'
    path.write_text(path.read_text().replace('binary: 0 264.2686','binary: 0 264.2686 '+str(shifted)).replace('binary: 0 100.0','binary: 0 100.0 20.0'))
    config=MS2PipelineConfig(raw,tmp_path/'out',file_indices=(1,),ms2_library_path=bank,fragment_ppm=ppm)
    result=pd.read_excel(run_one_file(config,1))
    assert result.matched_fragment_count.tolist()==[expected]
    assert result.ms2_tolerance_ppm.tolist()==[ppm]

@pytest.mark.parametrize('ppm,expected',[(10.,1),(20.,2)])
def test_current_targeted_tolerance_reaches_fragment_matching(current_inputs,tmp_path,ppm,expected):
    raw,bank=current_inputs
    config=TargetedMzMLConfig(tmp_path,None,tmp_path/'out',ms2_library_path=bank,fragment_ppm=ppm)
    spectrum=MS2Spectrum(raw/'sample.mzML','scan=1',538.5199,1.25,60.,np.array([264.2686,282.2791*(1+15e-6)]),np.array([100.,20.]))
    feature=pd.Series(dict(feature_id='target-1',feature_mz=538.5199,ms1_rt=1.25))
    result=annotate_one_spectrum(spectrum,feature,pd.DataFrame(),config)
    assert result.matched_fragment_count.tolist()==[expected]
    assert result.ms2_tolerance_ppm.tolist()==[ppm]

def test_current_tolerance_defaults_preserve_both_workflow_defaults(tmp_path):
    assert MS2PipelineConfig(tmp_path,tmp_path/'out').fragment_ppm==TEXT_FRAGMENT_PPM==20.
    assert TargetedMzMLConfig(tmp_path,None,tmp_path/'out').fragment_ppm==TARGETED_FRAGMENT_PPM==10.
