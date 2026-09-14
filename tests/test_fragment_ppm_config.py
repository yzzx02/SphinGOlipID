import io
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sphingolipid_toolkit import cli
from sphingolipid_toolkit import targeted_mzml_pipeline as targeted
from sphingolipid_toolkit.config import LEGACY_FRAGMENT_PPM, TARGETED_FRAGMENT_PPM, SphinGOlipIDConfig
from sphingolipid_toolkit.ms2_pipeline import MS2PipelineConfig, run_one_file


def synthetic_library():
    return pd.DataFrame({"理论值":[500.],"name":["So(d18:1)"],"classy":["So"],"structure":[""]})


def capture_logger():
    stream = io.StringIO()
    logger = logging.Logger("synthetic-ppm",level=logging.DEBUG)
    logger.addHandler(logging.StreamHandler(stream))
    return logger,stream


@pytest.mark.parametrize("ppm,expected", [(10,1),(20,2)])
def test_config_reaches_legacy_bounds_and_final_matches(tmp_path, ppm, expected):
    raw = tmp_path / "raw"
    raw.mkdir()
    offset_peak = 481.9894*(1+15/1e6)
    (raw/"hilic-msms-1.txt").write_text("\n".join([
        "spectrum:","index:1","scan start time,5.0,min","target m/z,500.0, m/z",
        "binaryDataArray:",f"binary: 0 60.044 {offset_peak} ",
        "binaryDataArray:","binary: 0 50.0 30.0 ",
    ]),encoding="utf-8")
    pd.DataFrame({"target":[500.],"下限":[499.9],"上限":[500.1],"RT":[5.],"Abund":[1000.]}).to_excel(raw/"hilic-msms-1.xlsx",index=False)
    library = tmp_path/"synthetic_library.xlsx"
    synthetic_library().to_excel(library,index=False)
    config = SphinGOlipIDConfig.from_single_input_dir(raw,tmp_path/"out",library,
        fragment_ppm=ppm,min_matched_fragments=1,min_match_score=0.,save_intermediate=True,file_indices=(1,))
    logger, log = capture_logger()
    result = pd.read_excel(run_one_file(config,1,logger=logger))
    assert result.matched_fragment_count.tolist() == [expected]
    bounds = pd.read_excel(config.output_dir/"intermediate"/"DB_sheet-1.xlsx")
    assert bounds["上限"].to_numpy() == pytest.approx(bounds.mz*(1+ppm/1e6))
    assert f"Effective fragment_ppm={ppm}" in log.getvalue()


@pytest.mark.parametrize("ppm,expected", [(10,1),(20,2)])
def test_config_reaches_targeted_generator_matcher_and_summary(tmp_path, ppm, expected):
    config = targeted.TargetedMzMLConfig(tmp_path,tmp_path/"unused.xlsx",tmp_path/"out",
        fragment_ppm=ppm,min_matched_fragments=1,min_match_score=0.)
    spectrum = targeted.MS2Spectrum(Path("synthetic.mzML"),"scan=1",500.,5.,25.,
        np.array([60.044,481.9894*(1+15/1e6)]),np.array([50.,30.]))
    feature = pd.Series({"feature_id":"synthetic","feature_mz":500.,"ms1_rt":5.})
    logger,log = capture_logger()
    result = targeted.annotate_one_spectrum(spectrum,feature,synthetic_library(),config,tmp_path/"temp",logger=logger)
    assert result.matched_fragment_count.tolist() == [expected]
    assert f"Effective fragment_ppm={ppm}" in log.getvalue()


@pytest.mark.parametrize("ppm", [None,10.,20.])
def test_cli_defaults_and_overrides_have_one_config_source(monkeypatch,tmp_path,ppm):
    configs = []
    monkeypatch.setattr(cli,"run_batch",lambda config,**kwargs: configs.append(config) or [])
    monkeypatch.setattr(targeted,"run_targeted_mzml_batch",lambda config,**kwargs: configs.append(config) or {})
    extra = [] if ppm is None else ["--fragment-ppm",str(ppm)]
    cli.main(["--input-dir",str(tmp_path),"--ms1-db","unused.xlsx","--output-dir",str(tmp_path)]+extra)
    targeted.main(["--data-dir",str(tmp_path),"--ms1-db","unused.xlsx","--output-dir",str(tmp_path)]+extra)
    assert configs[0].fragment_ppm == (LEGACY_FRAGMENT_PPM if ppm is None else ppm)
    assert configs[1].fragment_ppm == (TARGETED_FRAGMENT_PPM if ppm is None else ppm)
    assert MS2PipelineConfig(tmp_path,tmp_path,tmp_path/"unused.xlsx").fragment_ppm == LEGACY_FRAGMENT_PPM == 20.
    assert TARGETED_FRAGMENT_PPM == 10.
