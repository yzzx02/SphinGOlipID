"""Synthetic observed spectra exercise the real legacy query and finalizer."""
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.ms2_pipeline import _finalize_results
from sphingolipid_toolkit.scoring import add_standard_score_columns, match_fragments, score_fragment_matches
from sphingolipid_toolkit.targeted_mzml_pipeline import (
    MS2Spectrum, match_observed_to_theoretical, finalize_spectrum_matches,
)


def theory(masses, name="A"):
    return pd.DataFrame({"idf":[f"f{i}" for i in range(len(masses))], "mz":masses,
        "anno":name,"total":"synthetic", "target":500.,"lenDB":len(set(masses)),
        "RT":5.,"Abund":1000.})


def run_matchers(tmp_path, monkeypatch, masses, intensities, theoretical, ppm=20, top_n=3):
    observed = pd.DataFrame({"idx":"scan=1","A":masses,"B":intensities})
    theory_sheet = theoretical.assign(
        下限=theoretical.mz*(1-ppm/1e6),上限=theoretical.mz*(1+ppm/1e6))
    monkeypatch.setattr(core,"folder",str(tmp_path)+os.sep)
    monkeypatch.setattr(core,"MIN_FRAGMENT_INTENSITY",20.)
    # These in-memory readers stand in for Excel, but execute process_sheet/query.
    core.process_sheet("synthetic", SimpleNamespace(parse=lambda _:observed.copy()),
                       SimpleNamespace(parse=lambda _:theory_sheet.copy()))
    legacy = _finalize_results(tmp_path/"out-query.csv",tmp_path/"result.xlsx",1,0.,top_n)
    spectrum = MS2Spectrum(Path("synthetic.mzML"),"scan=1",500.,5.,25.,np.array(masses),np.array(intensities))
    raw = match_observed_to_theoretical(spectrum,theoretical,ppm/1e6,20.)
    targeted = add_standard_score_columns(finalize_spectrum_matches(raw,1,0.,top_n))
    api = match_fragments(
        observed.rename(columns={"A":"fragment_mz","B":"fragment_intensity"}),
        theoretical.rename(columns={"mz":"theoretical_mz"}),ppm,20.)
    return legacy,targeted,api


def test_one_observed_peak_two_nearby_theories_is_legacy_one_to_many(tmp_path, monkeypatch):
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[100.],[50.],theory([100.0001,100.0002]))
    for result in [legacy,targeted]:
        assert result.matched_fragment_count.tolist() == [2]
        assert result["实际mz"].iloc[0] == [100.,100.]
        assert result["强度总和"].iloc[0] == 100.
        assert result["匹配度分数"].iloc[0] == 1.
    assert score_fragment_matches(api)["matched_fragment_count"] == 2
    assert score_fragment_matches(api)["matched_intensity_sum"] == 100.


def test_multiple_observed_centroids_keep_highest_intensity_per_theory(tmp_path, monkeypatch):
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[100.,100.0005],[50.,80.],theory([100.0001]))
    for result in [legacy,targeted]:
        assert result.matched_fragment_count.tolist() == [1]
        assert result["强度总和"].iloc[0] == 80.
        assert result["实际mz"].iloc[0] == [100.0005]
    assert score_fragment_matches(api)["matched_fragment_count"] == 1
    assert score_fragment_matches(api)["matched_intensity_sum"] == 80.


def test_identical_theory_masses_are_counted_once(tmp_path, monkeypatch):
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[100.],[50.],theory([100.,100.]))
    assert legacy.matched_fragment_count.tolist() == targeted.matched_fragment_count.tolist() == [1]
    assert len(api) == 1


def test_legacy_candidate_scoring_and_top_n_are_locked(tmp_path, monkeypatch):
    theoretical = pd.concat([theory([100.,200.,300.],"A"),
        theory([100.,200.,400.,500.],"B"), theory([100.],"C")],ignore_index=True)
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[100.,200.,300.],[50.,30.,20.],theoretical,top_n=2)
    for result in [legacy,targeted]:
        assert result["注释"].tolist() == ["A","C"]
        assert result["总分数"].tolist() == [100.,80.]
        assert result["匹配度分数"].tolist() == [1.,1.]
        assert result.matched_fragment_count.tolist() == [3,1]
    assert {name:score_fragment_matches(group)["matched_fragment_count"]
            for name,group in api.groupby("anno")} == {"A":3,"B":2,"C":1}


@pytest.mark.parametrize("ppm,expected", [(10,1),(20,2)])
def test_real_matchers_share_explicit_ppm_and_intensity_cutoff(tmp_path, monkeypatch, ppm, expected):
    theoretical = theory([100.,200.,300.])
    # 15 ppm offset distinguishes 10/20; the 300 peak fails intensity cutoff.
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[100.,200.003,300.],[50.,30.,19.],theoretical,ppm)
    assert legacy.matched_fragment_count.tolist() == targeted.matched_fragment_count.tolist() == [expected]
    assert len(api) == expected


@pytest.mark.parametrize("direction", [-1,1])
def test_inclusive_ppm_endpoints_use_production_arithmetic(tmp_path,monkeypatch,direction):
    theoretical = theory([200.1234])
    observed = 200.1234*(1+direction*20/1e6)
    legacy,targeted,api = run_matchers(tmp_path,monkeypatch,[observed],[50.],theoretical)
    assert legacy.matched_fragment_count.tolist() == targeted.matched_fragment_count.tolist() == [1]
    assert len(api) == 1
