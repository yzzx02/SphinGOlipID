"""Historical encodings, synthetic precursor masses; never loads a source library.

Templates: MS1 DB_new  3.0.xlsx (2025-07), Sheet1 structure/classy columns.
SHA256: 7719a35322067ceac1bca20447d1db2b3c092374c363c4975cab53c670e9c88e.
The leading space is significant: legacy split(' ') retains an empty token.
"""
import pandas as pd
import pytest

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.targeted_mzml_pipeline import (
    generate_legacy_fragments_for_candidates, _derive_legacy_classy, _derive_legacy_structure,
)


@pytest.fixture(autouse=True)
def isolate_legacy_fragment_globals(monkeypatch):
    monkeypatch.setattr(core,"frag_id1",[])
    monkeypatch.setattr(core,"frag_mz1",[])


@pytest.mark.parametrize("encoding,losses", [
    ("-Glc", [162.0528]),
    (" -Gal -Glc", [0.,162.0528,324.1056]),
])
def test_linear_glycan_losses_and_water_partners(encoding, losses):
    core.GSL_fragments(encoding,"",1000.,"synthetic")
    expected = [mz for loss in losses for mz in (1000-loss,1000-loss-18.0106)]
    assert core.frag_mz1 == pytest.approx(expected)


def test_gm3_sequential_losses_and_neuac_diagnostics():
    core.GSL_fragments(" -NeuAc -Gal -Glc","",1000.,"GM3")
    expected = [1000.,981.9894,708.91,690.8994,292.10,274.09,
                546.8572,528.8466,384.8044,366.7938]
    assert core.frag_mz1 == pytest.approx(expected)
    assert core.frag_id1[4:6] == ["NeuAc+H","NeuAc+H-H2O"]


def test_historical_gm1_branch_encoding_changes_generated_losses():
    encoding = " -Gal -GalNAc -NeuAc -Gal -Glc"
    core.GSL_fragments(encoding,"0,3",1500.,"GM1")
    losses = [0.,162.0528,365.1328,656.2228,818.2756,980.3284]
    expected = []
    for i, loss in enumerate(losses):
        expected.extend([1500-loss,1500-loss-18.0106])
        if i == 3:
            expected.extend([292.10,274.09])
    # Legacy location 0 uses [: -1] and re-emits sequential prefixes. Lock
    # those raw duplicates too; production deduplicates theoretical m/z later.
    for loss in losses[:-1]:
        expected.extend([1500-loss,1500-loss-18.0106])
    # Token 3 is NeuAc; legacy first loses the branch, then terminal Gal.
    expected.extend([1208.91,1190.8994,1046.8572,1028.8466])
    assert core.frag_mz1 == pytest.approx(expected)
    assert core.frag_id1[-4:] == ["GM1-NeuAc","GM1-NeuAc-H2O", "GM1-NeuAc-Gal","GM1-NeuAc-Gal-H2O"]


def test_branch_metadata_reaches_production_generator(tmp_path):
    # Actual template, but artificial mass/RT and a small chain candidate.
    candidates = pd.DataFrame({"name":["GM1(d34:1)"],
        "structure":[" -Gal -GalNAc -NeuAc -Gal -Glc"], "classy":["0,3"]})
    fragments = generate_legacy_fragments_for_candidates(candidates,1500.,5.,1.,tmp_path)
    molecular = fragments[fragments.anno.eq("GM1(d18:1/16:0)")]
    assert not molecular.empty
    assert molecular[molecular.idf.eq("GM1(d34:1)-NeuAc")].mz.iloc[0] == pytest.approx(1208.91)
    assert molecular[molecular.idf.eq("GM1(d34:1)-NeuAc-Gal")].mz.iloc[0] == pytest.approx(1046.8572)


def test_lcb_d181_diagnostics_reach_production_generator(tmp_path):
    candidates = pd.DataFrame({"name":["Cer(d34:1)"], "structure":[""], "classy":["Cer"]})
    fragments = generate_legacy_fragments_for_candidates(candidates,600.,5.,1.,tmp_path)
    molecular = fragments[fragments.anno.eq("Cer(d18:1/16:0)")]
    assert core.LCB_mz["d18:1"] == [264.2686,282.2792,300.2898]
    assert molecular[molecular.idf.isin(core.LCB_fragment["d18:1"])].mz.tolist() == [264.2686,282.2792,300.2898]


def test_hash_branch_syntax_is_not_supported_by_legacy():
    with pytest.raises((KeyError,ValueError)):
        core.GSL_fragments("-Gal #-NeuAc -Glc","",1500.,"synthetic")


def test_csv_name_fallback_lacks_branch_metadata_gap_is_explicit():
    assert _derive_legacy_structure("GM1(d34:1)") == "-Gal -GalNAc -NeuAc -Gal -Glc"
    assert _derive_legacy_classy("GM1(d34:1)") == ""
