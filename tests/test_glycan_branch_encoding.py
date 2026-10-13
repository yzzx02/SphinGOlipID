"""Final-manuscript structures with artificial precursor masses only."""
import pandas as pd
import pytest

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.glycan_encoding import parse_glycan_encoding
from sphingolipid_toolkit.targeted_mzml_pipeline import (
    _read_masshunter_compound_csv, generate_legacy_fragments_for_candidates,
)


def fragments(sequence, positions=""):
    core.frag_id1, core.frag_mz1 = [], []
    core.GSL_fragments(sequence, positions, 1500., "synthetic")
    return {round(x, 6) for x in core.frag_mz1}


def test_linear_sequence():
    parsed = parse_glycan_encoding("Gal-Glc")
    assert parsed.main_chain == ("Gal", "Glc")
    assert not parsed.branches
    assert fragments("Gal-Glc") == {1337.9472, 1319.9366, 1175.8944, 1157.8838}


def test_final_manuscript_branch_first_and_retaining_losses():
    sequence = "Gal-Gal(-Fuc)-GlcNAc-Gal-Glc"
    parsed = parse_glycan_encoding(sequence)
    assert parsed.main_chain == ("Gal", "Gal", "GlcNAc", "Gal", "Glc")
    assert (parsed.branches[0].position, parsed.branches[0].residues) == (1, ("Fuc",))
    values = fragments(sequence)
    # Terminal Gal leaves with Fuc retained; Fuc can leave first; when its
    # attachment Gal leaves, Fuc must leave too and must not be subtracted twice.
    for loss in [162.0528, 146.0579, 308.1107, 470.1635, 673.2435]:
        assert round(1500-loss, 6) in values
        assert round(1500-loss-18.0106, 6) in values


def test_gm1_preserves_legacy_core_fragments_and_neuac_diagnostics():
    old = fragments(" -Gal -GalNAc -NeuAc -Gal -Glc", "0,3")
    new = fragments("Gal-GalNAc-Gal(-NeuAc)-Glc")
    # Empty legacy sentinel additionally emits unfragmented precursor/water.
    assert old - {1500., 1481.9894} <= new
    assert {292.1, 274.09} <= new
    parsed = parse_glycan_encoding(" -Gal -GalNAc -NeuAc -Gal -Glc", "0,3")
    normalized = parse_glycan_encoding("Gal-GalNAc-Gal(-NeuAc)-Glc")
    assert parsed.main_chain == normalized.main_chain
    assert parsed.branches == normalized.branches


@pytest.mark.parametrize("sequence", ["Gal#Fuc-Glc", "Gal(-Fuc", "Gal(Fuc)-Glc",
    "Gal(-Fuc(-Gal))-Glc", "Gal--Glc", "Gal(-Unknown)-Glc", "Gal-",
    "Gal(-Fuc-Gal)-Glc", "Gal(-Fuc)(-NeuAc)-Glc"])
def test_undefined_or_malformed_syntax_fails(sequence):
    with pytest.raises(ValueError):
        parse_glycan_encoding(sequence)


@pytest.mark.parametrize("sequence,positions", [
    ("Gal-GalNAc-Gal(-NeuAc)-Glc", ""),
    (" -Gal -GalNAc -NeuAc -Gal -Glc", "0,3")])
def test_csv_explicit_metadata_reaches_production(tmp_path, sequence, positions):
    path = tmp_path / "synthetic.csv"
    data = pd.DataFrame({"Formula":["C34H67NO3"],"RT":[5.],"Mass":[1498.99],
        "Cpd":["GM1(d34:1)"],"Comments":["synthetic"],
        "glycan_encoding":[sequence],"branch_positions":[positions]})
    path.write_text("# " + data.to_csv(index=False), encoding="utf-8")
    candidates = _read_masshunter_compound_csv(path)
    assert candidates.structure.iloc[0] == sequence
    assert candidates.classy.iloc[0] == positions
    result = generate_legacy_fragments_for_candidates(candidates,1500.,5.,1.,tmp_path)
    selected = result[result.anno.eq("GM1(d18:1/16:0)")]
    assert not selected.empty
    assert not selected.mz.duplicated().any()
    assert any(abs(selected.mz-1208.91) < 1e-6)
