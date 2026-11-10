from dataclasses import replace
import json
from pathlib import Path

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION, load_standard_msp
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import (
    C, H, N, O, PROTON, WATER, HEXOSE_RESIDUE, DIHEXOSYL_HG,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher


@pytest.mark.parametrize("compound_class", ["Hex2Cer", "LacCer"])
@pytest.mark.parametrize("series,carbon,db", [("d", 14, 0), ("d", 18, 1), ("d", 24, 2), ("t", 18, 0), ("t", 20, 1)])
def test_dihexosylceramide_lcb_pools_masses_gates_and_resolution(tmp_path, compound_class, series, carbon, db):
    oxygen = 3 if series == "t" else 2
    intact = carbon * C + (2 * carbon + 3 - 2 * db) * H + N + oxygen * O + PROTON
    precursor = (carbon + 16) * C + (2 * (carbon + 16) + 1 - 2 * db) * H + N + (oxygen + 1) * O + 2 * HEXOSE_RESIDUE + PROTON
    name = f"{compound_class}({series}{carbon}:{db}/16:0)"
    species = f"{compound_class}({series}{carbon + 16}:{db})"
    losses = {
        "M+H-C6H10O5": HEXOSE_RESIDUE, "M+H-C6H12O6": HEXOSE_RESIDUE + WATER,
        "M+H-C12H18O9": 2 * HEXOSE_RESIDUE - WATER, "M+H-C12H20O10": 2 * HEXOSE_RESIDUE,
        "M+H-C12H20O10-H2O": 2 * HEXOSE_RESIDUE + WATER,
        "M+H-C12H20O10-2H2O": 2 * HEXOSE_RESIDUE + 2 * WATER,
    }
    ions = [(precursor, "[M+H]+", "Precursor Ion"), (precursor - WATER, "M+H-H2O", "C类碎片")]
    ions += [(precursor - loss, ion, "Diagnostic_HG") for ion, loss in losses.items()]
    ions += [(100. + i, ion, "LCB碎片") for i, ion in enumerate(("LCB-H2O", "LCB-2H2O", "LCB-3H2O"))]
    source = tmp_path / "dihex.msp"
    lines = [f"Name: {name}", f"PrecursorMZ: {precursor:.4f}", "PrecursorType: [M+H]+",
             f"CompoundClass: {compound_class}", f"Comment: MS1_name={species};polarity=+", f"Num Peaks: {len(ions)}"]
    lines += [f'{mz:.4f} 100 "{ion}" "{role}"' for mz, ion, role in ions]
    source.write_text("\n".join(lines) + "\n", encoding="utf-8")
    record = next(r for r in load_standard_msp(source) if r.lipid_chain_name == name)
    fragments = {f.name: f for f in record.fragments}
    lcb_names = {"LCB-H2O", "LCB-2H2O"}
    assert {f.name for f in record.fragments if f.fragment_type == "LCB碎片"} == lcb_names
    assert {f.name for f in record.fragments if f.fragment_type == "Diagnostic_HG"} == DIHEXOSYL_HG
    assert fragments["M+H-H2O"].fragment_type == "Common"
    for ion, loss in (("LCB", 0), ("LCB-H2O", 1), ("LCB-2H2O", 2), ("LCB-3H2O", 3)):
        if ion in fragments:
            assert fragments[ion].mz == round(intact - loss * WATER, 4)
    if series == "d":
        hydroxy = next(r for r in load_standard_msp(source) if r.lipid_chain_name == name.replace("/16:0", "/h16:0"))
        assert next(f.mz for f in hydroxy.fragments if f.name == "M+H-H2O") == pytest.approx(precursor - WATER + O, abs=.0001)
    alias = FragmentRecord(fragments["LCB-H2O"].mz, "shared_lcb_support", "Common")
    record = replace(record, fragments=[*record.fragments, alias])
    bank = tmp_path / "dihex.sqlite"
    write_library_index(bank, [record], dict(version=LIBRARY_CACHE_VERSION,
                        rules_sha256=code_fingerprint(Path(__file__).resolve().parents[2] / "src/sphingolipid_toolkit/_vendor/lipidgate/ms2")))
    fast = LipidMS2Searcher(bank, native_min_candidates=0)
    reference = LipidMS2Searcher(bank, use_native_engine=False, use_score_prefilter=False)
    hg = [(fragments[ion].mz, 1000.) for ion in sorted(DIHEXOSYL_HG)[:3]]
    lcbs = [fragments[ion].mz for ion in sorted(lcb_names)]
    other_lcb = [fragments[ion].mz for ion in ("LCB", "LCB-3H2O")] if series == "t" else []
    moved = (fragments["M+H-H2O"].mz, 1000.)
    partial = [(fragments[n].mz,1000.) for n in ('M+H-C6H10O5','M+H-C6H12O6','M+H-C12H18O9','M+H-C12H20O10')]
    try:
        for peaks, expected in (
            (partial + [moved] + [(mz, 1000.) for mz in lcbs + other_lcb], None),
            ([(mz, 1.) for mz, _ in hg] + [(999., 1000.)], None),
            (hg, species),
            ([hg[0]], species),
            ([hg[1]], species),
            ([hg[0],(hg[1][0],10.)], species),
            ([(hg[0][0],10.),hg[1]], species),
            (hg + [(mz, 9.9) for mz in lcbs], species),
            (hg + [(mz, 1000.) for mz in other_lcb], species),
            (hg + [(lcbs[0], 10.)], name),
            (hg + [(lcbs[1], 10.)], name),
            (hg + [(fragments["LCB-H2O"].mz, 1000.)], name),
        ):
            spectrum = ExperimentalSpectrum("dihex", record.precursor_mz, 1., "+", normalize_peaks(peaks))
            rows = fast.score_spectrum(spectrum)
            assert rows == reference.score_spectrum(spectrum)
            assert [r["matched_name"] for r in rows] == ([expected] if expected else [])
            if rows:
                [row] = rows
                assert row["hg_score"] == 60.
                evidence = {f["name"]: f for f in json.loads(row["ms2_evidence_json"])["fragments"]}
                assert evidence["M+H-H2O"]["scoring_pool"] == "other"
                for ion in lcb_names:
                    assert evidence[ion]["scoring_pool"] == "lcb" and evidence[ion]["role"] == "chain"
                for ion in ("LCB", "LCB-3H2O") if series == "t" else ():
                    assert evidence[ion]["scoring_pool"] == "other" and evidence[ion]["role"] == "other"
                assert evidence["shared_lcb_support"]["observed_mz"] is None
                if peaks == hg + [(fragments["LCB-H2O"].mz, 1000.)]:
                    assert (row["lcb_score"], row["other_score"]) == (20., 0.)
    finally:
        fast.close()
        reference.close()
