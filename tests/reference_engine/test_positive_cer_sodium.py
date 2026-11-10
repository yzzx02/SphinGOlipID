from dataclasses import replace
import gzip
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.evidence_export import evidence_json
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION, load_standard_msp
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value, iter_blocks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_cer_sodium import CerSodiumComposition, GATE_NAME, OTHER_NAMES
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _pool_for_fragment, _pool_weights_for_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.rules import DEFAULT_RULES
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from scripts.extend_positive_cer_sodium_library import extend, seed_compositions, sodium_block


def record(oh=2):
    c = CerSodiumComposition(34, 1, oh)
    return LibraryRecord(0, "Cer", c.name, c.name, c.precursor_mz, "[M+Na]+",
                         formula=c.formula, polarity="+", fragments=c.fragments())


def spectrum(r, peaks):
    return ExperimentalSpectrum("CerNa", r.precursor_mz, 1., "+", normalize_peaks(peaks))


def score(r, peaks):
    search = object.__new__(LipidMS2Searcher)
    search.fragment_tolerance_da = .001
    s = spectrum(r, peaks)
    fast = search._score_sphingo_candidate(s, r)
    with patch("sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge", return_value=None):
        slow = search._score_sphingo_candidate(s, r)
    assert fast == slow
    return fast


@pytest.mark.parametrize("oh,p,formula", [(1,544.5064,"C34H67NO2"), (2,560.5013,"C34H67NO3"),
                                          (3,576.4962,"C34H67NO4"), (4,592.4911,"C34H67NO5")])
def test_cer_sodium_formula_losses_and_total_only_msp_roundtrip(tmp_path, oh, p, formula):
    c = CerSodiumComposition(34, 1, oh)
    assert c.precursor_mz == pytest.approx(p, abs=.00005)
    assert c.formula == formula
    f = {f.name:f for f in c.fragments()}
    assert set(f) == {GATE_NAME, *OTHER_NAMES}
    assert c.precursor_mz-f[GATE_NAME].mz == pytest.approx(18.01056468403)
    assert c.precursor_mz-f["M+Na-2H2O"].mz == pytest.approx(36.02112936806)
    assert c.precursor_mz-f["M+Na-CH2O-H2O"].mz == pytest.approx(48.02112936806)
    path = tmp_path/"cer_na.msp"
    path.write_text("\n".join(sodium_block(c))+"\n", encoding="utf-8")
    rows = load_standard_msp(path)
    assert len(rows)==1 and rows[0].lipid_name==rows[0].lipid_chain_name==c.name
    assert "/" not in rows[0].lipid_chain_name
    actual = score(rows[0], [(f[GATE_NAME].mz,1000.)])
    assert actual.passed_required_gates and actual.resolution_level=="species_level"
    assert actual.total_score==60.


@pytest.mark.parametrize("oh", [1,2,3,4])
def test_only_water_loss_admits_other_never_gates_or_confirms_chains(oh):
    r = record(oh)
    water = next(f for f in r.fragments if f.name==GATE_NAME)
    assert _pool_weights_for_record(r,DEFAULT_RULES.get("Cer"))=={"fah":0.,"hg":60.,"other":40.}
    assert all(_pool_for_fragment(r,f)==("hg" if f==water else "other") for f in r.fragments)
    gate_index = r.fragments.index(water)
    assert evidence_program(r)==(((1<<gate_index,1),),)
    only_gate = score(r,[(water.mz,1000.)])
    assert only_gate.passed_required_gates and only_gate.total_score==60.
    assert only_gate.pool_scores["other"].matched_count==0
    for other in [f for f in r.fragments if f!=water]:
        rejected = score(r,[(other.mz,1000.)])
        assert not rejected.passed_required_gates and rejected.resolution_level=="class_level"
    rejected = score(r,[(f.mz,1000.) for f in r.fragments if f!=water])
    assert not rejected.passed_required_gates
    full = score(r,[(f.mz,1000.) for f in r.fragments])
    assert full.passed_required_gates and full.total_score==100.
    assert full.resolution_level=="species_level" and full.downgrade_reason=="sum_composition_identity"
    audit = json.loads(evidence_json(spectrum(r,[(f.mz,1000.) for f in r.fragments]),full))
    assert sum(f["admission_gate"] for f in audit["fragments"])==1
    assert all(f["scoring_pool"]==("hg" if f["name"]==GATE_NAME else "other") for f in audit["fragments"])
    assert all(f["role"]!="chain" for f in audit["fragments"])


@pytest.mark.parametrize("strength,count", [(4.999,0),(5.,0),(5.001,1),(10.,1)])
def test_cer_sodium_other_threshold_uses_raw_intensity_and_has_no_gate(strength,count):
    r = record();f = {f.name:f for f in r.fragments}
    result = score(r,[(f[GATE_NAME].mz,1000.),(f["M+Na-2H2O"].mz,strength)])
    assert result.passed_required_gates and result.resolution_level=="species_level"
    assert result.pool_scores["other"].matched_count==count
    assert result.pool_scores["other"].pool_score==pytest.approx(count*40/3)


def test_extra_lcb_evidence_never_upgrades_sodium_resolution():
    r=record()
    r=replace(r,fragments=[*r.fragments,FragmentRecord(264.2686,"LCB-2H2O","LCB碎片")])
    result=score(r,[(f.mz,1000.) for f in r.fragments])
    assert result.passed_required_gates and result.resolution_level=="species_level"


def test_indexed_native_search_reports_one_sum_even_with_lcb_like_experimental_peaks(tmp_path):
    rows=[record(oh) for oh in [1,2,3,4]]
    for i,r in enumerate(rows):r.record_id=i
    path=tmp_path/"cer_na.sqlite"
    write_library_index(path,rows,dict(version=LIBRARY_CACHE_VERSION,
        rules_sha256=code_fingerprint(Path(__file__).resolve().parents[2]/"src/sphingolipid_toolkit/_vendor/lipidgate/ms2")))
    fast=LipidMS2Searcher(path,native_min_candidates=0)
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False)
    try:
        for r in rows:
            for include_other in [False,True]:
                peaks=[(f.mz,1000.) for f in r.fragments if include_other or f.name==GATE_NAME]
                spec=spectrum(r,[*peaks,(264.2686,500.),(282.2791,500.)])
                actual=fast.score_spectrum(spec)
                assert actual==slow.score_spectrum(spec)
                assert len(actual)==1 and actual[0]["matched_name"]==r.lipid_name
                assert actual[0]["resolution_level"]=="species_level"
                assert actual[0]["total_score"]==(100. if include_other else 60.)
            spec=spectrum(r,[(f.mz,1000.) for f in r.fragments if f.name!=GATE_NAME])
            assert fast.score_spectrum(spec)==slow.score_spectrum(spec)==[]
    finally:fast.close();slow.close()


def seed(name, cls="Cer", adduct="[M+H]+"):
    return [f"Name: {name}","PrecursorMZ: 600",f"PrecursorType: {adduct}",f"CompoundClass: {cls}",
            "Num Peaks: 1",'300 100 "existing" "Common"']


def test_curation_deduplicates_oxygen_isomers_keeps_low_lcb_excludes_short_acyl_and_esters(tmp_path):
    a=seed("Cer(d18:1/16:0)");b=seed("Cer(t18:0/16:1)")
    assert CerSodiumComposition(34,1,3) in seed_compositions(a)
    assert CerSodiumComposition(34,1,3) in seed_compositions(b)
    assert seed_compositions(seed("Cer(d8:1/14:0)"))
    assert not seed_compositions(seed("Cer(d18:1/13:0)"))
    assert not seed_compositions(seed("Cerd18:1/16:0(O-18:0)"))
    assert not seed_compositions(seed("Cer(d34:1)"))  # Never infer acyls from sums.
    assert not seed_compositions(seed("HexCer(d18:1/16:0)","HexCer"))
    blocks=[a,b,seed("Cer(d18:1/13:0)"),seed("Cerd18:1/16:0(O-18:0)")]
    source,out,repeat=[tmp_path/n for n in ["source.gz","sodium.gz","repeat.gz"]]
    with gzip.open(source,"wt",encoding="utf-8") as stream:
        stream.write("\n\n".join("\n".join(b) for b in blocks)+"\n")
    audit=extend(source,out)
    result=list(iter_blocks(out))
    assert result[:len(blocks)]==blocks
    assert audit["added_records"]==3
    added=result[len(blocks):]
    assert len({header_value(b,"Formula") for b in added})==3
    assert all(header_value(b,"PrecursorType")=="[M+Na]+" and "/" not in header_value(b,"Name") for b in added)
    assert extend(out,repeat)["added_records"]==0
    assert out.read_bytes()==repeat.read_bytes()
    with pytest.raises(ValueError,match="separate destination"):extend(source,source)
