from dataclasses import replace
import json
from pathlib import Path

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import (
    BLOCK_SIZE, MAX_CACHED_BLOCKS, IndexedLibrary, validate_precursor_range,
    write_library_index, repack_library_index,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION, load_standard_msp, normalize_imported_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_gb4 import B2_NAME, PositiveGb4
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher


def identity():
    import sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library as module
    return dict(version=LIBRARY_CACHE_VERSION, source_sha256="fixture",
                rules_sha256=code_fingerprint(Path(module.__file__).parent))


def record(index=0, mz=500, lipid_class="PE-P", adduct="[M+H]+"):
    return LibraryRecord(index, lipid_class, "PE(P-38:5)", "PE(P-18:0/20:5)", mz, adduct,
                         polarity="+", fragments=[FragmentRecord(294.3155, "Common", "Common")])


def test_stable_ties_filtering_and_exact_fields(tmp_path):
    records = [record(10, 750), record(20, 500, "LPS", "[M+NH4]+"),
               record(30, 500, "PE-O"), record(40, 500), record(50, 500, adduct="[M+Na]+")]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    indexed = IndexedLibrary(bank)
    assert list(indexed) == [records[2], records[3], records[4], records[0]]
    assert indexed[::-1] == [records[0], records[4], records[3], records[2]]
    assert indexed[-1] == records[0]
    assert list(indexed.precursors) == [500, 500, 500, 750]
    assert "LPS" in indexed.available_classes
    indexed.close()
    selected = IndexedLibrary(bank, allowed_class_keys={"PEP"}, allowed_adducts={"[M+H]+"}, mz_min=500, mz_max=600)
    assert list(selected) == [records[3]]
    selected.close()
    absent = IndexedLibrary(bank, allowed_class_keys={"DOESNOTEXIST"})
    assert len(absent) == 0 and absent._blocks == {}
    absent.close()


def test_candidate_cache_is_bounded_and_lazily_decoded(tmp_path):
    records = [record(i, 500+i/1000) for i in range(BLOCK_SIZE*(MAX_CACHED_BLOCKS+2))]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    indexed = IndexedLibrary(bank)
    assert len(indexed) == len(records)
    assert not indexed._blocks
    assert indexed[0] == records[0]
    assert len(indexed._blocks[0].records) == 1
    for index in range(0, len(records), BLOCK_SIZE):
        assert indexed[index] == records[index]
        assert len(indexed._blocks) <= MAX_CACHED_BLOCKS
    assert 0 not in indexed._blocks
    assert indexed[0] == records[0]
    assert indexed[1] is indexed[1]
    indexed.close()
    assert not indexed._blocks


@pytest.mark.parametrize("bad", [(0, None), (float("nan"), 1000), (1000, 500), (100, float("inf"))])
def test_invalid_precursor_range_is_rejected(bad):
    with pytest.raises(ValueError, match="m/z"):
        validate_precursor_range(*bad)


@pytest.mark.parametrize("unit", ["ppm", "da"])
def test_range_boundary_preserves_candidates_inside_mass_tolerance(tmp_path, unit):
    bank = tmp_path / "library.sqlite"
    records = [record(0, 499.999), record(1, 500), record(2, 600.001), record(3, 700)]
    write_library_index(bank, records, identity())
    options = dict(precursor_tolerance_ppm=5, precursor_tolerance_da=.005 if unit == "da" else None)
    full = LipidMS2Searcher(bank, **options)
    selected = LipidMS2Searcher(bank, precursor_mz_min=500, precursor_mz_max=600, **options)
    assert selected.find_candidates(500) == full.find_candidates(500)
    assert selected.find_candidates(600) == full.find_candidates(600)
    assert selected.find_candidates(499.999) == []
    assert selected.find_candidates(600.001) == []
    full.close()
    selected.close()


def test_disk_search_preserves_pep_support_and_tie_order(tmp_path, monkeypatch):
    import sphingolipid_toolkit._vendor.lipidgate.ms2.search as search_module
    pe = normalize_imported_record(replace(record(1, 750.5432), fragments=[
        FragmentRecord(285.2213, "(R=O)+(20:5)", "FA_Frag"),
        FragmentRecord(359.2581, "[M-C2H6O3NP-(RCH2=CH-OH)+H]+(P-18:0)", "Diagnostic_FA_Loss"),
        FragmentRecord(392.2924, "[M-C3H4O-(ROOH)+H]+(20:5)", "Diagnostic_FA_Loss"),
        FragmentRecord(609.5241, "[M-C2H8O4NP+H]+", "Diagnostic_HG"),
        FragmentRecord(750.5432, "[M+H]+", "Precursor Ion"),
    ]))
    records = [pe, replace(pe, record_id=2), record(3, 900)]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    monkeypatch.setattr(search_module, "load_library", lambda _: records)
    options = dict(min_total_score=0, fragment_tolerance_da=None, fragment_tolerance_ppm=15)
    full = LipidMS2Searcher(tmp_path / "fixture.msp", **options)
    indexed = LipidMS2Searcher(bank, **options)
    spectrum = ExperimentalSpectrum("scan_3768", 750.5440323345842, 10.60465, "+", normalize_peaks([
        (359.25775869833893, 1285.9898681640625), (392.2947607640448, 562.2628173828125),
        (609.5304463578401, 362.0347595214844), (294.31788640500923, 110.01720428466797),
    ]), precursor_charge=1)
    expected = full.score_spectrum(spectrum, top_n=3)
    actual = indexed.score_spectrum(spectrum, top_n=3)
    assert actual == expected
    assert actual[0]["final_score"] == 87.8328
    assert actual[0]["downgrade_reason"] == ""
    # At 10 ppm the measured HG falls outside tolerance: the explicit PE-P
    # dual gate must reject it identically in both storage paths.
    full.fragment_tolerance_ppm = indexed.fragment_tolerance_ppm = 10
    assert full.score_spectrum(spectrum, top_n=3) == indexed.score_spectrum(spectrum, top_n=3) == []
    full.close()
    indexed.close()


def test_parallel_capacity_uses_bounded_index_budget(tmp_path, monkeypatch):
    from sphingolipid_toolkit._vendor.lipidgate.ms2 import file_search
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, [record()], identity())
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 8 * 1024**3)
    assert file_search._estimated_worker_memory_bytes(bank) == 768 * 1024**2
    assert file_search.validate_parallel_capacity(4, 4, bank) == 4
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 4 * 1024**3)
    with pytest.raises(MemoryError, match="请减少进程数"):
        file_search.validate_parallel_capacity(4, 4, bank)


@pytest.mark.parametrize("lipid_class", ["Hex2Cer", "LacCer"])
def test_dihexosylceramide_index_keeps_hg_only_sum_composition(tmp_path, lipid_class):
    name = f"{lipid_class}(d18:1/16:0)"
    candidate = normalize_imported_record(LibraryRecord(0, lipid_class, name, name, 862.625, "[M+H]+", fragments=[
        FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片"),
        FragmentRecord(282.2791, "LCB-H2O", "LCB碎片"),
        FragmentRecord(502.4982, "M+H-C12H20O10-2H2O", "Diagnostic_HG"),
        FragmentRecord(520.5088, "M+H-C12H20O10-H2O", "Diagnostic_HG"),
        FragmentRecord(538.5194, "M+H-C12H20O10", "Diagnostic_HG"),
        FragmentRecord(556.5299, "M+H-C12H18O9", "Diagnostic_HG"),
        FragmentRecord(664.5511, "M+H-C6H12O6-H2O", "Diagnostic_HG"),
        FragmentRecord(682.5616, "M+H-C6H12O6", "Diagnostic_HG"),
        FragmentRecord(700.5722, "M+H-C6H10O5", "Diagnostic_HG"),
        FragmentRecord(844.6145, "M+H-H2O", "C类碎片"),
    ], polarity="+"))
    assert all(fragment.mz != 664.5511 for fragment in candidate.fragments)
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, [candidate], identity())
    fast = LipidMS2Searcher(bank)
    reference = LipidMS2Searcher(bank, use_native_engine=False, use_score_prefilter=False)
    try:
        hg = [(502.4982, 1000), (520.5088, 1000), (682.5616, 1000), (844.6145, 1000)]
        for peaks, expected in (
            (hg[2:] + [(664.5511, 1000), (264.2686, 1000)], None),
            (hg[1:] + [(664.5511, 1000), (264.2686, 1000)], name),
            (hg, f"{lipid_class}(d34:1)"),
            (hg + [(264.2686, 10)], name),
        ):
            spectrum = ExperimentalSpectrum("dihexosyl", 862.625, 1.0, "+", normalize_peaks(peaks))
            rows = fast.score_spectrum(spectrum)
            assert rows == reference.score_spectrum(spectrum)
            assert [row["matched_name"] for row in rows] == ([expected] if expected else [])
            if expected:
                assert rows[0]["hg_score"] == 60.0
                assert rows[0]["passed_required_gates"]
    finally:
        fast.close()
        reference.close()


def test_hex3cer_index_uses_half_yz_gate_and_hex2cer_lcb_reporting(tmp_path):
    record = normalize_imported_record(LibraryRecord(
        0, "Hex3Cer", "Hex3Cer(d32:1)", "Hex3Cer(d18:1/14:0)", 996.6465, "[M+H]+", polarity="+",
        fragments=[
            FragmentRecord(252.2686, "LCB-CH2O-H2O", "LCB碎片"),
            FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片"),
            FragmentRecord(282.2791, "LCB-H2O", "LCB碎片"),
            FragmentRecord(474.4669, "M+H-3Hex-2H2O", "Diagnostic_HG"),
            FragmentRecord(492.4775, "M+H-3Hex-H2O", "Common"),
            FragmentRecord(510.4881, "M+H-3Hex", "Diagnostic_HG"),
            FragmentRecord(636.5198, "M+H-2Hex-2H2O", "Common"),
            FragmentRecord(654.5303, "M+H-2Hex-H2O", "Common"),
            FragmentRecord(672.5409, "M+H-2Hex", "Common"),
            FragmentRecord(798.5726, "M+H-1Hex-2H2O", "Common"),
            FragmentRecord(816.5832, "M+H-1Hex-H2O", "Common"),
            FragmentRecord(834.5937, "M+H-1Hex", "Common"),
            FragmentRecord(978.6360, "M+H-H2O", "Common"),
            FragmentRecord(996.6465, "[M+H]+", "Precursor Ion"),
        ],
    ))
    assert sum(f.fragment_type == "Diagnostic_HG" for f in record.fragments) == 6
    assert [f.name for f in record.fragments if f.fragment_type == "LCB碎片"] == ["LCB-2H2O"]
    assert [f.name for f in record.fragments if f.name.startswith("LCB")] == ["LCB-2H2O", "LCB-H2O"]
    assert all(f.name not in {"M+H-1Hex-2H2O", "M+H-2Hex-2H2O"} for f in record.fragments)
    bank = tmp_path / "hex3cer.sqlite"
    write_library_index(bank, [record], identity())
    fast = LipidMS2Searcher(bank, native_min_candidates=0)
    reference = LipidMS2Searcher(bank, use_native_engine=False, use_score_prefilter=False)
    try:
        hg = [(834.5937, 1000), (816.5832, 1000), (672.5409, 50)]
        precursor, water_loss = (996.6465, 1000), (978.6360, 1000)
        z0_water = (474.4669, 1000)
        for peaks, expected, other_score in (
            (hg[:2] + [precursor, water_loss, z0_water, (798.5726, 1000), (264.2686, 1000)], None, None),
            (hg, "Hex3Cer(d32:1)", 0.),
            (hg + [z0_water], "Hex3Cer(d32:1)", 20. / 3),
            (hg + [(264.2686, 9.9)], "Hex3Cer(d32:1)", 0.),
            (hg + [(264.2686, 10)], "Hex3Cer(14:0/d18:1)", 0.),
        ):
            spectrum = ExperimentalSpectrum("hex3cer", record.precursor_mz, 1.0, "+", normalize_peaks(peaks))
            rows = fast.score_spectrum(spectrum)
            assert rows == reference.score_spectrum(spectrum)
            assert [row["matched_name"] for row in rows] == ([expected] if expected else [])
            if expected:
                [row] = rows
                assert row["hg_score"] == 60. and row["other_score"] == pytest.approx(other_score)
                assert (row["lcb_score"] > 0.) == any(mz == 264.2686 for mz, _ in peaks)
                assert row["passed_required_gates"]
                reason = "missing_lcb_chain_evidence" if "/" not in expected else ""
                assert row["downgrade_reason"] == reason
                evidence = {f["name"]: f for f in json.loads(row["ms2_evidence_json"])["fragments"]}
                assert evidence["M+H-3Hex-2H2O"]["scoring_pool"] == "other"
                assert evidence["LCB-2H2O"]["scoring_pool"] == "lcb"
    finally:
        fast.close()
        reference.close()


def test_gb4_masses_other_scoring_and_soft_lcb_confirmation(tmp_path):
    from scripts.curate_positive_gb4_library import gb4_block

    model = PositiveGb4(18, 2, 16, 0)
    assert model.formula == "C60H108N2O23"
    assert round(model.precursor_mz, 4) == 1225.7416
    assert round(PositiveGb4(18, 1, 14, 0).precursor_mz, 4) == 1199.7259
    fragments = model.fragments()
    masses = {f.name: round(f.mz, 4) for f in fragments}
    assert masses[B2_NAME] == 366.1395
    assert masses["M+H-HexNAc"] == 1022.6622
    assert masses["M+H-HexNAc-3Hex-2H2O"] == 500.4826
    assert masses["LCB-2H2O"] == 262.2529
    assert sum(f.fragment_type == "Diagnostic_HG" for f in fragments) == 9
    assert {f.name: f.fragment_type for f in fragments if f.name in {B2_NAME, "LCB-2H2O"}} == {
        B2_NAME: "Common", "LCB-2H2O": "Common",
    }

    source = tmp_path / "gb4.msp"
    source.write_text("\n".join(gb4_block(model)) + "\n", encoding="utf-8")
    [candidate] = load_standard_msp(source)
    assert (candidate.lipid_name, candidate.polarity) == ("Gb4(d34:2)", "+")
    bank = tmp_path / "gb4.sqlite"
    write_library_index(bank, [candidate], identity())
    fast = LipidMS2Searcher(bank, native_min_candidates=0)
    reference = LipidMS2Searcher(bank, use_native_engine=False, use_score_prefilter=False)
    try:
        hg = [(masses[name], 1000) for name in (
            "M+H-HexNAc", "M+H-HexNAc-H2O", "M+H-HexNAc-1Hex", "M+H-HexNAc-1Hex-H2O", "M+H-HexNAc-2Hex",
        )]
        other = [(masses[name], 1000) for name in ("[M+H]+", "M+H-H2O", B2_NAME)]
        lcb = masses["LCB-2H2O"]
        for peaks, expected, other_score in (
            (hg[:4] + other + [(lcb, 1000)], None, None),
            ([(mz, 1) for mz, _ in hg] + other, None, None),
            (hg, "Gb4(d34:2)", 0.),
            (hg + [(lcb, 9.9)], "Gb4(d34:2)", 25. / 3),
            (hg + [(lcb, 10)], "Gb4(d18:2/16:0)", 25. / 3),
            (hg + other + [(lcb, 10)], "Gb4(d18:2/16:0)", 25.),
        ):
            spectrum = ExperimentalSpectrum("gb4", candidate.precursor_mz, 1.0, "+", normalize_peaks(peaks))
            rows = fast.score_spectrum(spectrum)
            assert rows == reference.score_spectrum(spectrum)
            assert [row["matched_name"] for row in rows] == ([expected] if expected else [])
            if expected:
                [row] = rows
                assert row["passed_required_gates"] and row["hg_score"] == 75.
                assert row["other_score"] == pytest.approx(other_score) and row["lcb_score"] == 0.
                reason = "missing_lcb_chain_evidence" if "/" not in expected else ""
                assert row["downgrade_reason"] == reason
                evidence = {f["name"]: f for f in json.loads(row["ms2_evidence_json"])["fragments"]}
                assert evidence["LCB-2H2O"]["scoring_pool"] == "other"
                assert evidence["LCB-2H2O"]["role"] == "chain"
    finally:
        fast.close()
        reference.close()


def test_fragment_templates_preserve_occurrences_types_and_provenance(tmp_path):
    import math
    import pickle
    import sqlite3
    import zlib

    repeated = FragmentRecord(294.3155, "Common", "Common", 100, 1, None)
    positive_zero = FragmentRecord(300.0, "zero", "Common", 0.0)
    negative_zero = replace(positive_zero, intensity=-0.0)
    original = replace(record(), fragments=[repeated, repeated, positive_zero, negative_zero],
                       metadata={"source_chain_name": "PE(P-18:0/20:5)", "custom": "audit"})
    bank = tmp_path / "templates.sqlite"
    write_library_index(bank, [original, replace(original, record_id=2, metadata=dict(original.metadata))], identity())
    with sqlite3.connect(bank) as connection:
        payload = connection.execute("SELECT payload FROM blocks WHERE block_id=0").fetchone()[0]
        templates, values, columns = pickle.loads(zlib.decompress(payload))
        assert len(templates) == 3
        assert columns["occurrences"][0] == columns["occurrences"][1]
        assert list(columns["counts"]) == [4, 4]
    indexed = IndexedLibrary(bank)
    first, second = indexed[0], indexed[1]
    assert first == original
    assert first.fragments[0] is not first.fragments[1]
    assert first.fragments[0] is not second.fragments[0]
    assert type(first.fragments[0].intensity) is int
    assert math.copysign(1, first.fragments[2].intensity) == 1
    assert math.copysign(1, first.fragments[3].intensity) == -1
    first.metadata["custom"] = "edited"
    assert second.metadata["custom"] == "audit"
    indexed.close()


def test_stream_repack_preserves_legacy_data_and_rejects_wrong_identity(tmp_path):
    import json
    import pickle
    import sqlite3
    import zlib
    from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import _pack_record

    originals = [record(1, 500), replace(record(2, 500), metadata={"origin": "fixture"})]
    source = tmp_path / "old.sqlite"
    expected = dict(version=LIBRARY_CACHE_VERSION, source_sha256="fixture", rules_sha256="old")
    metadata = dict(expected, index_format=1, record_encoding="tuples-v1", block_size=1024,
                    record_count=2, classes=["PE-P"], adducts=["[M+H]+"])
    with sqlite3.connect(source) as connection:
        connection.executescript("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT);"
                                 "CREATE TABLE blocks (block_id INTEGER PRIMARY KEY, payload BLOB);"
                                 "CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, precursor REAL, class_id INTEGER, adduct_id INTEGER);")
        connection.execute("INSERT INTO metadata VALUES ('library', ?)", (json.dumps(metadata),))
        connection.execute("INSERT INTO blocks VALUES (0, ?)",
                           (zlib.compress(pickle.dumps([_pack_record(r) for r in originals])),))
        connection.executemany("INSERT INTO entries VALUES (?, ?, 0, 0)", [(0, 500), (1, 500)])
    legacy = IndexedLibrary(source)
    assert list(legacy) == originals
    legacy.close()
    destination = tmp_path / "new.sqlite"
    migrated = repack_library_index(source, destination, identity(), expected_identity=expected)
    assert migrated["normalized_from_rules_sha256"] == "old"
    assert len(migrated["normalized_from_index_sha256"]) == 64
    indexed = IndexedLibrary(destination)
    assert list(indexed) == originals
    indexed.close()
    with pytest.raises(ValueError, match="identity"):
        repack_library_index(source, tmp_path / "wrong.sqlite", identity(),
                             expected_identity=dict(expected, source_sha256="different"))
    with pytest.raises(ValueError, match="version or source"):
        repack_library_index(source, tmp_path / "wrong_version.sqlite", dict(identity(), version=999),
                             expected_identity=expected)
    with sqlite3.connect(source) as connection:
        connection.execute("UPDATE entries SET precursor=600 WHERE ordinal=1")
    with pytest.raises(ValueError, match="precursor index"):
        repack_library_index(source, destination, identity(), expected_identity=expected)
    indexed = IndexedLibrary(destination)
    assert list(indexed) == originals
    indexed.close()
