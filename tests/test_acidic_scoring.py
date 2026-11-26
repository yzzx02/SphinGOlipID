import json

import pytest

from sphingolipid_toolkit.engine import IdentificationConfig, SphingolipidIdentifier

from sphingolipid_toolkit.acidic_scoring import AcidicSearcher, SCORING_PROFILES, normalize_acidic_weights
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION, normalize_imported_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import (
    ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_gm3 import PositiveGm3
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_ganglioside import GangliosideComposition
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_gd2 import PositiveGd2
from sphingolipid_toolkit._vendor.lipidgate.ms2 import positive_catalog_glycosphingolipid as glycan_catalog
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2 import indexed_library
from pathlib import Path


def write_record(path, name, adduct, model, fragments, precursor):
    record = normalize_imported_record(LibraryRecord(0, name, model.species_name, model.name, precursor,
                                                      adduct, formula=model.formula, polarity='+', fragments=fragments))
    write_library_index(path, [record], {'version': LIBRARY_CACHE_VERSION,
                                         'rules_sha256': code_fingerprint(Path(indexed_library.__file__).parent)})
    return record


def test_protonated_gm3_pools_and_weights(tmp_path):
    bank = tmp_path / 'gm3.sqlite'
    model = PositiveGm3(18, 1, 18, 0)
    record = write_record(bank, 'GM3', '[M+H]+', model, model.fragments(), model.precursor_mz)
    observed = [(548.5374455047629, 100.), (872.6489713745731, 59.84195194405749),
                (710.5920102500127, 45.), (292.1010729528085, 30.), (274.0923934275803, 8.7),
                (264.26771648055944, 12.603661540045552), (1163.7397576881472, 42.9),
                (454.15537086165386, 8.65), (1181.7469575306252, 22.1), (530.5336127593374, 1.)]
    spectrum = ExperimentalSpectrum('958592', record.precursor_mz, 15.9764, '+', normalize_peaks(observed))
    search = AcidicSearcher(bank, acidic_weights={}, fragment_tolerance_ppm=20)
    try:
        row, = search.score_spectrum(spectrum)
        weighted = json.loads(row['ms2_evidence_json'])
        assert row['matched_name'] == model.name
        assert row['final_score'] == pytest.approx(93.6956, abs=0.0002)
        assert weighted['score_profile']['weights_pct'] == {'hg': 60, 'lcb': 20, 'other': 20}
        assert {f['name'] for f in weighted['fragments'] if f['scoring_pool'] == 'hg'} == {
            f.name for f in model.fragments() if f.fragment_type == 'Diagnostic_HG'
        } | {'[Neu5Ac+H-H2O]+', '[Neu5Ac+H-2H2O]+'}
        assert any(f['name'] == 'Z0-H2O' and f['observed_mz'] for f in weighted['fragments'])
        search.acidic_weights = normalize_acidic_weights({'GM3': (8, 1, 1)})
        changed, = search.score_spectrum(spectrum)
        quality = weighted['score_profile']['pool_scores']
        assert changed['final_score'] == pytest.approx(80*quality['hg']+10*quality['lcb']+10*quality['other'], abs=.00006)
        rejected = ExperimentalSpectrum('loss', record.precursor_mz, 1., '+', normalize_peaks(
            [pair for pair in observed if pair[0] not in (292.1010729528085, 274.0923934275803)]))
        assert search.score_spectrum(rejected) == []
    finally:
        search.close()


def test_gm3_half_headgroup_boundary_and_separate_chain_confirmation(tmp_path):
    bank = tmp_path / 'gm3.sqlite'
    model = PositiveGm3(18, 1, 18, 0)
    record = write_record(bank, 'GM3', '[M+H]+', model, model.fragments(), model.precursor_mz)
    yz = [f.mz for f in record.fragments if f.fragment_type == 'Diagnostic_HG']
    lcb = next(f.mz for f in record.fragments if f.name == 'LCB-2H2O')
    precursor = record.precursor_mz
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank)) as engine:
        def identify(count, lcb_percent):
            peaks = [(mz, 100.) for mz in yz[:count]] + [(precursor, 99.), (lcb, lcb_percent)]
            return engine.identify(precursor, peaks, precursor_charge=1)
        assert not identify(3, 99.)
        total, = identify(4, .99)
        chain, = identify(4, 1.)
        gate = json.loads(total['ms2_evidence_json'])['hg_admission']
        assert gate == {'matched_channels': 4, 'library_channels': 8,
                        'required_channels': 4, 'minimum_fraction': .5}
        assert total['resolution_level'] == 'species_level'
        assert chain['resolution_level'] == 'chain_level'
        assert chain['matched_name'] == model.name


@pytest.mark.parametrize('cls', ['GM1', 'GD1', 'GT1', 'GD3'])
@pytest.mark.parametrize('adduct', ['[M+H]+', '[M+2NH4]2+'])
def test_ganglioside_headgroup_and_lcb_scores(tmp_path, cls, adduct):
    assert SCORING_PROFILES[cls] == (('hg', 'lcb', 'other'), (60, 20, 20))
    bank = tmp_path / 'ganglioside.sqlite'
    model = GangliosideComposition(cls, 'd', 16, 1, False, 18, 0)
    record = write_record(bank, cls, adduct, model, model.fragments(adduct), model.precursor_mz(adduct))
    products = {f.name: f.mz for f in record.fragments}
    gates = {'Y0': 40., 'Z0': 20., 'B-NeuAc': 100.,
             'B-NeuAc-H2O' if cls == 'GD3' else 'B-HexHexNAc': 30.}
    def identify(engine, ions):
        peaks = [(products[name], intensity) for name, intensity in ions.items()]
        return engine.identify(record.precursor_mz, peaks, precursor_charge=2 if '2+' in adduct else 1)

    with SphingolipidIdentifier(IdentificationConfig(library_path=bank)) as engine:
        total, = identify(engine, gates)
        assert total['passed_required_gates'] and total['resolution_level'] == 'species_level'
        assert (total['hg_score'], total['lcb_score'], total['other_score'], total['final_score']) == (60., 0., 0., 60.)
        assert 'b_score' not in total
        for missing in ('Y0', 'Z0', 'B-NeuAc'):
            ions = {key: val for key, val in gates.items() if key != missing}
            assert not identify(engine, ions | {'LCB-H2O': 100.})
        if cls != 'GD3':
            no_extended_hexnac = {'Y0': 40., 'Z0': 20., 'B-NeuAc': 100., 'B-NeuAc-H2O': 30.}
            assert not identify(engine, no_extended_hexnac)

        below, = identify(engine, gates | {'LCB-2H2O': .99, 'Z0-H2O': .5})
        assert below['final_score'] == total['final_score']
        chain, = identify(engine, gates | {'LCB-2H2O': 1., 'Z0-H2O': .51})
        assert chain['lcb_score'] == pytest.approx(20 * 1.05 * .01 / (.01 + .05))
        assert chain['resolution_level'] == 'chain_level'
        assert chain['other_score'] == pytest.approx(20 / 3)
        assert chain['final_score'] == pytest.approx(60 + chain['lcb_score'] + chain['other_score'], abs=.00006)
        audit = json.loads(chain['ms2_evidence_json'])
        assert audit['score_profile']['weights_pct'] == {'hg': 60, 'lcb': 20, 'other': 20}
        fragments = {f['name']: f for f in audit['fragments']}
        assert {fragments[name]['scoring_pool'] for name in gates} == {'hg'}
        assert fragments['B-NeuAc']['role'] == 'headgroup'
        assert fragments['LCB-2H2O']['role'] == 'chain'
        assert {f['name'] for f in audit['fragments'] if f['contributes_to_score']} == {
            'B-NeuAc', 'LCB-2H2O', 'Z0-H2O'}

    weights = {cls: (65, 35, 0)}
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank, acidic_weights=weights)) as engine:
        row, = identify(engine, gates | {'LCB-H2O': 100.})
        assert (row['hg_score'], row['lcb_score'], row['final_score']) == (65., 35., 100.)


def test_positive_profiles_use_no_separate_glycan_b_score():
    assert all('b' not in pools for pools, _ in SCORING_PROFILES.values())
    for cls in glycan_catalog.CLASSES | {'GD2', 'Gal-Fuc-GM1', 'Gal-Fuc-GD1'}:
        assert SCORING_PROFILES[cls] == (('hg', 'lcb', 'other'), (60, 20, 20))


@pytest.mark.parametrize('cls', ['GD2', 'Gal-Fuc-GM1', 'Gal-Fuc-GD1'])
def test_glycan_b_ions_cannot_supply_lcb_points(tmp_path, cls):
    bank = tmp_path / 'glycan.sqlite'
    adduct = '[M+H]+'
    if cls == 'GD2':
        model = PositiveGd2(PositiveGm3(18, 1, 18, 0))
        precursor, products = model.precursor_mz, model.fragments()
    else:
        model = GangliosideComposition(cls, 'd', 18, 1, False, 18, 0)
        precursor, products = model.precursor_mz(adduct), model.fragments(adduct)
    record = write_record(bank, cls, adduct, model, products, precursor)
    head = [(f.mz, 100.) for f in record.fragments if f.fragment_type in {'Diagnostic_HG', 'Diagnostic_B'}]
    lcb = next(f.mz for f in record.fragments if f.name == 'LCB-2H2O')
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank)) as engine:
        first, = engine.identify(precursor, head)
        second, = engine.identify(precursor, head + [(lcb, 100.)])
        assert (first['hg_score'], first['lcb_score'], first['other_score']) == (60., 0., 0.)
        assert (second['hg_score'], second['lcb_score'], second['other_score']) == (60., 20., 0.)
        fragments = json.loads(second['ms2_evidence_json'])['fragments']
        assert all(f['scoring_pool'] == 'hg' for f in fragments if f['fragment_type'] == 'Diagnostic_B')
        assert 'b_score' not in second


def test_sum_composition_without_lcb_channels_keeps_zero_lcb_points(tmp_path):
    bank = tmp_path / 'sum.sqlite'
    precursor = 1186.731
    name = 'Hex4Cer(d34:1)'
    record = normalize_imported_record(LibraryRecord(
        0, 'Hex4Cer', name, name, precursor, '[M+H]+', 'C58H107NO23', '+',
        glycan_catalog.fragments('Hex4Cer', precursor)))
    write_library_index(bank, [record], {'version': LIBRARY_CACHE_VERSION,
        'rules_sha256': code_fingerprint(Path(indexed_library.__file__).parent)})
    observed = [(f.mz, 100.) for f in record.fragments]
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank)) as engine:
        row, = engine.identify(precursor, observed + [(264.2686, 100.)])
        assert (row['hg_score'], row['lcb_score'], row['other_score'], row['final_score']) == (60., 0., 20., 80.)
        assert row['resolution_level'] == 'species_level'
        assert 'b_score' not in row
        audit = json.loads(row['ms2_evidence_json'])
        assert audit['score_profile']['weights_pct'] == {'hg': 60, 'lcb': 20, 'other': 20}
        assert not any(f['scoring_pool'] == 'lcb' for f in audit['fragments'])
        no_b = [(f.mz, 100.) for f in record.fragments if f.fragment_type != 'Diagnostic_B']
        assert not engine.identify(precursor, no_b)


def test_sm_phosphocholine_loss_scores_only_in_other(tmp_path):
    bank = tmp_path / 'sm.sqlite'
    precursor_mz = 649.5279
    loss_mz = precursor_mz - 183.06604448
    record = normalize_imported_record(LibraryRecord(
        0, 'SM', 'SM(d30:0)', 'SM(d14:0/16:0)', precursor_mz, '[M+H]+', polarity='+',
        fragments=[FragmentRecord(precursor_mz, '[M+H]+', 'Precursor Ion'),
                   FragmentRecord(184.0733, '[C5H15NO4P]+', 'Diagnostic_HG'),
                   FragmentRecord(210.2216, 'LCB-2H2O', 'LCB碎片')],
    ))
    write_library_index(bank, [record], {'version': LIBRARY_CACHE_VERSION,
                                         'rules_sha256': code_fingerprint(Path(indexed_library.__file__).parent)})
    primary_peaks = [(184.0733, 1000.), (210.2216, 800.)]
    search = AcidicSearcher(bank, acidic_weights={}, fragment_tolerance_ppm=20)
    try:
        def score(peaks):
            spectrum = ExperimentalSpectrum('sm', precursor_mz, 1., '+', normalize_peaks(peaks))
            return search.score_spectrum(spectrum)

        reference, = score(primary_peaks)
        row, = score([*primary_peaks, (loss_mz, 50.)])
        evidence = json.loads(row['ms2_evidence_json'])
        loss, = (f for f in evidence['fragments'] if f['name'] == 'M+H-183')
        assert loss['fragment_type'] == 'NL'
        assert loss['scoring_pool'] == 'other'
        assert loss['contributes_to_score'] is True
        assert loss['observed_mz'] == pytest.approx(loss_mz)
        assert row['matched_fragment_count'] == reference['matched_fragment_count'] + 1
        assert row['hg_score'] == reference['hg_score']
        assert row['resolution_level'] == reference['resolution_level']
        assert row['other_score'] - reference['other_score'] == pytest.approx(20 / 3)
        assert row['lcb_score'] == reference['lcb_score'] > 0

        boundary, = score([*primary_peaks, (loss_mz, 5.)])
        outside_window, = score([*primary_peaks, (loss_mz + .03, 50.)])
        assert boundary['other_score'] == outside_window['other_score'] == reference['other_score']
        assert score([(210.2216, 800.), (loss_mz, 50.)]) == []
    finally:
        search.close()


@pytest.mark.parametrize('weights', [(0, 0, 0), (-1, 2, 3), (1, 2), (1, 2, float('nan'))])
def test_invalid_pool_weights(weights):
    with pytest.raises(ValueError):
        normalize_acidic_weights({'GM3': weights})
