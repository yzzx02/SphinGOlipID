"""Native-CID literature support must distinguish DMS from ordinary water loss."""
from dataclasses import replace
import gzip
from pathlib import Path

import pytest

from sphingolipid_toolkit.engine import IdentificationConfig,SphingolipidIdentifier
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import LibraryRecord,FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import (
    LIBRARY_CACHE_VERSION,normalize_imported_record,_encode_fragment_fields,
    _expand_hydroxy_fa_sphingolipids,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index,IndexedLibrary
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_legacy_sphingolipid import fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from scripts.build_sphingolipid_name_catalog import reviewed_reference_formula
from scripts.extend_verified_modified_sphingoid_library import (
    append_dms,rebuild_verified_index,DMS_CLASS,DMS_NAME,DMS_FORMULA,scientific_payload,
)


def identity():
    import sphingolipid_toolkit._vendor.lipidgate.ms2 as ms2
    return dict(version=LIBRARY_CACHE_VERSION,source_sha256='synthetic',
        rules_sha256=code_fingerprint(Path(ms2.__file__).parent))


def dms_record(record_id=99):
    return LibraryRecord(record_id,DMS_CLASS,DMS_NAME,DMS_NAME,328.321006,'[M+H]+',
        DMS_FORMULA,'+',fragments(DMS_CLASS,328.321006,2))


@pytest.fixture
def dms_index(tmp_path):
    bank=tmp_path/'dms.sqlite'
    write_library_index(bank,[normalize_imported_record(dms_record())],identity())
    return bank


def identify(bank,peaks):
    with SphingolipidIdentifier(IdentificationConfig(library_path=bank)) as identifier:
        return identifier.identify(328.321006,peaks)


@pytest.mark.parametrize('peaks',[
    [(310.3104,100.),(292.2999,80.),(280.2999,70.),(328.321006,20.)],
    [(110.0964,100.),(310.3104,80.),(292.2999,70.)],
    [(58.0651,100.),(84.0808,80.),(310.3104,70.)],
])
def test_dms_requires_110_and_another_distinct_n_fragment(dms_index,peaks):
    assert identify(dms_index,peaks)==[]


@pytest.mark.parametrize('second',[58.0651,84.0808])
def test_two_dms_head_fragments_admit_without_other_gate(dms_index,second):
    rows=identify(dms_index,[(110.0964,100.),(second,2.)])
    assert len(rows)==1 and rows[0]['matched_name']==DMS_NAME
    assert rows[0]['passed_required_gates']
    assert rows[0]['hg_score']==60. and rows[0]['other_score']==0.


def test_dms_uses_strongest_hg_and_existing_strict_other_threshold(dms_index):
    def score(second,other):
        return identify(dms_index,[(110.0964,100.),(84.0808,second),(310.3104,other)])[0]
    low=score(2.,.5);high=score(80.,.5);supported=score(2.,.5001)
    assert low['hg_score']==high['hg_score']==supported['hg_score']==60.
    assert low['other_score']==high['other_score']==0.
    assert supported['other_score']>0.
    assert all(r['lcb_score']==0. for r in [low,high,supported])


def test_weak_110_does_not_turn_a_wrong_head_combination_into_admission(dms_index):
    # The normal global 0.2% preprocessing cutoff remains in force.
    assert identify(dms_index,[(58.0651,100.),(84.0808,90.),(110.0964,.1)])==[]


def test_elemental_assignments_match_the_published_three_decimal_dms_peaks():
    data={f.name:f.mz for f in fragments(DMS_CLASS,328.321006,2)}
    published={'DMS-C3H8N+':58.066,'DMS-C5H10N+':84.082,'DMS-C7H12N+':110.097}
    for name,observed in published.items():
        assert abs(data[name]-observed)/observed*1e6<20.
    assert data['M+H-H2O']==310.3104 and data['M+H-2H2O']==292.2999


def test_lm1_formula_correction_adds_exactly_the_missing_neuac_residue():
    formula,note=reviewed_reference_formula('LM1','C60H110N2O23')
    assert formula=='C71H127N3O31' and 'missing_NeuAc' in note
    assert reviewed_reference_formula('Cer','C34H67NO3')==('C34H67NO3','')
    with pytest.raises(ValueError):reviewed_reference_formula('LM1','C60H110N2O23Na')
    with pytest.raises(ValueError):reviewed_reference_formula('LM1','C71H127N3O31')


def test_general_catalog_extender_cannot_open_an_unvalidated_dms_grid():
    from scripts.extend_positive_catalog_glycosphingolipid_library import legacy_candidates
    rows=[['C20H41NO2','','','N,N-dimethylSo(d18:1)'],
          ['C22H45NO2','','','N,N-dimethylSo(d20:1)']]
    assert len(list(legacy_candidates(rows)))==1


def write_msp(path,records):
    with gzip.open(path,'wt',encoding='utf-8',newline='\n') as handle:
        for r in records:
            handle.write(f'Name: {r.lipid_chain_name}\nPrecursorMZ: {r.precursor_mz}\nPrecursorType: {r.adduct}\nCompoundClass: {r.compound_class}\nFormula: {r.formula}\n')
            handle.write(f'Comment: MS1_name={r.lipid_name};polarity=+\nNum Peaks: {len(r.fragments)}\n')
            for f in r.fragments:handle.write(f'{f.mz} {f.intensity} {_encode_fragment_fields(f)}\n')
            handle.write('\n')


@pytest.mark.parametrize('parent_formula',['C18H37NO2',''])
def test_fresh_rebuild_preserves_ids_derived_records_and_order(tmp_path,parent_formula):
    spb=normalize_imported_record(LibraryRecord(5,'SPB','SPB(d18:1)','SPB(d18:1)',
        300.2897,'[M+H]+',parent_formula,'+',[FragmentRecord(282.2791,'M+H-H2O','LCB碎片')]))
    cer=normalize_imported_record(LibraryRecord(88,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',
        538.5199,'[M+H]+','C34H67NO3','+',[FragmentRecord(264.2686,'LCB-2H2O','LCB碎片')]))
    derived=replace(_expand_hydroxy_fa_sphingolipids([cer])[-1],record_id=999)
    old=tmp_path/'old.sqlite';write_library_index(old,[spb,cer,derived],identity())
    source=tmp_path/'old.msp.gz';write_msp(source,[spb,cer])
    extended=tmp_path/'new.msp.gz';append_dms(source,old,extended)
    destination=tmp_path/'new.sqlite'
    report=rebuild_verified_index(extended,old,destination,allow_dms=True)
    assert report['preserved_original_records']==3
    assert report['new_record_ids']==[1000]
    reader=IndexedLibrary(destination)
    try:
        original=[r for r in reader if r.record_id!=1000]
        assert [r.record_id for r in original]==[5,88,999]
        assert [scientific_payload(r) for r in original]==[scientific_payload(r) for r in [spb,cer,derived]]
    finally:reader.close()


def test_rebuild_refuses_changed_historical_fragments(tmp_path):
    original=normalize_imported_record(dms_record(77))
    old=tmp_path/'old.sqlite';write_library_index(old,[original],identity())
    altered=replace(original,fragments=[replace(f,mz=f.mz+.1) for f in original.fragments])
    source=tmp_path/'changed.msp.gz';write_msp(source,[altered])
    with pytest.raises(ValueError,match='Fresh normalization differs'):
        rebuild_verified_index(source,old,tmp_path/'new.sqlite')


def test_missing_formula_does_not_bypass_parent_precursor_validation(tmp_path):
    wrong=LibraryRecord(5,'SPB','SPB(d18:1)','SPB(d18:1)',100.,'[M+H]+','','+',[])
    old=tmp_path/'old.sqlite';write_library_index(old,[wrong],identity())
    source=tmp_path/'old.msp.gz';write_msp(source,[wrong])
    with pytest.raises(ValueError,match='Missing eligible'):
        append_dms(source,old,tmp_path/'new.msp.gz')


def test_packaged_reference_catalog_retains_corrected_unconfirmed_names_outside_ms2():
    from sphingolipid_toolkit.library_catalog import iter_reference_names,reference_name_candidates
    wanted={'LM1(d34:1)','N,N,N-trimethylSo(d18:1)','N,N-dimethylSo(d18:1)'}
    rows={row['reference_name']:row for row in iter_reference_names() if row['reference_name'] in wanted}
    lm=rows['LM1(d34:1)']
    assert lm['formula']=='C71H127N3O31' and lm['original_formula']=='C60H110N2O23'
    assert abs(float(lm['precursor_mz'])-1518.85263)<1e-5
    assert not lm['canonical_class'] and lm['coverage_status'].startswith('reference_only')
    tms=rows['N,N,N-trimethylSo(d18:1)']
    assert tms['adduct']=='[M]+' and abs(float(tms['precursor_mz'])-342.336656)<1e-5
    assert tms['coverage_status'].startswith('reference_only')
    dms=rows['N,N-dimethylSo(d18:1)']
    assert dms['coverage_status']=='ms2_composition_covered'
    assert dms['canonical_class']==DMS_CLASS
    assert all(row['reference_name']!='LM1(d34:1)' for row in reference_name_candidates(float(lm['precursor_mz'])))


def test_packaged_library_dms_is_searchable_with_the_current_engine():
    from sphingolipid_toolkit.engine import library_path
    with SphingolipidIdentifier(IdentificationConfig(top_n=10)) as identifier:
        rows=identifier.identify(328.321006,[(110.0964,100.),(84.0808,70.),(58.0651,40.),(310.3104,20.),(280.2999,10.)])
    dms=[r for r in rows if r['compound_class']==DMS_CLASS]
    assert len(dms)==1 and dms[0]['matched_name']==DMS_NAME
    assert dms[0]['resolution_level']=='species_level'
    assert dms[0]['engine_revision']=='sphingolipid-literature-20261110-r1'


@pytest.mark.parametrize('tolerance',[{'fragment_ppm':20.},{'fragment_ppm':None,'fragment_da':.3}])
def test_published_dms_ion_positions_pass_while_spb20_control_does_not(dms_index,tolerance):
    # Published m/z positions from Fig. 7b; heights below are synthetic,
    # deliberately not presented as digitized standard-spectrum intensities.
    positive=[(110.097,100.),(84.082,70.),(58.066,40.),(310.311,10.)]
    ordinary=[(310.311,100.),(292.300,80.),(280.300,50.),(95.086,20.),(55.055,30.)]
    with SphingolipidIdentifier(IdentificationConfig(library_path=dms_index,**tolerance)) as identifier:
        assert any(r['compound_class']==DMS_CLASS for r in identifier.identify(328.321,positive))
        assert identifier.identify(328.321,ordinary)==[]



def test_dms_other_threshold_uses_raw_base_peak_even_when_precursor_dominates(dms_index):
    def score(optional):
        peaks=[(328.321006,1000.),(110.0964,100.),(84.0808,70.)]
        if optional is not None:peaks.append((310.3104,optional))
        return identify(dms_index,peaks)[0]
    absent=score(None);equal=score(5.);above=score(5.001)
    assert equal['other_score']==absent['other_score']
    assert above['other_score']>equal['other_score']
