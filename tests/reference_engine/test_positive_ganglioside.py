from dataclasses import replace
import gzip
import importlib.util
from pathlib import Path
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import normalize_imported_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum,LibraryRecord,normalize_peaks,FragmentRecord
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_ganglioside import (
    ADDUCTS,LEGACY_CLASSES as CLASSES,GangliosideComposition,NH3,WATER,from_name,normalize_fragments,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher

@pytest.mark.parametrize('cls,expected',[
    ('GM1',773.945603694331),('GD1',919.493311950501),
    ('GT1',1065.041020206671),('GD3',736.927213977641),
])
def test_charge_and_ammonia_offsets_preserve_protonated_products(cls,expected):
    model=GangliosideComposition(cls,'d',18,1,False,18,0)
    assert model.precursor_mz('[M+2H]2+')==pytest.approx(expected,abs=1e-9)
    assert model.precursor_mz('[M+H+NH4]2+')-expected==pytest.approx(NH3/2)
    assert model.precursor_mz('[M+2NH4]2+')-expected==pytest.approx(NH3)
    assert model.precursor_mz('[M+2H-H2O]2+')-expected==pytest.approx(-WATER/2)
    for adduct in ADDUCTS:
        products={f.name:f for f in model.fragments(adduct)}
        assert products['B-NeuAc'].mz==292.1027
        assert products['B-NeuAc-H2O'].mz==274.0921
        assert products['Y0'].mz==566.5507 and products['Z0'].mz==548.5401
        assert products['LCB-H2O'].mz==282.2791 and products['LCB-2H2O'].mz==264.2686
    assert any(f.mz>expected and f.fragment_type=='Diagnostic_HG' for f in model.fragments('[M+2H]2+'))

@pytest.mark.parametrize('cls',sorted(CLASSES))
@pytest.mark.parametrize('adduct',list(ADDUCTS))
@pytest.mark.parametrize('series',['d','t'])
def test_ganglioside_default_requires_both_fragment_directions_and_shared_one_percent(tmp_path,cls,adduct,series):
    model=GangliosideComposition(cls,series,18,1 if series=='d' else 0,False,18,0)
    record=normalize_imported_record(LibraryRecord(0,cls,model.species_name,model.name,model.precursor_mz(adduct),
        adduct,formula=model.formula,polarity='+',fragments=model.fragments(adduct)))
    bank=tmp_path/'ganglioside.sqlite'
    write_library_index(bank,[record],dict(version=LIBRARY_CACHE_VERSION,
        rules_sha256=code_fingerprint(Path(__file__).resolve().parents[2]/'src/sphingolipid_toolkit/_vendor/lipidgate/ms2')))
    fast=LipidMS2Searcher(bank,native_min_candidates=0,fragment_tolerance_da=.001,min_total_score=0.)
    slow=LipidMS2Searcher(bank,use_native_engine=False,use_score_prefilter=False,fragment_tolerance_da=.001,min_total_score=0.)
    names={f.name:f for f in record.fragments}
    b2='B-NeuAc-H2O' if cls=='GD3' else 'B-HexHexNAc'
    structural=[(names[n].mz,1000.) for n in ('Y0','Z0','B-NeuAc',b2)]
    try:
        for peaks,passed in ((structural,True),(structural[:2],False),(structural[2:],False),
                             (structural[1:],False),
                             ([structural[0],structural[0],*structural[2:]],False),
                             ([*structural[:2],structural[2]],False),
                             ([(structural[0][0],10.),*structural[1:]],True),
                             (structural[1:]+[(f.mz,1000.) for f in record.fragments if f.fragment_type in {'Common','Precursor Ion'}],False),
                             ([(f.mz,1000.) for f in record.fragments if f.fragment_type in {'Common','Precursor Ion'}],False)):
            observed=ExperimentalSpectrum('gate',record.precursor_mz,1.,'+',normalize_peaks(peaks))
            actual=fast._score_sphingo_candidate(observed,record)
            with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
                assert actual==slow._score_sphingo_candidate(observed,record)
            assert actual.passed_required_gates==passed
            assert fast.score_spectrum(observed)==slow.score_spectrum(observed)
            if passed:
                assert actual.resolution_level=='species_level'
                assert actual.pool_scores['hg'].pool_score==50.
                assert actual.pool_scores['b'].pool_score==25.
        if cls!='GD3':
            # Two NeuAc B peaks cannot replace the extended HexNAc B condition.
            peaks=[*structural[:2],(names['B-NeuAc'].mz,1000.),(names['B-NeuAc-H2O'].mz,1000.)]
            assert not fast._score_sphingo_candidate(ExperimentalSpectrum(
                'missing-extended',record.precursor_mz,1.,'+',normalize_peaks(peaks)),record).passed_required_gates
            # Extended glycan evidence without the mandatory NeuAc pair also fails.
            peaks=[*structural[:2],(names['B-HexHexNAc'].mz,1000.),(names['B-HexHexNAc-H2O'].mz,1000.)]
            assert not fast._score_sphingo_candidate(ExperimentalSpectrum(
                'missing-neuac',record.precursor_mz,1.,'+',normalize_peaks(peaks)),record).passed_required_gates
        for label in ('LCB-H2O','LCB-2H2O'):
            for value,level in ((9.999,'species_level'),(10.,'chain_level')):
                observed=ExperimentalSpectrum('confirm',record.precursor_mz,1.,'+',normalize_peaks([
                    *structural,(names[label].mz,value)]))
                actual=fast._score_sphingo_candidate(observed,record)
                assert actual.passed_required_gates and actual.resolution_level==level
                assert fast.score_spectrum(observed)==slow.score_spectrum(observed)
    finally:fast.close();slow.close()

def test_legacy_gm1_peak_identities_masses_and_intensities_survive_role_upgrade():
    old=[FragmentRecord(264.2686,'LCB-2H2O','Common',45.),
         FragmentRecord(282.2791,'LCB-H2O','LCB碎片',51.),
         FragmentRecord(366.1400,'Reference-HG-366.1400','Diagnostic_HG',100.)]
    new=normalize_fragments('GM1','GM1(d18:1/18:0)','[M+2NH4]2+',old)
    assert [(f.name,f.mz,f.intensity) for f in new[:3]]==[(f.name,f.mz,f.intensity) for f in old]
    assert [f.fragment_type for f in new[:3]]==['Common','Common','Diagnostic_B']
    assert not any(f.name=='B-HexHexNAc' for f in new)
    assert normalize_fragments('GM1','GM1(d18:1/18:0)','[M+2NH4]2+',new)==new

def test_glycan_family_source_does_not_claim_a_linkage_isomer():
    for cls in CLASSES:
        model=GangliosideComposition(cls,'t',18,0,False,24,0)
        assert model.name.startswith(cls+'(')
        assert from_name(model.name)==model
        assert not from_name(model.name,cls+'a')

def test_curator_preserves_unrelated_blocks_and_existing_peak_data(tmp_path):
    from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import iter_blocks,parse_peaks,header_value
    script=Path(__file__).resolve().parents[2]/'scripts/curate_positive_ganglioside_library.py'
    spec=importlib.util.spec_from_file_location('ganglioside_curator',script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def block(cls,name,adduct,peaks):
        return [f'Name: {name}','PrecursorMZ: 790.9722',f'PrecursorType: {adduct}',f'CompoundClass: {cls}',
                f'Num Peaks: {len(peaks)}',*peaks]
    blocks=[block('GM3','GM3(d18:1/18:0)','[M+H]+',['292.1027 100.00 "sialic" "Common"']),
            block('Cer','Cer(t18:0/18:0)','[M+H]+',['300.2897 90.00 "LCB-H2O" "LCB碎片"']),
            block('GM1','GM1(d18:1/18:0)','[M+2NH4]2+',[
                '264.2686 45.00 "LCB-2H2O" "Common"','282.2791 51.00 "LCB-H2O" "LCB碎片"',
                '366.1400 100.00 "Reference-HG-366.1400" "Diagnostic_HG"']),
            block('GM1','GM1(d18:1/18:0)','[M-H]-',['366.1400 100.00 "original negative" "Common"']),
            block('PC','PC(18:0/20:6)','[M+H]+',['184.0733 100.00 "original PC" "Diagnostic_HG"'])]
    source=tmp_path/'source.msp.gz';destination=tmp_path/'curated.msp.gz'
    with gzip.open(source,'wt',encoding='utf-8') as output:output.write('\n\n'.join('\n'.join(b) for b in blocks)+'\n\n')
    audit=module.curate(source,destination)
    revised=list(iter_blocks(destination))
    for index in (0,1,3,4):assert revised[index]==blocks[index]
    old,new=parse_peaks(blocks[2]),parse_peaks(revised[2])
    assert [{k:p[k] for k in ('name','mz','intensity')} for p in new[:3]]==[{k:p[k] for k in ('name','mz','intensity')} for p in old]
    assert audit['added_records']==55
    assert not any('/12:' in header_value(b,'Name') for b in revised[len(blocks):])
    assert module.curated_block(revised[2])==revised[2]
    with pytest.raises(ValueError,match='separate destination'):module.curate(source,source)
