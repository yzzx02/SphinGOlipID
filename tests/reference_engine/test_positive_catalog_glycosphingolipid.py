import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.library import write_standard_msp
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_glycosphingolipid import fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher


def catalog_record():
    return LibraryRecord(0,'Hex4Cer','Hex4Cer(d34:1)','Hex4Cer(d34:1)',1186.731,
        '[M+H]+','C58H107NO23','+',fragments('Hex4Cer',1186.731))


def run_record(tmp_path, removed=(), other_intensity=0.):
    record = catalog_record()
    path = write_standard_msp([record], tmp_path/'catalog.msp')
    values = [(f.mz,other_intensity if f.fragment_type in {'Common','Precursor Ion'} else 100.)
              for f in record.fragments if f.name not in removed]
    spectrum = ExperimentalSpectrum('test',record.precursor_mz,1.,'+',normalize_peaks(values))
    with_search = LipidMS2Searcher(path,precursor_tolerance_ppm=10,fragment_tolerance_ppm=20)
    try:return with_search.score_spectrum(spectrum,top_n=3)
    finally:with_search.close()


def test_catalog_totals_require_both_core_hg_and_two_b_products(tmp_path):
    rows=run_record(tmp_path)
    assert len(rows)==1
    assert rows[0]['passed_required_gates']
    assert rows[0]['resolution_level']=='species_level'
    assert rows[0]['matched_name']=='Hex4Cer(d34:1)'
    assert rows[0]['lcb_score']==0
    assert rows[0]['hg_score']>0 and rows[0]['b_score']>0
    assert rows[0]['other_score']==0
    assert run_record(tmp_path,removed={'Z0'})==[]
    assert run_record(tmp_path,removed={'B-glycan-H2O'})==[]


def test_catalog_other_threshold_is_strict_and_does_not_admit(tmp_path):
    boundary=run_record(tmp_path,other_intensity=.5)
    above=run_record(tmp_path,other_intensity=.50001)
    assert boundary[0]['other_score']==0
    assert above[0]['other_score']>0
    assert run_record(tmp_path,removed={'Z0'},other_intensity=100)==[]


def test_catalog_extension_keeps_old_blocks_and_requires_eligible_parent(tmp_path):
    import gzip
    from scripts.extend_positive_catalog_glycosphingolipid_library import extend
    from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import iter_blocks
    parent=LibraryRecord(0,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',538.5199,
        '[M+H]+','C34H67NO3','+',[])
    short=LibraryRecord(1,'Cer','Cer(d30:1)','Cer(d18:1/12:0)',482.4573,
        '[M+H]+','C30H59NO3','+',[])
    plain=write_standard_msp([parent,short],tmp_path/'parents.msp')
    source=tmp_path/'parents.msp.gz'
    with gzip.open(source,'wb') as handle:handle.write(plain.read_bytes())
    catalog=tmp_path/'catalog.csv'
    catalog.write_text('# Formula,RT,Mass,Cpd,Comments\nC58H107NO23,,,GalGalGalGlcCer(d34:1),\nC54H99NO23,,,GalGalGalGlcCer(d30:1),\n',encoding='utf-8')
    dest=tmp_path/'extended.msp.gz'
    result=extend(source,dest,catalog)
    assert result['added_records']==1
    assert list(iter_blocks(dest))[:2]==list(iter_blocks(source))
    assert result['preserved_records_and_order']
    from sphingolipid_toolkit._vendor.lipidgate.ms2.library import iter_standard_msp_records
    added=list(iter_standard_msp_records(dest))[-1]
    assert added.polarity=='+' and added.lipid_name=='Hex4Cer(d34:1)'


@pytest.mark.parametrize('cls,head_mz',[('CAEP',126.0315),('N-methylCAEP',140.0471)])
def test_phosphonate_requires_head_and_loss_and_reports_totals(tmp_path,cls,head_mz):
    from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_phosphoceramide import fragments as phospho_fragments
    record=LibraryRecord(1,cls,f'{cls}(d34:1)',f'{cls}(d34:1)',645.5335,'[M+H]+',polarity='+',
        fragments=phospho_fragments(cls,645.5335))
    head=next(f for f in record.fragments if f.name=='HG-phosphonate')
    assert head.mz==pytest.approx(head_mz,abs=.0001)
    path=write_standard_msp([record],tmp_path/'phosphonate.msp')
    searcher=LipidMS2Searcher(path,precursor_tolerance_ppm=10,fragment_tolerance_ppm=20)
    try:
        values=[(f.mz,100.) for f in record.fragments if f.fragment_type=='Diagnostic_HG']
        spectrum=ExperimentalSpectrum('test',record.precursor_mz,1.,'+',normalize_peaks(values))
        rows=searcher.score_spectrum(spectrum)
        assert rows and rows[0]['resolution_level']=='species_level'
        assert rows[0]['lcb_score']==0
        spectrum.peaks=normalize_peaks([(f.mz,100.) for f in record.fragments if f.name!='HG-phosphonate'])
        assert searcher.score_spectrum(spectrum)==[]
    finally:searcher.close()


@pytest.mark.parametrize('cls',['Ac1HexCer','Ac2HexCer','Ac5HexCer','Ac6HexCer','MIP2C','HexSPB','Hex3SPB','SHexSPB','KetoSPB'])
def test_legacy_type_templates_require_their_structural_evidence(tmp_path,cls):
    from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_catalog_legacy_sphingolipid import fragments as factory,required_names
    record=LibraryRecord(1,cls,f'{cls}(d34:1)',f'{cls}(d34:1)',1000.,'[M+H]+',polarity='+',
        fragments=factory(cls,1000.))
    path=write_standard_msp([record],tmp_path/'legacy_type.msp')
    searcher=LipidMS2Searcher(path,precursor_tolerance_ppm=10,fragment_tolerance_ppm=20)
    try:
        values=[(f.mz,100.) for f in record.fragments if f.fragment_type=='Diagnostic_HG']
        spectrum=ExperimentalSpectrum('test',1000.,1.,'+',normalize_peaks(values))
        rows=searcher.score_spectrum(spectrum)
        assert len(rows)==1 and rows[0]['resolution_level']=='species_level'
        assert rows[0]['other_score']==0 and rows[0]['lcb_score']==0
        absent=sorted(required_names(cls))[0]
        spectrum.peaks=normalize_peaks([(f.mz,100.) for f in record.fragments if f.name!=absent])
        assert searcher.score_spectrum(spectrum)==[]
    finally:searcher.close()


def test_reference_hydroxyl_position_suffix_does_not_add_two_oxygens(tmp_path):
    import gzip
    from scripts.extend_positive_catalog_glycosphingolipid_library import extend
    from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import iter_blocks,header_value
    parent=LibraryRecord(0,'Cer','Cer(t34:1)','Cer(d18:1/h16:0)',554.5143,'[M+H]+','C34H67NO4','+',[])
    plain=write_standard_msp([parent],tmp_path/'parent.msp')
    source=tmp_path/'parent.msp.gz'
    with gzip.open(source,'wb') as handle:handle.write(plain.read_bytes())
    catalog=tmp_path/'catalog.csv'
    catalog.write_text('# Formula,RT,Mass,Cpd,Comments\nC36H73N2O6P,,,CAEP(d34:1)(2OH),\n',encoding='utf-8')
    destination=tmp_path/'extended.msp.gz'
    result=extend(source,destination,catalog)
    assert result['added_records']==1
    added=list(iter_blocks(destination))[-1]
    assert header_value(added,'Name')=='CAEP(t34:1)'
    assert header_value(added,'Formula')=='C36H73N2O6P'
