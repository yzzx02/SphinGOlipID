import json
import re
from pathlib import Path
from dataclasses import replace

import pytest
from sphingolipid_toolkit.engine import IdentificationConfig, SphingolipidIdentifier
from sphingolipid_toolkit.acidic_scoring import normalize_acidic_weights
from sphingolipid_toolkit._vendor.lipidgate.ms2 import indexed_library
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION, normalize_imported_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint

BANK = Path(__file__).resolve().parents[1] / 'src/sphingolipid_toolkit/data/ms2/sphingolipid_positive.sqlite'

def specimen(cls, adduct, name, mz):
    library = IndexedLibrary(BANK, allowed_class_keys=[re.sub(r'[^A-Z0-9]','',cls.upper())],
                             allowed_adducts=[adduct], mz_min=mz-.002, mz_max=mz+.002)
    try:
        return normalize_imported_record(next(r for r in library if r.lipid_chain_name == name))
    finally:
        library.close()

def identify_one(tmp_path, record, weights=None):
    path = tmp_path / 'positive.sqlite'
    write_library_index(path, [replace(record,record_id=0)], {
        'version': LIBRARY_CACHE_VERSION, 'rules_sha256': code_fingerprint(Path(indexed_library.__file__).parent)})
    return SphingolipidIdentifier(IdentificationConfig(library_path=path, min_total_score=0, acidic_weights=weights or {}))

def test_lysosphingomyelin_channels_and_distinct_peaks(tmp_path):
    r = specimen('LSM','[M+H]+','LSM(d18:1)',465.3452)
    mz = {f.name:f.mz for f in r.fragments}
    expected = {'[C5H15NO4P]+':'hg','[C2H6O4P]+':'hg','LCB-H2O':'lcb','LCB-2H2O':'lcb',
                'M+H-H2O':'other','M+H-trimethylamine(-59)':'other','[C5H12N]+':'other','[C5H14NO]+':'other'}
    assert set(mz)==set(expected)
    assert [mz[k] for k in ('[C2H6O4P]+','[C5H14NO]+','[C5H12N]+','M+H-H2O')]==pytest.approx([124.9998,104.1070,86.0964,447.3346],abs=.00015)
    with identify_one(tmp_path,r) as app:
        assert app.identify(r.precursor_mz, [(mz['[C2H6O4P]+'],100.),(mz['LCB-2H2O'],30.)])==[]
        gate=[(mz['[C5H15NO4P]+'],100.),(mz['LCB-2H2O'],1.)]
        row,=app.identify(r.precursor_mz,gate+[(mz['[C5H14NO]+'],.5)])
        assert row['fah_score']==row['other_score']==0 and row['lcb_score']>0
        row,=app.identify(r.precursor_mz,[(m,100.) for m in mz.values()])
        ev=json.loads(row['ms2_evidence_json'])
        assert (row['hg_score'],row['lcb_score'],row['other_score'],row['fah_score'],row['final_score'])==(60.,20.,20.,0.,100.)
        assert {x['name']:x['scoring_pool'] for x in ev['fragments']}==expected
        assert len({f['observed_mz'] for f in ev['fragments'] if f['contributes_to_score']})==5

@pytest.mark.parametrize(('cls','add','name','prec'),[
    ('Hex-IPC','[M+H]+','Hex-IPC(d14:0/14:1)',840.4869),
    ('GD1a','[M+2NH4]2+','GD1a(d18:1/16:0)',922.5042), ('GD1b','[M+2NH4]2+','GD1b(d18:1/16:0)',922.5042),
    ('GT1b','[M+2NH4]2+','GT1b(d18:1/16:0)',1068.0519), ('GQ1b','[M+2NH4]2+','GQ1b(d18:1/16:0)',1213.5996),
    ('Gb4','[M+H]+','Gb4(d18:1/16:0)',1227.7572),
    ('Cer1P','[M+H]+','Cer1P(d18:1/16:0)',618.4857)])
def test_positive_chain_products_do_not_score_in_fah(tmp_path,cls,add,name,prec):
    r=specimen(cls,add,name,prec)
    with identify_one(tmp_path,r,{cls:(55,30,15)}) as app:
        row,=app.identify(r.precursor_mz,[(f.mz,100.) for f in r.fragments])
        ev=json.loads(row['ms2_evidence_json'])
        assert ev['score_profile']['weights_pct']==pytest.approx({'hg':55,'lcb':30,'other':15})
        assert (row['hg_score'],row['lcb_score'],row['fah_score'])==pytest.approx((55,30,0))
        for f in ev['fragments']:
            if f['name'].startswith(('LCB-','Ceramide fragment')):
                assert (f['scoring_pool'],f['role'])==('lcb','chain')
            if f['name'].startswith('Reference-HG-'):
                assert f['scoring_pool']=='hg'

@pytest.mark.parametrize(('adduct','prec'),[('[M+H]+',460.3091),('[M+NH4]+',477.3357)])
def test_sulfonolipid_acyl_loss_and_headgroup_weights(tmp_path,adduct,prec):
    r=specimen('SL',adduct,'SL(m10:1/14:1)',prec)
    with identify_one(tmp_path,r,{'SL':(50,35,15)}) as app:
        row,=app.identify(r.precursor_mz,[(f.mz,100.) for f in r.fragments])
        ev=json.loads(row['ms2_evidence_json'])
        assert ev['score_profile']['weights_pct']=={'fah':50,'hg':35,'other':15}
        assert (row['fah_score'],row['hg_score'],row['lcb_score'])==(50,35,0)

def test_documented_caep_chain_requires_an_independent_lcb_product(tmp_path):
    r=specimen('CAEP','[M+H]+','CAEP(d18:3/16:0)',641.5017)
    mz={f.name:f.mz for f in r.fragments}
    head=[(mz['HG-phosphonate'],100.),(mz['M+H-head'],70.)]
    with identify_one(tmp_path,r) as app:
        assert app.identify(r.precursor_mz,head)==[]
        assert app.identify(r.precursor_mz,head+[(mz['LCB-2H2O'],.99)])==[]
        assert app.identify(r.precursor_mz,[(mz['LCB-2H2O'],100.)])==[]
        row,=app.identify(r.precursor_mz,head+[(mz['LCB-2H2O'],1.)])
        assert row['resolution_level']=='chain_level' and row['lcb_score']>0
        assert json.loads(row['ms2_evidence_json'])['score_profile']['weights_pct']=={'hg':60,'lcb':20,'other':20}

def test_old_two_pool_overrides_preserve_explicit_weights_and_adduct_scope():
    assert normalize_acidic_weights({'SM':(80,20)})=={'SM':(80,0,20),'SM [M+Na]+':(80,20)}
    normalized = normalize_acidic_weights({'SM':(80,20),'SM [M+Na]+':(70,30)})
    assert normalized['SM']==(80,0,20)
    assert normalized['SM [M+Na]+']==pytest.approx((70,30))

