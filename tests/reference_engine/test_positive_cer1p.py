from dataclasses import replace
import gzip
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.esterified_ceramide import C,H,N,O,PROTON,WATER
from sphingolipid_toolkit._vendor.lipidgate.ms2.evidence_export import evidence_json
from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.library import LIBRARY_CACHE_VERSION,load_standard_msp
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum,FragmentRecord,LibraryRecord,normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.msp_tools import header_value,iter_blocks
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_cer1p import Cer1PSodiumComposition,GATE_NAMES,PHOSPHORIC_ACID,SODIUM_PHOSPHATE,normalize_h_fragments
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_hexosylceramide import has_chain_evidence
from sphingolipid_toolkit._vendor.lipidgate.ms2.provenance import code_fingerprint
from sphingolipid_toolkit._vendor.lipidgate.ms2.rules import DEFAULT_RULES
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _pool_for_fragment,_pool_weights_for_record
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from sphingolipid_toolkit._vendor.lipidgate.ms2.search_plan import evidence_program
from scripts.curate_positive_cer1p_library import curate,sodium_block,sodium_composition


def h_record(series='d',carbon=18,db=1):
    oh={'m':1,'d':2,'t':3}[series];p=620.5014
    intact=carbon*C+(2*carbon+3-2*db)*H+N+oh*O+PROTON
    fragments=normalize_h_fragments(f'Cer1P({series}{carbon}:{db}/16:0)',p,[
        FragmentRecord(p-PHOSPHORIC_ACID,GATE_NAMES['[M+H]+'],'Diagnostic_HG'),
        FragmentRecord(p-WATER,'M+H-H2O','Common'),
        FragmentRecord(round(intact-2*WATER,4),'LCB-2H2O','LCB碎片')])
    return LibraryRecord(0,'Cer1P',f'Cer1P({series}{carbon+16}:{db})',f'Cer1P({series}{carbon}:{db}/16:0)',
                         p,'[M+H]+',polarity='+',fragments=fragments)


def na_record(oh=2):
    c=Cer1PSodiumComposition(34,0,oh)
    return LibraryRecord(0,'Cer1P',c.name,c.name,c.precursor_mz,'[M+Na]+',formula=c.formula,polarity='+',fragments=c.fragments())


def spec(r,peaks):return ExperimentalSpectrum('Cer1P',r.precursor_mz,1.,'+',normalize_peaks(peaks))


def score(r,peaks):
    search=object.__new__(LipidMS2Searcher);search.fragment_tolerance_da=.0001
    fast=search._score_sphingo_candidate(spec(r,peaks),r)
    with patch('sphingolipid_toolkit._vendor.lipidgate.ms2.native_bridge.load_bridge',return_value=None):
        slow=search._score_sphingo_candidate(spec(r,peaks),r)
    assert fast==slow
    return fast


@pytest.mark.parametrize('series',['m','d','t'])
def test_protonated_cer1p_only_phosphate_admits_other_has_no_gate(series):
    r=h_record(series);gate=next(f for f in r.fragments if f.name==GATE_NAMES[r.adduct])
    assert normalize_h_fragments(r.lipid_chain_name,r.precursor_mz,r.fragments)==r.fragments
    assert _pool_weights_for_record(r,DEFAULT_RULES.get('Cer1P'))=={'fah':0.,'hg':60.,'other':40.}
    assert all(_pool_for_fragment(r,f)==('hg' if f==gate else 'other') for f in r.fragments)
    mask=1<<r.fragments.index(gate)
    assert evidence_program(r)==(((mask,1),),)
    passed=score(r,[(gate.mz,1000.)])
    assert passed.passed_required_gates and passed.total_score==60. and passed.resolution_level=='species_level'
    assert passed.pool_scores['fah'].pool_score==0. and passed.pool_scores['other'].pool_score==0.
    rejected=score(r,[(f.mz,1000.) for f in r.fragments if f!=gate])
    assert not rejected.passed_required_gates and rejected.resolution_level=='class_level'
    full=score(r,[(f.mz,1000.) for f in r.fragments])
    assert full.total_score==100. and full.resolution_level==('species_level' if series=='m' else 'chain_level')
    audit=json.loads(evidence_json(spec(r,[(gate.mz,1000.)]),passed))
    assert sum(f['admission_gate'] for f in audit['fragments'])==1
    assert all(f['role']=='other' and f['scoring_pool']=='other' for f in audit['fragments'] if f['name'].startswith('LCB'))


@pytest.mark.parametrize('series,carbon,db',[('d',18,1),('d',20,2),('t',18,0),('t',19,1)])
@pytest.mark.parametrize('name',['LCB-H2O','LCB-2H2O'])
def test_cer1p_shared_formula_channels_raw_one_percent_is_only_a_reporting_check(series,carbon,db,name):
    r=h_record(series,carbon,db);f={x.name:x for x in r.fragments};gate=f[GATE_NAMES[r.adduct]]
    for intensity,expected in [(0.,'species_level'),(4.999,'species_level'),(5.,'species_level'),
                               (5.001,'species_level'),(9.999,'species_level'),(10.,'chain_level'),(10.001,'chain_level')]:
        result=score(r,[(gate.mz,1000.),(f[name].mz,intensity)])
        assert result.passed_required_gates and result.resolution_level==expected
        assert has_chain_evidence(r,result.matched_fragments)==(expected=='chain_level')
        assert result.pool_scores['other'].matched_count==int(intensity>5.)
        assert result.pool_scores['other'].pool_score==pytest.approx(40/3 if intensity>5. else 0.)


@pytest.mark.parametrize('oh',[1,2,3])
def test_sodium_formula_phosphate_losses_and_sum_only(tmp_path,oh):
    c=Cer1PSodiumComposition(34,0,oh);r=na_record(oh)
    assert c.formula==f'C34H70NO{oh+4}P'
    if oh==2:assert c.precursor_mz==pytest.approx(642.48329669,abs=.0000001)
    f={f.name:f for f in r.fragments}
    assert c.precursor_mz-f['M+Na-H3PO4'].mz==pytest.approx(97.97689557)
    assert c.precursor_mz-f['M+Na-H2NaPO4'].mz==pytest.approx(119.95883982)
    assert len(r.fragments)==5 and all('LCB' not in f.name for f in r.fragments)
    path=tmp_path/'na.msp';path.write_text('\n'.join(sodium_block(c))+'\n',encoding='utf-8')
    rows=load_standard_msp(path)
    assert len(rows)==1 and rows[0].lipid_name==rows[0].lipid_chain_name==c.name
    result=score(rows[0],[(f['M+Na-H3PO4'].mz,1000.)])
    assert result.passed_required_gates and result.resolution_level=='species_level' and result.total_score==60.
    rejected=score(rows[0],[(f.mz,1000.) for f in rows[0].fragments if f.name!='M+Na-H3PO4'])
    assert not rejected.passed_required_gates
    with_lcb=replace(rows[0],fragments=[*rows[0].fragments,FragmentRecord(264.2686,'LCB-2H2O','Common')])
    result=score(with_lcb,[(f.mz,1000.) for f in with_lcb.fragments])
    assert result.resolution_level=='species_level' and result.total_score<=100.


@pytest.mark.parametrize('strength,count',[(4.999,0),(5.,0),(5.001,1)])
def test_cer1p_sodium_other_raw_cutoff(strength,count):
    r=na_record();f={x.name:x for x in r.fragments}
    result=score(r,[(f['M+Na-H3PO4'].mz,1000.),(f['M+Na-H2O'].mz,strength)])
    assert result.passed_required_gates and result.pool_scores['other'].matched_count==count
    assert result.resolution_level=='species_level'


def test_indexed_native_cer1p_scoring_gate_and_resolution_match_reference(tmp_path):
    rows=[h_record('d'),h_record('t',18,0),h_record('m'),na_record()]
    for i,r in enumerate(rows):r.record_id=i
    path=tmp_path/'cer1p.sqlite'
    write_library_index(path,rows,dict(version=LIBRARY_CACHE_VERSION,
        rules_sha256=code_fingerprint(Path(__file__).resolve().parents[2]/'src/sphingolipid_toolkit/_vendor/lipidgate/ms2')))
    fast=LipidMS2Searcher(path,native_min_candidates=0)
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False)
    try:
        for r in rows:
            gate=next(f for f in r.fragments if f.name==GATE_NAMES[r.adduct])
            for strength in [0.,5.,10.]:
                peaks=[(gate.mz,1000.)]+[(f.mz,strength) for f in r.fragments if f.name.startswith('LCB')]
                actual=fast.score_spectrum(spec(r,peaks))
                assert actual==slow.score_spectrum(spec(r,peaks))
                assert actual and all(x['passed_required_gates'] for x in actual)
                assert any(x['resolution_level']=='chain_level' for x in actual)==(strength>=10. and r.lipid_chain_name.startswith(('Cer1P(d','Cer1P(t')) and r.adduct=='[M+H]+')
            peaks=[(f.mz,1000.) for f in r.fragments if f!=gate]
            # Other records in this fixture may legitimately have another loss
            # at the omitted position; check the scoped candidate itself.
            rejected=score(r,peaks)
            assert not rejected.passed_required_gates
    finally:fast.close();slow.close()


def seed(name,adduct='[M+H]+'):
    return [f'Name: {name}','PrecursorMZ: 620.5014',f'PrecursorType: {adduct}','CompoundClass: Cer1P',
            f'Comment: MS1_name={name};polarity=+','Num Peaks: 3',
            '266.2842 100 "LCB-2H2O" "LCB碎片"',
            '522.5245 100 "M+H-H3PO4" "Diagnostic_HG"','602.4908 100 "M+H-H2O" "Common"']


def test_curation_changes_only_protonated_cer1p_appends_deduplicated_na_and_is_idempotent(tmp_path):
    blocks=[seed('Cer1P(d18:0/16:0)'),seed('Cer1P(d16:0/18:0)'),seed('Cer1P(m18:1/16:0)'),seed('Cer1P(t18:0/16:0)')]
    negative=seed('Cer1P(d18:0/16:0)','[M-H]-');blocks.insert(1,negative)
    assert sodium_composition(seed('Cer1P(d8:0/14:0)')) is not None
    assert sodium_composition(seed('Cer1P(d18:0/13:0)')) is None
    assert sodium_composition(seed('Cer1P(d34:0)')) is None
    source,out,repeat=[tmp_path/n for n in ['in.gz','out.gz','repeat.gz']]
    with gzip.open(source,'wt',encoding='utf-8') as f:f.write('\n\n'.join('\n'.join(b) for b in blocks)+'\n')
    # Use valid LCB masses for the two d chain isomers, so preservation checks
    # distinguish reclassification from an unrelated formula correction.
    blocks[2]=[x.replace('266.2842','238.2529') for x in blocks[2]]
    blocks[3]=[x.replace('266.2842','248.2736') for x in blocks[3]]
    blocks[4]=[x.replace('266.2842','282.2791') for x in blocks[4]]
    with gzip.open(source,'wt',encoding='utf-8') as f:f.write('\n\n'.join('\n'.join(b) for b in blocks)+'\n')
    audit=curate(source,out);revised=list(iter_blocks(out))
    assert audit['changed_protonated_by_series']=={'d':2,'m':1,'t':1}
    assert audit['added_sodium_records']==3 and revised[1]==negative
    assert all(header_value(a,'Name')==header_value(b,'Name') for a,b in zip(blocks,revised))
    assert all('/' not in header_value(b,'Name') for b in revised[len(blocks):])
    repeated=curate(out,repeat)
    assert repeated['changed_protonated_by_series']=={} and repeated['added_sodium_records']==0
    assert out.read_bytes()==repeat.read_bytes()
    with pytest.raises(ValueError,match='separate destination'):curate(source,source)
