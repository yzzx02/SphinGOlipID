from dataclasses import replace

import pytest

from sphingolipid_toolkit._vendor.lipidgate.ms2.indexed_library import write_library_index
from sphingolipid_toolkit._vendor.lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from sphingolipid_toolkit._vendor.lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity


def pair():
    common = [FragmentRecord(mz, name, 'Common') for mz, name in (
        (646.6496,'M+H-CH2O-H2O'),(658.6496,'M+H-2H2O'),(676.6602,'M+H-H2O'))]
    common.append(FragmentRecord(694.6708,'[M+H]+','Precursor Ion'))
    h = LibraryRecord(0,'Cer','Cer(t44:1)','Cer(d18:1/h26:0)',694.6708,'[M+H]+',fragments=[
        FragmentRecord(252.2686,'LCB-CH2O-H2O','Common'),
        FragmentRecord(264.2686,'LCB-2H2O','LCB碎片'),
        FragmentRecord(282.2791,'LCB-H2O','LCB碎片'),*common])
    t = replace(h,record_id=1,lipid_chain_name='Cer(t18:0/26:1)',fragments=[
        FragmentRecord(264.2686,'LCB-3H2O','Common'),
        FragmentRecord(282.2791,'LCB-2H2O','LCB碎片'),
        FragmentRecord(300.2897,'LCB-H2O','LCB碎片'),
        FragmentRecord(318.3003,'LCB','Common'),*common])
    return h,t


@pytest.mark.parametrize('reverse',[False,True])
def test_shared_282_pair_keeps_original_fragment_count_ranking_without_t252(tmp_path,reverse):
    h,t = pair()
    source = tmp_path/'ties.sqlite'
    write_library_index(source,[t,h] if reverse else [h,t],identity())
    spectrum = ExperimentalSpectrum('shared282',694.6708,1.,'+',normalize_peaks([(f.mz,1000.) for f in h.fragments]))
    fast = LipidMS2Searcher(source,native_min_candidates=0)
    slow = LipidMS2Searcher(source,use_native_engine=False,use_score_prefilter=False)
    try:
        rows = fast.score_spectrum(spectrum,top_n=2)
        assert rows == slow.score_spectrum(spectrum,top_n=2)
        scores=[fast._score_sphingo_candidate(spectrum,r) for r in (h,t)]
        assert fast._resolve_trihydroxy_lcb_isomers(scores)==scores
        from sphingolipid_toolkit._vendor.lipidgate.ms2.ranking_policy import build_original_rank_tiers
        tiers=build_original_rank_tiers(scores,scores_tied=fast._scores_tied,class_key=fast._normal_class_key)
        assert [len(tier) for tier in tiers]==[1,1]
        assert {r['matched_name'] for r in rows} == {h.lipid_chain_name,t.lipid_chain_name}
        assert all(r['total_score']==100. and r['resolution_level']=='chain_level' for r in rows)
        assert sorted(r['result_rank'] for r in rows)==[1,2]
        assert sorted(r['matched_fragment_count'] for r in rows)==[6,7]
    finally:
        fast.close();slow.close()
