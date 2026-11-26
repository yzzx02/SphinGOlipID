import json
import pytest

from sphingolipid_toolkit.engine import SphingolipidIdentifier
from sphingolipid_toolkit._vendor.lipidgate.ms2.headgroup_core import GM2_CORE, curate_headgroup
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_ganglioside import ADDUCTS, CLASSES, GangliosideComposition


@pytest.mark.parametrize('family', sorted(CLASSES))
@pytest.mark.parametrize('adduct', sorted(ADDUCTS))
def test_ganglioside_hg_cap_preserves_auxiliary_ion_masses(family, adduct):
    model=GangliosideComposition(family,'d',18,1,False,18,0)
    original=model.fragments(adduct)
    curated=curate_headgroup(family,adduct,original)
    assert len([p for p in curated if p.fragment_type=='Diagnostic_HG'])<=8
    assert {(p.name,p.mz) for p in original}=={(p.name.removeprefix('Aux '),p.mz) for p in curated}
    assert all(p.fragment_type=='Common' for p in curated if p.name.startswith('Aux '))
    assert curated==curate_headgroup(family,adduct,curated)


def test_bundled_gm2_ammonium_reports_total_and_credits_other_separately():
    with SphingolipidIdentifier() as engine:
        parent, = [r for r in engine._searcher('+').find_candidates(1401.85765)
                   if r.compound_class=='GM2' and r.adduct=='[M+NH4]+']
        core={f.name for f in parent.fragments if f.fragment_type in ('Diagnostic_HG','Diagnostic_B')}
        assert core==GM2_CORE and len(core)==8
        assert parent.fragments==curate_headgroup('GM2',parent.adduct,parent.fragments)
        byname={f.name:f.mz for f in parent.fragments}
        make=lambda names:engine.identify(parent.precursor_mz,[(byname[n],100.) for n in names],precursor_charge=1)
        required=['B1alpha','B1alpha-H2O','B2Y2alpha','B2']
        rows=make(required+['Aux Y1'])
        row,=[r for r in rows if r['compound_class']=='GM2']
        assert row['resolution_level']=='species_level' and '/' not in row['matched_name']
        proof=json.loads(row['ms2_evidence_json'])
        item,=[f for f in proof['fragments'] if f['name']=='Aux Y1']
        assert item['scoring_pool']=='other' and item['contributes_to_score']
        assert not any(r['compound_class']=='GM2' for r in make(required[:3]+['Aux Y1']))
