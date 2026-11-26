from sphingolipid_toolkit.engine import SphingolipidIdentifier
from sphingolipid_toolkit.library_catalog import reference_name_candidates
from sphingolipid_toolkit._vendor.lipidgate.ms2.positive_sphingolipid import (
    chain_double_bond_violations, impossible_two_chain_sum_double_bonds,
)


def test_positive_lcb_trienes_use_the_same_cap_at_every_carbon_length():
    for carbon in (16, 18, 20, 22):
        for series in 'mdt':
            assert not chain_double_bond_violations(f'Cer({series}{carbon}:3/18:3)', 'Cer')
            assert chain_double_bond_violations(f'Cer({series}{carbon}:4/18:0)', 'Cer')
    assert chain_double_bond_violations('Cer(d22:3/18:4)', 'Cer')
    assert not impossible_two_chain_sum_double_bonds('SM(d40:6)')
    assert impossible_two_chain_sum_double_bonds('SM(d40:7)')


def test_packaged_lyso_triene_is_searchable_only_with_both_head_losses():
    precursor = 594.33064448
    desulfated, base = (514.37383, 100.), (352.32101, 35.)
    dehydration = (576.32008, 80.)
    def matches(peaks):
        return [r for r in identifier.identify(precursor, peaks) if r['compound_class'] == 'SHexSPB']
    with SphingolipidIdentifier() as identifier:
        result, = matches([desulfated, base, dehydration])
        assert result['matched_name'] == 'SHexSPB(d22:3)'
        assert result['passed_required_gates']
        assert not matches([desulfated, dehydration])
        assert not matches([base, dehydration])
    assert any(r['reference_name'] == 'Lyso-sulfo(d22:3)' for r in reference_name_candidates(precursor))
