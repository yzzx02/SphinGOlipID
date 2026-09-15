import pytest
from sphingolipid_toolkit.evidence_config import (
    CLASS_RULES, STRUCTURAL_EVIDENCE_REGISTRY, POOL_WEIGHTS, ClassEvidenceRule,
    EvidenceScoringConfig, class_mapping_records, get_class_rule,
)


def test_fixed_weights_and_registered_two_chain_policies():
    assert dict(POOL_WEIGHTS) == {"primary":60,"secondary":20,"support":20}
    assert CLASS_RULES is STRUCTURAL_EVIDENCE_REGISTRY
    for rule in CLASS_RULES.values():
        if rule.chain_count == 2:
            assert rule.scoring_policy == "structural_60_20_20"


@pytest.mark.parametrize("lipid,primary,secondary", [
    ("SM","HG","LCB"),("GM3","HG","LCB"),("PE_cer","HG","LCB"),
    ("N_CAEP","HG","LCB"),("CerP","HG","LCB"),
    ("Cer","LCB","diagnostic_nl"),("HexCer","diagnostic_nl","LCB"),
    ("Hex2Cer","diagnostic_nl","LCB"),("LacCer","diagnostic_nl","LCB"),
    ("Gb3Cer","diagnostic_nl","LCB"),("type-I-B","diagnostic_nl","LCB")])
def test_subclass_mapping(lipid,primary,secondary):
    rule = get_class_rule(lipid)
    assert (rule.primary_evidence,rule.secondary_evidence) == (primary,secondary)


def test_all_class_rows_have_explicit_gates_or_unspecified_status():
    rows = class_mapping_records()
    assert len(rows) == len(CLASS_RULES)
    for row in rows:
        assert "common" not in row["species_gate"]
        assert "common" not in row["molecular_species_gate"]
        assert row["scoring_policy"] in {"structural_60_20_20","UNSPECIFIED"}
    assert get_class_rule("MysteryCer").status == "UNSPECIFIED"


def test_single_chain_rules_are_explicit_not_cer_defaults():
    assert get_class_rule("So").species_groups == ("LCB",)
    assert get_class_rule("S1P").species_groups == ("HG",)
    assert get_class_rule("S1P").molecular_groups == ("HG","LCB")
    assert get_class_rule("Glu_So").molecular_groups == ("NL","LCB")
    for name in ("So","S1P","Lyso_SM","Lyso_sulfo","Glu_So","Gb3_So"):
        assert get_class_rule(name).chain_count == 1
        assert get_class_rule(name).scoring_policy == "UNSPECIFIED"


def test_registry_can_be_explicitly_extended():
    custom = ClassEvidenceRule("Example",(),False,("NL",),("NL","LCB"),"diagnostic_nl","LCB")
    assert get_class_rule("Example",{"Example":custom}) is custom


def test_custom_registry_classification_and_no_common_gate():
    from sphingolipid_toolkit.evidence import classify_fragment
    custom = ClassEvidenceRule("ReviewedAlias",("SM",),True,("HG",),("HG","LCB"),"HG","LCB")
    record = classify_fragment("ReviewedAlias","[phosphocholine+H]+",registry={"ReviewedAlias":custom})
    assert record.fragment_type.value == "HG"
    with pytest.raises(ValueError):
        ClassEvidenceRule("Invalid",(),False,("common",),("common","LCB"),None,None)
    with pytest.raises(ValueError):
        ClassEvidenceRule("Invalid",(),True,("HG",),("LCB",),None,None)


@pytest.mark.parametrize("kwargs", [{"primary_half_intensity":0},{"secondary_half_intensity":-1},{"primary_half_intensity":float("nan")}, {"precursor_cluster_exclusion_da":-1},
    {"minimum_normalized_intensity":2},{"minimum_group_fraction":0}])
def test_invalid_config_rejected(kwargs):
    with pytest.raises(ValueError):
        EvidenceScoringConfig(**kwargs)


def test_support_and_precursor_cannot_be_structural_pools():
    for evidence in ("Precursor", "common", "supporting_nl"):
        with pytest.raises(ValueError):
            ClassEvidenceRule("bad",(),False,("NL",),("NL","LCB"),evidence,"LCB")


def test_single_chain_scoring_policy_cannot_be_invented():
    with pytest.raises(ValueError):
        ClassEvidenceRule("So",(),False,("LCB",),("LCB",),"LCB","HG",1)
