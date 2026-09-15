import pytest
from sphingolipid_toolkit.evidence_config import (
    CLASS_RULES, SCORING_TEMPLATES, ClassEvidenceRule, EvidenceScoringConfig,
    ScoringTemplate, class_mapping_records, get_class_rule,
)


def test_exact_author_weights_and_sum():
    assert dict(SCORING_TEMPLATES["HG-dominant"].weights) == {"HG":60,"LCB":20,"NL":15,"common":5}
    assert dict(SCORING_TEMPLATES["NL-dominant glycosphingolipid"].weights) == {"LCB":35,"NL":55,"common":10}
    assert dict(SCORING_TEMPLATES["LCB-dominant"].weights) == {"LCB":60,"NL":30,"common":10}
    assert all(sum(dict(template.weights).values()) == 100 for template in SCORING_TEMPLATES.values())


@pytest.mark.parametrize("lipid,template", [("SM","HG-dominant"),("GM3","HG-dominant"),
    ("PE_cer","HG-dominant"),("N_CAEP","HG-dominant"),("CerP","HG-dominant"),
    ("Cer","LCB-dominant"),("HexCer","NL-dominant glycosphingolipid"),
    ("LacCer","NL-dominant glycosphingolipid"),("Gb3Cer","NL-dominant glycosphingolipid")])
def test_subclass_template_mapping(lipid,template):
    assert get_class_rule(lipid).scoring_template == template


def test_all_class_rows_have_explicit_gates_or_unspecified_status():
    rows = class_mapping_records()
    assert len(rows) == len(CLASS_RULES)
    for row in rows:
        assert "common" not in row["species_gate"]
        assert "common" not in row["molecular_species_gate"]
        assert row["scoring_template"] in {*SCORING_TEMPLATES,"UNSPECIFIED"}
    assert get_class_rule("MysteryCer").status == "UNSPECIFIED"


def test_single_chain_rules_are_explicit_not_cer_defaults():
    assert get_class_rule("So").species_groups == ("LCB",)
    assert get_class_rule("S1P").species_groups == ("HG",)
    assert get_class_rule("S1P").molecular_groups == ("HG","LCB")
    assert get_class_rule("Glu_So").molecular_groups == ("NL","LCB")
    for name in ("So","S1P","Lyso_SM","Lyso_sulfo","Glu_So","Gb3_So"):
        assert get_class_rule(name).chain_count == 1
        assert get_class_rule(name).scoring_template is None


def test_registry_can_be_explicitly_extended():
    custom = ClassEvidenceRule("Example",(),False,("NL",),("NL","LCB"),"LCB-dominant")
    assert get_class_rule("Example",{"Example":custom}) is custom


def test_custom_registry_classification_and_no_common_gate():
    from sphingolipid_toolkit.evidence import classify_fragment
    custom = ClassEvidenceRule("ReviewedAlias",("SM",),True,("HG",),("HG","LCB"),"HG-dominant")
    record = classify_fragment("ReviewedAlias","[phosphocholine+H]+",registry={"ReviewedAlias":custom})
    assert record.fragment_type.value == "HG"
    with pytest.raises(ValueError):
        ClassEvidenceRule("Invalid",(),False,("common",),("common","LCB"),None)
    with pytest.raises(ValueError):
        ClassEvidenceRule("Invalid",(),True,("HG",),("LCB",),None)


@pytest.mark.parametrize("kwargs", [{"k_hg":0},{"k_nl":-1},{"k_lcb":float("nan")},
    {"minimum_normalized_intensity":2},{"minimum_group_fraction":0}])
def test_invalid_config_rejected(kwargs):
    with pytest.raises(ValueError):
        EvidenceScoringConfig(**kwargs)


def test_precursor_weights_and_wrong_totals_are_rejected():
    for weights in [(("Precursor",100),),(("HG",90),),(("HG",50),("HG",50))]:
        with pytest.raises(ValueError):
            ScoringTemplate("bad",weights)
