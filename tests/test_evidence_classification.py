import pandas as pd
import pytest

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.evidence import (
    EvidenceType, FragmentOrigin, SOURCE_FUNCTIONS, annotate_theoretical_fragments,
    capture_fragment_origins, classify_fragment,
)
from sphingolipid_toolkit.multi_evidence_adapter import generate_evidence_fragments_for_candidates


@pytest.mark.parametrize("lipid,label,kind", [
    ("SM","[phosphocholine+H]+","HG"), ("SM","d18:1","LCB"),
    ("SM","SM(d34:1)-H2O","common"), ("SM","SM(d34:1)-C5H14NO4P","NL"),
    ("GM3","NeuAc+H","HG"), ("GM3","NeuAc+H-H2O","HG"),
    ("GM3","d18:1+OH","LCB"), ("GM3","GM3(d34:1)-NeuAc","NL"),
    ("GM3","GM3(d34:1)-NeuAc-Gal-Glc","NL"), ("GM3","GM3(d34:1)-H2O","common"),
    ("Cer","d18:1","LCB"), ("Cer","M+H-H2O","NL"), ("Cer","M+H-2H2O","NL"),
    ("Cer","Cer(d34:1)-1H2O","NL"), ("Cer","Cer(d34:1)-2H2O","NL"),
    ("Cer","Cer(t34:1)-3H2O","common"),
    ("HexCer","HexCer(d34:1)-Glc","NL"), ("LacCer","LacCer(d34:1)-Gal-Glc","NL"),
    ("GM1","GM1(d34:1)-NeuAc [branch-first 1]","NL"),
    ("GM1","GM1(d34:1)-Gal [main 1; removed branches []; retained branches [0]]-H2O","NL"),
    ("N-CAEP","[C3H10NO3P+H]+]","HG"),
    ("EO-Cer","EO-Cer(d34:1)-E 18:2 -H2O","NL"),
])
def test_existing_labels(lipid,label,kind):
    assert classify_fragment(lipid,label).fragment_type.value == kind


def test_cer_override_is_class_specific_and_author_attributed():
    record = classify_fragment("Cer","M+H-H2O")
    assert record.gate_eligible and record.diagnostic_strength == "diagnostic"
    assert record.classification_source == "author:Cer_dehydration"
    assert classify_fragment("SM","M+H-H2O").fragment_type == EvidenceType.COMMON


def test_unknown_and_mass_only_labels_do_not_invent_evidence():
    for label in ("184.0733", "FA16:0", "d99:99", "unverified diagnostic"):
        record = classify_fragment("SM",label)
        assert record.classification_status == "UNSPECIFIED"
        assert not record.scoring_eligible and not record.gate_eligible
    assert not classify_fragment("Unknown","[phosphocholine+H]+").gate_eligible


def test_source_is_preferred_over_label_fallback():
    record = classify_fragment("Unknown","[phosphocholine+H]+",FragmentOrigin("SM"))
    assert record.fragment_type == EvidenceType.HG
    assert record.classification_source == "generator:SM"


def test_precursor_is_not_a_scoring_fragment():
    record = classify_fragment("SM","[M+H]+",FragmentOrigin("MS1"))
    assert record.fragment_type == EvidenceType.PRECURSOR
    assert not record.scoring_eligible and not record.gate_eligible


def test_all_class_specific_functions_are_observed_without_changing_output(monkeypatch):
    for source in SOURCE_FUNCTIONS:
        function = getattr(core,source)
        assert hasattr(function,"__wrapped__")
        if source == "OCer":
            continue  # file-producing path has a dedicated compatibility test
        args = ("Gal-Gal(-Fuc)-GlcNAc-Gal-Glc","",1500.,"synthetic") if source == "GSL_fragments" else (1500.,"synthetic")
        monkeypatch.setattr(core,"frag_id1",[])
        monkeypatch.setattr(core,"frag_mz1",[])
        function.__wrapped__(*args)
        old = (core.frag_id1.copy(),core.frag_mz1.copy())
        core.frag_id1,core.frag_mz1 = [],[]
        with capture_fragment_origins() as origins:
            function(*args)
        assert (core.frag_id1,core.frag_mz1) == old
        assert len(origins) == len(set(zip(*old)))


def test_actual_gm1_generator_metadata_preserves_labels_and_masses(tmp_path):
    from sphingolipid_toolkit.targeted_mzml_pipeline import generate_legacy_fragments_for_candidates
    candidates = pd.DataFrame({"name":["GM1(d34:1)"],"classy":[""],"structure":["Gal-GalNAc-Gal(-NeuAc)-Glc"]})
    old = generate_legacy_fragments_for_candidates(candidates,1500.,5.,1.,tmp_path)
    enriched = generate_evidence_fragments_for_candidates(candidates,1500.,5.,1.)
    assert set(zip(old.anno,old.mz)) == set(zip(enriched.anno,enriched.mz))
    selected = enriched[enriched.anno.eq("GM1(d18:1/16:0)")]
    hg = selected[selected.fragment_type.eq("HG")]
    assert set(hg.mz) == {292.1,274.09}
    nl = selected[selected.fragment_type.eq("NL")]
    assert len(nl[nl.scoring_eligible]) == 2
    assert len(nl) > 2
    assert nl.classification_source.eq("generator:GSL_fragments").all()
    assert selected[selected.fragment_type.eq("LCB")].classification_source.eq("generator:LCB_fragment").all()
    assert all(set(row.idf.split(" | ")) <= set(enriched.loc[enriched.anno.eq(row.anno) & enriched.mz.eq(row.mz),"fragment_name"].iloc[0].split(" | ")) for row in old.itertuples())


def test_identical_mass_priority_is_independent_of_order():
    data = pd.DataFrame({"idf":["[phosphocholine+H]+","d18:1","SM(d34:1)-C5H14NO4P","SM(d34:1)-H2O"],"mz":[100.]*4})
    for seed in range(4):
        result = annotate_theoretical_fragments(data.sample(frac=1,random_state=seed),"SM")
        assert len(result) == 1
        assert result.fragment_type.iloc[0] == "HG"
        assert len(result.fragment_name.iloc[0].split(" | ")) == 4


def test_unverified_glycan_fallback_is_explanation_until_configured():
    table = pd.DataFrame({"id":["LacCer(d34:1)-Gal-Glc"],"mz":[600.]})
    result = annotate_theoretical_fragments(table,"LacCer")
    assert not result.scoring_eligible.iloc[0]
    explicit = annotate_theoretical_fragments(table,"LacCer",eligibility_overrides={table.id.iloc[0]:True})
    assert explicit.scoring_eligible.iloc[0] and explicit.gate_eligible.iloc[0]


def test_captured_and_uncaptured_three_chain_generation_are_identical(tmp_path,monkeypatch):
    import os
    monkeypatch.setattr(core,"folder",str(tmp_path)+os.sep)
    core.frag_id1,core.frag_mz1 = [],[]
    core.OCer.__wrapped__(50,1500.,"1-O-acetyl-Cer(d50:1)",5.,1.)
    before = (tmp_path/"out_DB.csv").read_bytes()
    (tmp_path/"out_DB.csv").unlink()
    core.frag_id1,core.frag_mz1 = [],[]
    with capture_fragment_origins() as origins:
        core.OCer(50,1500.,"1-O-acetyl-Cer(d50:1)",5.,1.)
    assert (tmp_path/"out_DB.csv").read_bytes() == before
    assert all(origin.function == "OCer" for origin in origins.values())
