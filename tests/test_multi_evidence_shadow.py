import importlib.util
from pathlib import Path

from sphingolipid_toolkit.evidence_config import EvidenceScoringConfig


def test_real_generation_synthetic_spectrum_explains_rank_change_and_sensitivity():
    path = Path(__file__).resolve().parents[1] / "scripts" / "multi_evidence_shadow.py"
    spec = importlib.util.spec_from_file_location("synthetic_shadow",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    scenario = next(module.scenarios())
    result = module.evaluate_scenario(*scenario,compare_baseline=True).set_index("annotation")
    assert result.loc["SM(d18:1/16:0)","legacy_rank"] == 1
    assert result.loc["SM(d16:1/18:0)","bd428de_rank"] == 1
    assert result.loc["SM(d18:1/16:0)","multi_evidence_rank"] == 1
    assert result.loc["SM(d16:1/18:0)","multi_evidence_rank"] == 2
    assert result.precursor_pass.all()
    assert result.annotation_level.eq("molecular_species").all()
    sensitive = module.evaluate_scenario(*scenario,config=EvidenceScoringConfig(minimum_normalized_intensity=.01)).set_index("annotation")
    assert sensitive.loc["SM(d18:1/16:0)","annotation_level"] == "species"
    assert sensitive.loc["SM(d16:1/18:0)","annotation_level"] == "molecular_species"
    assert result.multi_evidence_score.to_dict() == sensitive.multi_evidence_score.to_dict()
