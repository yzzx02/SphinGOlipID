"""Synthetic-only evidence scoring/eligibility audit. Never opens a sample file.

Run from any working directory: python scripts/multi_evidence_shadow.py
Outputs are confined to audit/multi_evidence* and generated documentation tables.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "src"))

import pandas as pd

from sphingolipid_toolkit import ms2_legacy_core as core
from sphingolipid_toolkit.evidence import SOURCE_FUNCTIONS, annotate_theoretical_fragments, capture_fragment_origins
from sphingolipid_toolkit.evidence_config import EvidenceScoringConfig, class_mapping_records
from sphingolipid_toolkit.multi_evidence_adapter import generate_evidence_fragments_for_candidates
from sphingolipid_toolkit.multi_evidence_scoring import PrecursorEvidence, evaluate_candidate_evidence, rank_evidence_candidates
from sphingolipid_toolkit.targeted_mzml_pipeline import MS2Spectrum, match_observed_to_theoretical, finalize_spectrum_matches


def generated(name, function, structure, mass):
    return generate_evidence_fragments_for_candidates(
        pd.DataFrame({"name":[name],"classy":[function],"structure":[structure]}),mass,5.,1.)


def evaluate_scenario(scenario, lipid_class, theory, observed, mass, config=EvidenceScoringConfig()):
    spectrum = MS2Spectrum(Path("synthetic_never_opened.mzML"),scenario,mass,5.,25.,
        observed.fragment_mz.to_numpy(),observed.fragment_intensity.to_numpy())
    # These unchanged production functions implement the 542fa7f legacy score.
    old = finalize_spectrum_matches(match_observed_to_theoretical(spectrum,theory,20e-6,0.),1,0.,100)
    results = []
    for annotation, candidate in theory.groupby("anno",sort=True):
        rows = old[old["注释"].eq(annotation)] if not old.empty else pd.DataFrame()
        old_row = rows.iloc[0] if not rows.empty else None
        result = evaluate_candidate_evidence(lipid_class,annotation,candidate,observed,
            PrecursorEvidence(True,0.,"[M+H]+","synthetic_declared_existing_match"),config=config,
            total_composition=candidate.total.iloc[0],
            legacy_match_score=float(old_row["匹配度分数"]) if old_row is not None else 0.,
            legacy_total_score=float(old_row["总分数"]) if old_row is not None else 0.)
        result.summary.update({"scope":"synthetic_only","scenario":scenario,"precursor":mass,
            "legacy_rank":old["注释"].tolist().index(annotation)+1 if old_row is not None else None})
        results.append(result.summary)
    return rank_evidence_candidates(results)


def scenarios():
    sm = generated("SM(d34:1)","SM","",700.)
    names = ["SM(d18:1/16:0)","SM(d16:1/18:0)"]
    sm = sm[sm.anno.isin(names)].copy()
    peaks = [(184.0733,300.),(700.-183.066,100.),(60.044,40.)]
    # Same precursor and class. Both chains cover two of three LCB fragments;
    # A has more raw intensity, B has more balanced relative LCB support.
    for name,intensities in zip(names,([1000.,1.],[100.,100.])):
        candidate = sm[sm.anno.eq(name)]
        masses = candidate[candidate.fragment_type.eq("LCB")].mz.tolist()[:2]
        peaks.extend(zip(masses,intensities))
    yield "balanced_vs_concentrated_LCB","SM",sm,pd.DataFrame(peaks,columns=["fragment_mz","fragment_intensity"]),700.
    for lipid,function,structure,mass in [
        ("GM3","","NeuAc-Gal-Glc",1200.),("Cer","Cer","",600.),
        ("LacCer","","Gal-Glc",1000.),("So","So","",300.2898)]:
        name = lipid + ("(d18:1)" if lipid == "So" else "(d34:1)")
        theory = generated(name,function,structure,mass)
        annotation = name if lipid == "So" else lipid+"(d18:1/16:0)"
        theory = theory[theory.anno.eq(annotation)].copy()
        selected = theory[theory.scoring_eligible].copy()
        intensities = [100. if kind == "common" else 30. for kind in selected.fragment_type]
        observed = pd.DataFrame({"fragment_mz":selected.mz.tolist(),"fragment_intensity":intensities})
        yield lipid+"_supported",lipid,theory,observed,mass
        if lipid == "GM3":
            absent_hg = observed[~observed.fragment_mz.isin(theory[theory.fragment_type.eq("HG")].mz)]
            yield "GM3_missing_HG","GM3",theory,absent_hg,mass


def function_audit():
    rows = []
    mapping = class_mapping_records()
    for function_name in SOURCE_FUNCTIONS:
        lipid = next((row["lipid_class"] for row in mapping if function_name in row["generator_functions"].split("; ")),"UNSPECIFIED")
        if function_name == "GSL_fragments":
            lipid = "GM3"
        core.frag_id1,core.frag_mz1 = [],[]
        with TemporaryDirectory(prefix="sphingolipid_function_audit_") as temporary:
            old_folder = core.folder
            try:
                core.folder = temporary+os.sep
                with capture_fragment_origins() as origins:
                    if function_name == "GSL_fragments":
                        core.GSL_fragments("NeuAc-Gal-Glc","",1500.,"GM3(d34:1)")
                    elif function_name == "OCer":
                        core.OCer(50,1500.,"1-O-acetyl-Cer(d50:1)",5.,1.)
                    else:
                        getattr(core,function_name)(1500.,lipid+"(d34:1)")
            finally:
                core.folder = old_folder
        if not core.frag_id1:
            rows.append({"scope":"synthetic_only","function":function_name,"lipid_class":lipid,
                "fragment_name":"(no fragments emitted)","classification_status":"UNSPECIFIED", "scoring_eligible":False})
            continue
        annotated = annotate_theoretical_fragments(pd.DataFrame({"id":core.frag_id1,"mz":core.frag_mz1}),lipid,origins)
        for record in annotated.to_dict("records"):
            rows.append({"scope":"synthetic_only","function":function_name,"lipid_class":lipid,**record})
    # shared generators outside class functions: existing LCB and precursor
    # dehydration are inventoried without changing their generation.
    for lipid, labels, masses in [("Cer",["Cer(d34:1)-1H2O","Cer(d34:1)-2H2O"],[600.-18.0106,600.-2*18.0106]),
        ("SM",list(core.LCB_fragment["d18:1"]),core.LCB_mz["d18:1"]),
        ("SM",["[C2H5NO+H]+","SM(d34:1)-1H2O"],[60.044,700.-18.0106])]:
        annotated = annotate_theoretical_fragments(pd.DataFrame({"id":labels,"mz":masses}),lipid)
        for record in annotated.to_dict("records"):
            rows.append({"scope":"synthetic_only","function":"learn_fuc / almost_fuc shared evidence","lipid_class":lipid,**record})
    return pd.DataFrame(rows)


def main():
    audit = ROOT / "audit"
    audit.mkdir(exist_ok=True)
    cases = list(scenarios())
    results = [evaluate_scenario(*case) for case in cases]
    pd.DataFrame([row for result in results for row in result.to_dict("records")]).to_csv(audit / "multi_evidence_scoring_examples.csv",index=False,encoding="utf-8-sig")
    sensitivity = []
    # Illustrative sensitivity grid, not fitted thresholds. Author's 50% group
    # requirement stays fixed; only normalized intensity minimum is varied.
    for minimum in (0.,.01,.05,.10):
        for case in cases:
            sensitivity.append(evaluate_scenario(*case,config=EvidenceScoringConfig(minimum_normalized_intensity=minimum)))
    pd.DataFrame([row for result in sensitivity for row in result.to_dict("records")]).to_csv(audit / "multi_evidence_gate_sensitivity.csv",index=False,encoding="utf-8-sig")
    function_audit().to_csv(audit / "multi_evidence_fragment_classification.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame(class_mapping_records()).to_csv(audit / "multi_evidence_class_mapping.csv",index=False,encoding="utf-8-sig")
    documentation = ROOT / "docs" / "MULTI_EVIDENCE_SCORING.md"
    if documentation.exists():
        content = documentation.read_text(encoding="utf-8")
        columns = ["lipid_class","has_HG","species_gate","molecular_species_gate","scoring_template",
            "HG examples","LCB evidence","NL examples","common examples"]
        table = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"]*len(columns)) + " |"]
        for row in class_mapping_records():
            table.append("| " + " | ".join(str(row[column]).replace("|","\\|") or "—" for column in columns) + " |")
        marker = "<!-- CLASS_MAPPING -->"
        content = content.split(marker)[0] + marker + "\n\n" + "\n".join(table) + "\n"
        documentation.write_text(content,encoding="utf-8")


if __name__ == "__main__":
    main()
