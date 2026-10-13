"""Reproducible synthetic-only comparisons against the committed audit baseline.

No raw-data paths or historical result files are read. Run from repository root:
python scripts/scientific_logic_shadow.py
"""
from pathlib import Path
import subprocess
import sys
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from sphingolipid_toolkit import ms2_legacy_core as core, rt_iup as rt
from sphingolipid_toolkit import targeted_mzml_pipeline as target

BASELINE = "4ceeae448bee9851c69ead63d323d7f5a8eb50a9"


def baseline_module(stem):
    name = f"sphingolipid_toolkit._shadow_old_{stem}"
    module = types.ModuleType(name)
    module.__file__ = str(ROOT / "src" / "sphingolipid_toolkit" / (stem + ".py"))
    sys.modules[name] = module
    source = subprocess.check_output(["git", "show", f"{BASELINE}:src/sphingolipid_toolkit/{stem}.py"], cwd=ROOT)
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    return module


def theory(masses, name):
    return pd.DataFrame({"idf":[f"f{i}" for i in range(len(masses))],"mz":masses,
        "anno":name,"total":"synthetic","target":500.,"lenDB":len(set(masses)),
        "RT":5.,"Abund":1000.})


def main():
    output = ROOT / "audit"
    output.mkdir(exist_ok=True)
    old_target = baseline_module("targeted_mzml_pipeline")
    old_rt = baseline_module("rt_iup")
    old_core = baseline_module("ms2_legacy_core")
    rows = []
    scenarios = [
        ("one_observed_two_theories",[100.,200.],[50.,40.],
            pd.concat([theory([100.0001,100.0002],"A"),theory([200.],"B")],ignore_index=True)),
        ("two_observed_two_theories",[100.,100.00005],[50.,80.],theory([100.0001,100.0002],"A")),
        ("identical_labels",[100.],[50.],theory([100.,100.],"A")),
    ]
    for scenario,masses,intensities,theoretical in scenarios:
        spectrum = target.MS2Spectrum(Path("synthetic_not_read.mzML"),"synthetic",500.,5.,25.,
            np.array(masses),np.array(intensities))
        results = []
        for implementation in (old_target,target):
            edges = implementation.match_observed_to_theoretical(spectrum,theoretical,20e-6,20.)
            results.append(implementation.finalize_spectrum_matches(edges,1,0.,10))
        for name in sorted(theoretical.anno.unique()):
            row = dict(scope="synthetic_only",baseline_commit=BASELINE,scenario=scenario,precursor=500.,candidate=name)
            for prefix,result in zip(("old","new"),results):
                candidate = result[result["注释"].eq(name)].iloc[0]
                rank = result["注释"].tolist().index(name)+1
                row.update({f"{prefix}_count":len(candidate["实际mz"]),
                    f"{prefix}_matched_intensity":candidate["强度总和"],f"{prefix}_score":candidate["总分数"],
                    f"{prefix}_rank":rank,f"{prefix}_top1":rank == 1})
            rows.append(row)
    pd.DataFrame(rows).to_csv(output / "fragment_matching_shadow_comparison.csv",index=False)

    rows = []
    for label,legacy,positions,normalized in [
        ("GM1"," -Gal -GalNAc -NeuAc -Gal -Glc","0,3","Gal-GalNAc-Gal(-NeuAc)-Glc"),
        ("Figure4_type_I_B"," -Gal -Fuc -Gal -GlcNAc -Gal -Glc","0,2","Gal-Gal(-Fuc)-GlcNAc-Gal-Glc")]:
        values = []
        for module,sequence,location in [(old_core,legacy,positions),(core,normalized,"")]:
            module.frag_id1,module.frag_mz1 = [],[]
            module.GSL_fragments(sequence,location,1500.,label)
            values.append({round(float(mz),6) for mz in module.frag_mz1})
        for mz in sorted(values[0] | values[1]):
            rows.append(dict(scope="synthetic_only",candidate=label,precursor=1500.,theoretical_mz=mz,
                old_present=mz in values[0],new_present=mz in values[1],
                old_encoding=legacy,old_positions=positions,new_encoding=normalized))
    pd.DataFrame(rows).to_csv(output / "glycan_branch_shadow_comparison.csv",index=False)

    rows = []
    x = np.array([16.,18.,20.,22.,24.])
    for curvature in [0.,.004,.01]:
        data = pd.DataFrame({"细类":"synthetic Cer","曲线不饱和度":1.,"x碳数":x,
            "归一化保留时间":.4*x+curvature*(x-20)**2})
        old,new = old_rt.fit_line(data),rt.fit_line(data)
        rows.append(dict(scope="synthetic_only",series=f"curvature_{curvature}",
            old_model=old["拟合类型"],new_model=new["拟合类型"],old_r2=old["R²"],new_r2=new["R²"],
            old_valid=old["拟合成功"],new_valid=new["拟合成功"]))
    pd.DataFrame(rows).to_csv(output / "rt_linear_first_shadow_comparison.csv",index=False)


if __name__ == "__main__":
    main()
