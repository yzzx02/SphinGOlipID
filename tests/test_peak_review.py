import gzip
import hashlib
import json

import numpy as np
import pandas as pd

from sphingolipid_toolkit.peak_review import _matched_ions, peak_figure


def test_peak_review_fills_the_exact_interval_and_retains_every_match():
    masses = [252.2680, 264.2692, 282.2758, 602.5852, 620.5982, 638.6110, 782.6571, 800.6732]
    intensities = [11, 82, 22, 19, 70, 13, 15, 2]
    row = pd.Series({"lipid_name": "HexCer(d18:1/h22:0)", "adduct": "[M+H]+", "polarity": "+",
                     "observed_mz": 800.65743, "apex_rt_min": 13.35, "apex_intensity": 72,
                     "rt_min": 13.0825, "rt_max": 13.8625, "mz_min": 800.6500, "mz_max": 800.6590,
                     "native_scan_id": "scan=4", "total_score": 96, "matched_fragment_count": 8,
                     "实际mz": str(masses), "强度": str(intensities), "碎片": [f"ion {i}" for i in range(8)]})
    trace = np.asarray([[13.0, 0], [13.1, 8], [13.35, 72], [13.7, 39], [13.8, 11], [14.0, 0]])
    spectrum = np.asarray(masses + [1000]), np.asarray(intensities + [2])
    figure = peak_figure(row, "source.mzML", trace, spectrum, 13.60)
    eic, ms2 = figure.axes

    interval, = eic.collections
    boundary_x, boundary_y = interval.get_paths()[0].vertices.T
    assert (boundary_x.min(), boundary_x.max()) == (row.rt_min, row.rt_max)
    assert boundary_y.min() == 0 and boundary_y.max() == 100
    assert len(eic.patches) == 0
    np.testing.assert_array_equal(eic.lines[0].get_xdata(), trace[:, 0])
    np.testing.assert_allclose(eic.lines[0].get_ydata(), 100 * trace[:, 1] / trace[:, 1].max())
    _, apex = eic.lines
    assert list(apex.get_xdata()) == [row.apex_rt_min]
    assert eic.get_legend() is None and ms2.get_legend() is None
    assert eic.get_xlabel() == "RT (min)"
    acquired, matched = ms2.collections
    assert len(acquired.get_segments()) == len(masses) + 1
    assert [line[0, 0] for line in matched.get_segments()] == masses
    assert 1 < len(ms2.texts) <= 4
    assert {label.xy for label in ms2.texts} <= set(zip(masses, 100 * np.array(intensities) / max(intensities)))
    assert all(label.arrow_patch is None for label in ms2.texts)
    visible_text = {label.get_text() for label in figure.texts}
    assert {"MS1", "MS2", "m/z 800.6500–800.6590  ·  RT 13.082–13.863 min", "RT 13.600 min"} <= visible_text
    assert not any("scan=4" in label or "XCMS integration" in label for label in visible_text)
    assert {f"ion {i}" for i in range(8)} <= visible_text


def test_review_uses_scoring_pools_from_the_verified_workbook_evidence(tmp_path):
    evidence = json.dumps({"fragments": [
        {"name": "LCB-2H2O", "observed_mz": 264.2692, "intensity": 120, "scoring_pool": "lcb"},
        {"name": "LCB-H2O", "observed_mz": 282.2758, "intensity": 15, "scoring_pool": "other",
         "library_pool": "lcb", "auxiliary_lcb_support": True},
        {"name": "unmatched", "observed_mz": None},
    ]})
    key = hashlib.sha256(evidence.encode("utf-8")).hexdigest()
    sidecar = tmp_path / "ms2-evidence-test.json.gz"
    sidecar.write_bytes(gzip.compress(json.dumps({key: evidence}).encode("utf-8")))
    row = {"ms2_evidence_json": json.dumps({"external_evidence_sha256": key}), "evidence_file": sidecar.name,
           "实际mz": [264.2692, 282.2758], "强度": [120, 15], "碎片": ["LCB-2H2O", "LCB-H2O"]}

    assert [ion.pool for ion in _matched_ions(row, tmp_path)] == ["LCB", "Other (aux. LCB)"]
    sidecar.write_bytes(gzip.compress(json.dumps({key: evidence + " "}).encode("utf-8")))
    fallback = _matched_ions(row, tmp_path)
    assert [ion.pool for ion in fallback] == ["", ""]
    assert [ion.mz for ion in fallback] == row["实际mz"]
