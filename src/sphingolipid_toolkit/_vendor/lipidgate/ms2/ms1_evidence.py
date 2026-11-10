"""Keep measured MS1 precursor evidence independent of feature-table linking."""

from __future__ import annotations

import math

from .config import DEFAULT_SEARCH_CONFIG


CONFIRMED_WITHOUT_FEATURE = "confirmed_ms1_precursor_without_feature"


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return math.nan
    return value


def _present(value):
    return str(value).strip().casefold() not in {"", "nan", "none", "<na>"}


def confirmed_ms1_precursor(row) -> bool:
    """Require the reader's unique, positive centroid and adjacent-survey proof.

    The refiner uses only this file's MS1 scans, without consulting a library.
    Verify its provenance, mass and raw acquisition time when reusing an audit;
    a flag alone, a plotted signal, or an unmatched feature is insufficient.
    """
    if str(row.get("precursor_refinement_status", "")) != "confirmed":
        return False
    if str(row.get("precursor_mz_source", "")) not in {
        "MS1_REFERENCE_CENTROID", "MS1_PRECEDING_CENTROID",
    }:
        return False
    if not (_present(row.get("precursor_ms1_scan_id")) or
            _present(row.get("precursor_ms1_native_id"))):
        return False
    count = _number(row.get("precursor_confirmation_count"))
    mz = _number(row.get("precursor_mz"))
    ms1_mz = _number(row.get("precursor_ms1_mz"))
    ms1_rt = _number(row.get("precursor_ms1_rt_raw_min"))
    ms2_rt = _number(row.get("ms2_rt_raw_min"))
    if not math.isfinite(ms2_rt):
        ms2_rt = _number(row.get("rt_minutes"))
    if not all(math.isfinite(v) for v in (count, mz, ms1_mz, ms1_rt, ms2_rt)):
        return False
    if count < 1 or not count.is_integer() or mz <= 0 or ms1_mz <= 0 or ms1_rt < 0:
        return False
    ppm = abs(mz - ms1_mz) / ms1_mz * 1e6
    gap = (ms2_rt - ms1_rt) * 60
    return (ppm <= DEFAULT_SEARCH_CONFIG.ms1_precursor_confirmation_ppm
            and -1e-7 <= gap <= DEFAULT_SEARCH_CONFIG.ms1_precursor_max_gap_seconds + 1e-7)


def has_ms1_evidence(row) -> bool:
    return (str(row.get("ms1_support_status", "")) == "MS1-supported"
            or confirmed_ms1_precursor(row))


def restore_ms1_support(data):
    """Preserve confirmed physical MS1 support after an unsuccessful feature link.

    Keep feature IDs, apexes, bounds and link failures as they were: confirmation
    is MS1 evidence, not a fabricated detected feature or chromatographic RT.
    """
    out = data.copy()
    if out.empty or "precursor_refinement_status" not in out:
        return out
    if "ms1_support_status" not in out:
        out["ms1_support_status"] = "MS2-only"
    if not (out["ms1_support_status"].ne("MS1-supported") &
            out["precursor_refinement_status"].eq("confirmed")).any():
        return out
    import pandas as pd

    confirmed = pd.Series([confirmed_ms1_precursor(row) for _, row in out.iterrows()],
                          index=out.index, dtype=bool)
    recover = out["ms1_support_status"].ne("MS1-supported") & confirmed
    if recover.any():
        if "ms1_support_reason" in out:
            out.loc[recover, "ms1_feature_link_reason"] = out.loc[recover, "ms1_support_reason"]
        out.loc[recover, "ms1_support_status"] = "MS1-supported"
        out.loc[recover, "ms1_support_reason"] = CONFIRMED_WITHOUT_FEATURE
        out.loc[recover, "ms1_evidence_source"] = "confirmed_ms1_precursor"
    return out
