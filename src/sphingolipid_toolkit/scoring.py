"""Scoring and normalized result-column helpers."""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def ppm_error(observed_mz: float, theoretical_mz: float) -> float:
    """Return signed ppm error for an observed/theoretical m/z pair."""

    return (float(observed_mz) - float(theoretical_mz)) / float(theoretical_mz) * 1_000_000


def match_fragments(
    observed_fragments: pd.DataFrame,
    theoretical_fragments: pd.DataFrame,
    ppm_tolerance: float,
    min_intensity: float = 0,
    observed_mz_col: str = "fragment_mz",
    observed_intensity_col: str = "fragment_intensity",
    theoretical_mz_col: str = "theoretical_mz",
) -> pd.DataFrame:
    """Match normalized observed and theoretical fragment tables.

    This utility is for API-level tests and future frontends. The production
    MS2 pipeline still delegates to the converted legacy matcher.
    """

    required_observed = {observed_mz_col, observed_intensity_col}
    required_theoretical = {theoretical_mz_col}
    missing_observed = required_observed - set(observed_fragments.columns)
    missing_theoretical = required_theoretical - set(theoretical_fragments.columns)
    if missing_observed:
        raise KeyError(f"Observed fragments missing required column(s): {sorted(missing_observed)}")
    if missing_theoretical:
        raise KeyError(f"Theoretical fragments missing required column(s): {sorted(missing_theoretical)}")

    observed = observed_fragments.copy()
    observed[observed_mz_col] = pd.to_numeric(observed[observed_mz_col], errors="coerce")
    observed[observed_intensity_col] = pd.to_numeric(observed[observed_intensity_col], errors="coerce")
    observed = observed.dropna(subset=[observed_mz_col, observed_intensity_col])
    observed = observed[observed[observed_intensity_col] >= min_intensity]

    theoretical = theoretical_fragments.copy()
    theoretical[theoretical_mz_col] = pd.to_numeric(theoretical[theoretical_mz_col], errors="coerce")
    theoretical = theoretical.dropna(subset=[theoretical_mz_col])

    rows: list[dict[str, object]] = []
    for _, obs in observed.iterrows():
        for _, theo in theoretical.iterrows():
            error = ppm_error(obs[observed_mz_col], theo[theoretical_mz_col])
            if abs(error) <= ppm_tolerance:
                row = {
                    "observed_mz": float(obs[observed_mz_col]),
                    "observed_intensity": float(obs[observed_intensity_col]),
                    "theoretical_mz": float(theo[theoretical_mz_col]),
                    "ppm_error": float(error),
                    "matched": True,
                }
                for col in ("fragment_name", "fragment_type", "fragment_formula", "evidence_level", "description"):
                    if col in theoretical.columns:
                        row[col] = theo[col]
                rows.append(row)
    return pd.DataFrame(
        rows,
        columns=[
            "observed_mz",
            "observed_intensity",
            "theoretical_mz",
            "fragment_name",
            "fragment_type",
            "fragment_formula",
            "ppm_error",
            "matched",
            "evidence_level",
            "description",
        ],
    ).dropna(axis=1, how="all")


def score_fragment_matches(matches: pd.DataFrame, total_fragment_intensity: float | None = None) -> dict[str, float | int]:
    """Summarize matched-fragment evidence without changing legacy scoring."""

    if matches.empty:
        matched_intensity_sum = 0.0
        count = 0
    else:
        count = int(matches.get("matched", True).sum()) if "matched" in matches else int(len(matches))
        intensity = pd.to_numeric(matches.get("observed_intensity", pd.Series(dtype=float)), errors="coerce").fillna(0)
        matched_intensity_sum = float(intensity.sum())

    total = float(total_fragment_intensity) if total_fragment_intensity is not None else matched_intensity_sum
    ratio = matched_intensity_sum / total if total else 0.0
    return {
        "matched_fragment_count": count,
        "matched_intensity_sum": matched_intensity_sum,
        "total_fragment_intensity": total,
        "matched_intensity_ratio": ratio,
    }


def add_standard_score_columns(result: pd.DataFrame) -> pd.DataFrame:
    """Add stable English aliases to a legacy MS2 result table."""

    out = result.copy()
    if out.empty:
        for col in _standard_result_columns():
            if col not in out.columns:
                out[col] = pd.Series(dtype="object")
        return out

    if "lipid_name" not in out.columns and "注释" in out.columns:
        out["lipid_name"] = out["注释"]
    if "observed_mz" not in out.columns and "target" in out.columns:
        out["observed_mz"] = out["target"]
    if "observed_rt" not in out.columns and "RT" in out.columns:
        out["observed_rt"] = out["RT"]
    if "signal_intensity" not in out.columns and "Abund" in out.columns:
        out["signal_intensity"] = out["Abund"]
    if "matched_fragment_count" not in out.columns and "实际mz" in out.columns:
        out["matched_fragment_count"] = out["实际mz"].apply(_sequence_length)
    if "matched_intensity_sum" not in out.columns and "强度总和" in out.columns:
        out["matched_intensity_sum"] = out["强度总和"]
    if "matched_intensity_ratio" not in out.columns and "相对强度" in out.columns:
        out["matched_intensity_ratio"] = out["相对强度"]
    if "ms2_match_score" not in out.columns and "匹配度分数" in out.columns:
        out["ms2_match_score"] = out["匹配度分数"]
    if "rank" not in out.columns and {"target", "总分数"}.issubset(out.columns):
        out["rank"] = out.groupby("target")["总分数"].rank(method="dense", ascending=False).astype(int)
    return out


def _standard_result_columns() -> Sequence[str]:
    return (
        "lipid_name",
        "observed_mz",
        "observed_rt",
        "signal_intensity",
        "matched_fragment_count",
        "matched_intensity_sum",
        "matched_intensity_ratio",
        "ms2_match_score",
        "rank",
    )


def _sequence_length(value: object) -> int:
    if isinstance(value, (list, tuple, set)):
        return len(value)
    if pd.isna(value):
        return 0
    return 1

