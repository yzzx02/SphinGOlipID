"""Existing result projection, kept independent of native MS readers for GUI use."""

import pandas as pd
from .confidence import identification_confidence
from .ms1_evidence import restore_ms1_support


MS2_RESULT_EXPORT_COLUMNS = (
    "source_file",
    "scan_id",
    "compound_class",
    "matched_name",
    "total_C",
    "total_DB",
    "adduct",
    "precursor_mz",
    "feature_rt",
    "rt_minutes",
    "ms1_support_status",
    "ppm_error",
    "result_rank",
    "final_score",
    "置信度",
    "注释水平",
    "matched_fragments",
)
MS2_RESULT_NUMBER_FORMATS = (
    ("rt_minutes", "0.000"),
    ("precursor_mz", "0.0000"),
    ("ppm_error", "0.00"),
    ("final_score", "0.00"),
)


def annotation_level_label(
    resolution_level: object,
    matched_name: object = "",
    compound_class: object = "",
) -> str:
    normalized = str(resolution_level or "").strip().lower()
    if normalized in {
        "double_bond_level",
        "tentative_double_bond_level",
        "chain_level",
        "tentative_chain_level",
    }:
        return "链水平"
    return "分子种类水平"


def prepare_ms2_result_export_df(combined: pd.DataFrame) -> pd.DataFrame:
    if combined.empty:
        return pd.DataFrame(columns=MS2_RESULT_EXPORT_COLUMNS)
    export_df = restore_ms1_support(combined)
    export_df["置信度"] = [
        identification_confidence(row) for _, row in export_df.iterrows()
    ]
    # The two display times share one coordinate system. Local experiments
    # supply normalized columns; ordinary runs use their acquisition times.
    for target, candidates in {
        "source_file": ("source_file", "selected_source_file"),
        "scan_id": ("scan_id", "selected_scan_id"),
        "rt_minutes": ("ms2_rt_normalized_min", "rt_minutes", "selected_ms2_rt"),
        "feature_rt": (
            "ms1_feature_rt_normalized_min",
            "ms1_feature_rt_raw_min",
            "feature_rt",
        ),
    }.items():
        for source in candidates:
            if source in export_df:
                export_df[target] = export_df[source]
                break
    if "total_score" in export_df.columns:
        export_df["final_score"] = export_df["total_score"]
    elif "final_score" not in export_df.columns:
        export_df["final_score"] = 0.0
    export_df["final_score"] = (
        pd.to_numeric(export_df["final_score"], errors="coerce")
        .clip(lower=0.0, upper=100.0)
        .round(2)
    )
    if "rt_minutes" in export_df.columns:
        export_df["rt_minutes"] = pd.to_numeric(
            export_df["rt_minutes"], errors="coerce"
        ).round(3)
    if "feature_rt" in export_df.columns:
        export_df["feature_rt"] = pd.to_numeric(
            export_df["feature_rt"], errors="coerce"
        ).round(3)
    if "precursor_mz" in export_df.columns:
        export_df["precursor_mz"] = pd.to_numeric(
            export_df["precursor_mz"], errors="coerce"
        ).round(4)
    if "feature_mz" in export_df.columns:
        export_df["feature_mz"] = pd.to_numeric(
            export_df["feature_mz"], errors="coerce"
        ).round(4)
    if "ppm_error" in export_df.columns:
        export_df["ppm_error"] = pd.to_numeric(
            export_df["ppm_error"], errors="coerce"
        ).round(2)
    if "total_C" in export_df.columns:
        export_df["total_C"] = pd.to_numeric(
            export_df["total_C"], errors="coerce"
        ).astype("Int64")
    if "total_DB" in export_df.columns:
        export_df["total_DB"] = pd.to_numeric(
            export_df["total_DB"], errors="coerce"
        ).astype("Int64")
    if "resolution_level" in export_df.columns:
        matched_names = export_df.get(
            "matched_name", pd.Series("", index=export_df.index)
        )
        compound_classes = export_df.get(
            "compound_class", pd.Series("", index=export_df.index)
        )
        export_df["注释水平"] = [
            annotation_level_label(level, matched_name, compound_class)
            for level, matched_name, compound_class in zip(
                export_df["resolution_level"],
                matched_names,
                compound_classes,
            )
        ]
    elif "注释水平" not in export_df.columns:
        export_df["注释水平"] = "分子种类水平"
    tentative = (
        export_df.get("resolution_level", pd.Series("", index=export_df.index))
        .fillna("")
        .astype(str)
        .str.startswith("tentative_")
    )
    if "passed_required_gates" in export_df:
        tentative |= export_df["passed_required_gates"].eq(False)
    labels = export_df["注释水平"].fillna("").astype(str)
    export_df.loc[tentative & ~labels.str.startswith("暂定"), "注释水平"] = (
        "暂定" + labels
    )
    if "ms1_support_status" not in export_df:
        export_df["ms1_support_status"] = "未评估"
    unsupported = export_df["ms1_support_status"].eq("MS2-only")
    if "feature_rt" in export_df:
        export_df.loc[unsupported, "feature_rt"] = float("nan")
    return export_df.reindex(columns=MS2_RESULT_EXPORT_COLUMNS)
