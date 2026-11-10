"""Lightweight public API; load search/workflow only when requested."""
from importlib import import_module
from .models import ExperimentalPeak, ExperimentalSpectrum, FragmentRecord, LibraryRecord
from .native_policy import install_native_policies

install_native_policies()

_LAZY = {
    "LipidMS2Searcher": ".search",
    "MS2FeatureAnnotationResult": ".workflow",
    "MS2SearchResult": ".workflow",
    "run_ms2_feature_annotation_result": ".workflow",
    "run_ms2_search": ".workflow",
    "run_ms2_search_result": ".workflow",
}
__all__ = ["ExperimentalPeak", "ExperimentalSpectrum", "FragmentRecord", "LibraryRecord", *_LAZY]


def __getattr__(name):
    if name not in _LAZY:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(_LAZY[name], __name__), name)
    globals()[name] = value
    return value
