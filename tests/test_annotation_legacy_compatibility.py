"""No sample data or production pipeline execution; artificial output only."""

import ast
from pathlib import Path

import pandas as pd
import pytest

from sphingolipid_toolkit import annotation_legacy_adapter as adapter
from sphingolipid_toolkit.annotation_evidence import EvidenceSummary
from sphingolipid_toolkit.annotation_resolution import RuleRegistry


@pytest.mark.parametrize("enabled", [False, True])
def test_sidecar_preserves_legacy_object_columns_values_scores_and_order(enabled, monkeypatch):
    # Deliberately unsorted, duplicate rows/index and missing values detect mutation.
    legacy = pd.DataFrame({
        "注释": ["mock-B", "mock-A", "mock-B"],
        "总分数": [12.0, 90.0, 12.0], "RT": [None, 1.0, None],
        "碎片": [["support"], ["headgroup"], ["support"]],
    }, index=[7, 2, 7])
    before = legacy.copy(deep=True)
    if not enabled:
        def forbidden(*args, **kwargs):
            pytest.fail("Disabled framework must not evaluate evidence")
        monkeypatch.setattr(adapter, "resolve_annotation", forbidden)
    returned, sidecar = adapter.annotate_sidecar(
        legacy, EvidenceSummary("mock", "unknown", "mock-scan"), RuleRegistry(),
        adapter.AnnotationResolutionConfig(enabled),
    )
    assert returned is legacy
    pd.testing.assert_frame_equal(legacy, before)
    assert (sidecar is not None) == enabled


def test_default_disabled():
    assert adapter.AnnotationResolutionConfig().annotation_resolution_enabled is False
    marker = object()
    returned, sidecar = adapter.annotate_sidecar(marker, None, None)
    assert returned is marker and sidecar is None


def test_existing_execution_modules_have_no_new_framework_imports():
    source = Path(__file__).resolve().parents[1] / "src" / "sphingolipid_toolkit"
    for name in ("ms2_legacy_core.py", "ms2_pipeline.py", "targeted_mzml_pipeline.py",
                 "scoring.py", "rt_iup.py", "config.py", "__init__.py"):
        tree = ast.parse((source / name).read_text(encoding="utf-8-sig"))
        imports = [node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        imports += [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names]
        assert not any("annotation_evidence" in item or "annotation_resolution" in item
                       or "annotation_legacy_adapter" in item for item in imports)
