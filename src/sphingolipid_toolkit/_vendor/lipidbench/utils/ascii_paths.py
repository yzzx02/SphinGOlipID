"""Temporary ASCII aliases for native readers that cannot open Unicode paths on Windows."""

from __future__ import annotations

import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path


def staging_parent(path: Path) -> Path | None:
    """Prefer an ASCII location on the source volume for cheap hard links."""
    path = Path(path).resolve()
    if str(path.parent).isascii():
        return path.parent
    anchor = Path(path.anchor)
    return anchor if str(anchor).isascii() else None


@contextmanager
def ascii_mzml_path(path):
    path = Path(path).resolve()
    if os.name != "nt" or str(path).isascii():
        yield path
        return
    parent = staging_parent(path)
    try:
        temporary = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_", dir=parent)
    except OSError:
        temporary = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_")
    with temporary as tmp:
        alias = Path(tmp) / "input.mzML"
        try:
            os.link(path, alias)
        except OSError:
            shutil.copy2(path, alias)
        yield alias
