"""Content identities for library caches and reproducible run checkpoints."""
from __future__ import annotations

import hashlib
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code_fingerprint(root: Path) -> str:
    """Include uncommitted effective code, not timestamps or unrelated outputs."""
    digest = hashlib.sha256()
    python_paths = list(root.rglob("*.py"))
    if not python_paths and getattr(sys, "frozen", False):
        if root.name == "ms2":
            rules_file = Path(sys._MEIPASS) / "libraries" / "ms2" / "prebuilt" / "cache_rules_sha256.txt"
            if rules_file.is_file():
                return rules_file.read_text(encoding="ascii").strip()
        return sha256(Path(sys.executable))
    paths = sorted(path for path in root.rglob("*")
                   if path.is_file() and path.suffix in {".py", ".cpp", ".h", ".dll", ".pyd", ".so", ".dylib"})
    for path in paths:
        digest.update(path.relative_to(root).as_posix().encode("utf-8") + b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()
