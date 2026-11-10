"""Load verified C++ compilations of unchanged policy source, or use Python.

The extensions keep Python object arithmetic (including NumPy float32) rather
than inferring C numeric types. Canonical module names preserve types, mutable
module globals, public APIs and audit paths.
"""
from __future__ import annotations

import hashlib
import importlib.abc
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys


POLICY_MODULES = (
    "scoring", "gate_policy", "resolution_policy", "oxidized_fragment_policy",
    "chain_utils", "matching", "record_facts", "indexed_library", "ranking_policy",
    "scoring_policy", "evidence_policy", "search", "search_plan", "native_plan",
    "sphingolipid_scoring",
    "equivalent_candidates", "native_engine", "prepared_scoring",
)


class _PolicyLoader(importlib.machinery.ExtensionFileLoader):
    def exec_module(self, module):
        super().exec_module(module)
        # Existing resource/provenance code uses the logical source directory.
        # __spec__.origin retains the actual verified extension location.
        module.__file__ = str(Path(__file__).parent / (module.__name__.rsplit(".", 1)[1] + ".py"))


class _PolicyFinder(importlib.abc.MetaPathFinder):
    def __init__(self):
        self.root = Path(__file__).parent
        self.directory = self.root / "native_backend" / "policy"
        try:
            self.manifest = json.loads((self.directory / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            self.manifest = {}

    def find_spec(self, fullname, path=None, target=None):
        if os.environ.get("LIPIDGATE_PYTHON_POLICY") == "1" or not fullname.startswith("sphingolipid_toolkit._vendor.lipidgate.ms2."):
            return None
        name = fullname.removeprefix("sphingolipid_toolkit._vendor.lipidgate.ms2.")
        manifest = self.manifest
        directives = manifest.get("directives", {})
        if (name not in POLICY_MODULES or manifest.get("python_abi") != list(sys.version_info[:2])
                or manifest.get("platform") != sys.platform or manifest.get("machine") != platform.machine()
                or manifest.get("cython") != "3.3.0"
                or any(directives.get(key) is not False for key in ("infer_types", "annotation_typing", "cdivision"))
                or "-ffp-contract=off" not in manifest.get("flags", ())
                or "-fno-fast-math" not in manifest.get("flags", ())):
            return None
        item = manifest.get("modules", {}).get(name, {})
        binary_name = item.get("binary_name", "")
        if (not binary_name or Path(binary_name).name != binary_name
                or not any(binary_name.endswith(suffix) for suffix in importlib.machinery.EXTENSION_SUFFIXES)):
            return None
        binary = self.directory / binary_name
        source = self.root / (name + ".py")
        try:
            if hashlib.sha256(binary.read_bytes()).hexdigest() != item.get("binary_sha256"):
                return None
            if source.is_file():
                if hashlib.sha256(source.read_bytes()).hexdigest() != item.get("source_sha256"):
                    return None
            elif not getattr(sys, "frozen", False):
                return None
        except OSError:
            return None
        return importlib.util.spec_from_file_location(fullname, binary, loader=_PolicyLoader(fullname, str(binary)))


def install_native_policies():
    if not any(isinstance(finder, _PolicyFinder) for finder in sys.meta_path):
        sys.meta_path.insert(0, _PolicyFinder())


def policy_backend_status():
    return {name: str(getattr(getattr(sys.modules.get("sphingolipid_toolkit._vendor.lipidgate.ms2." + name), "__spec__", None), "origin", ""))
            for name in POLICY_MODULES if "sphingolipid_toolkit._vendor.lipidgate.ms2." + name in sys.modules}
