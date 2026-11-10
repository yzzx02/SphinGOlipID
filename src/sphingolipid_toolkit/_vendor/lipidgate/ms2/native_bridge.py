"""Verify and load the optional CPython numeric adapter, with a full fallback."""
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys


@lru_cache(maxsize=1)
def load_bridge():
    folder=Path(__file__).parent/'native_backend'
    try:
        manifest=json.loads((folder/'search_bridge.json').read_text(encoding='utf-8'))
        if (manifest.get('api_version')!=1 or manifest.get('python_abi')!=list(sys.version_info[:2])
                or manifest.get('platform')!=sys.platform or manifest.get('machine')!=platform.machine()
                or '-ffp-contract=off' not in manifest.get('flags',())
                or '-fno-fast-math' not in manifest.get('flags',())):
            return None
        name=manifest['binary_name']
        if Path(name).name!=name:return None
        binary=folder/name
        if hashlib.sha256(binary.read_bytes()).hexdigest()!=manifest.get('binary_sha256'):return None
        for source,wanted in manifest['source_sha256'].items():
            if Path(source).name!=source:return None
            path=folder/source
            if path.is_file():
                if hashlib.sha256(path.read_bytes()).hexdigest()!=wanted:return None
            elif not getattr(sys,'frozen',False):return None
        if set(manifest['source_sha256'])!={'search_bridge.cpp','search_core.cpp'}:return None
        spec=importlib.util.spec_from_file_location('sphingolipid_toolkit._vendor.lipidgate.ms2._search_bridge',binary)
        if spec is None or spec.loader is None:return None
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module
    except (OSError,ValueError,TypeError,KeyError,ImportError,AttributeError):
        return None
