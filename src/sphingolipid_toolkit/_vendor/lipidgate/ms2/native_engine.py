"""Optional batched C++ matching over bounded, unchanged library blocks.

Only numeric overlap/matching moves across the C ABI. All evidence gates,
scores, candidate order and export provenance remain in the Python scorer.
Missing/incompatible binaries and unsupported custom inputs use Python.
"""
from __future__ import annotations

from array import array
import ctypes as ct
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import sys

import numpy as np

from .indexed_library import IndexedLibrary
from .models import FragmentMatch
from . import native_bridge
from .record_facts import get_native_query


def _binary_path():
    name = {"win32": "search_core.dll", "darwin": "search_core.dylib"}.get(sys.platform, "search_core.so")
    return Path(__file__).parent / "native_backend" / name


@lru_cache(maxsize=1)
def _kernel():
    binary = _binary_path()
    try:
        manifest = json.loads(binary.with_suffix(".json").read_text(encoding="utf-8"))
        if (manifest.get("api_version") != 1 or manifest.get("platform") != sys.platform
                or manifest.get("machine") != platform.machine()
                or manifest.get("binary_sha256") != hashlib.sha256(binary.read_bytes()).hexdigest()):
            return None
        source = binary.with_suffix(".cpp")
        if source.is_file() and manifest.get("source_sha256") != hashlib.sha256(source.read_bytes()).hexdigest():
            return None
        dll = ct.CDLL(str(binary))
        dll.lg_api_version.argtypes = []
        dll.lg_api_version.restype = ct.c_uint32
        if dll.lg_api_version() != 1:
            return None
        doubles, bytes_, uints = ct.POINTER(ct.c_double), ct.POINTER(ct.c_uint8), ct.POINTER(ct.c_uint32)
        dll.lg_batch_match.argtypes = [
            doubles, doubles, ct.c_uint32, doubles, bytes_, ct.c_uint32,
            uints, ct.c_uint32, uints, bytes_, ct.c_uint32,
            uints, ct.c_uint32, ct.c_double, ct.c_uint8, ct.c_uint8,
            bytes_, ct.POINTER(ct.c_int32),
        ]
        dll.lg_batch_match.restype = ct.c_int
        return dll
    except (OSError, ValueError, AttributeError, TypeError):
        return None


def native_available():
    return _kernel() is not None


def _pointer(values, c_type):
    if not values:
        return ct.POINTER(c_type)()
    return ct.cast((c_type * len(values)).from_buffer(values), ct.POINTER(c_type))


def _supported_number(value):
    # Preserve custom numeric semantics by falling back rather than coercing
    # float32/Decimal values and silently changing an endpoint comparison.
    return type(value) in (float, int) and math.isfinite(value) and float(value) == value


def _supported_intensities(peaks):
    # pyOpenMS supplies float32 relative intensities. Widening float32 values
    # preserves their ordering exactly; no intensity arithmetic moves to C++.
    # NumPy's mixed Python-float/float32 weak promotion can round comparisons,
    # so keep that custom mixed case on Python rather than changing tie order.
    types = {type(p.relative_intensity) for p in peaks}
    if not types.issubset({float, int, np.float32, np.float64}):
        return False
    if np.float32 in types and not types.issubset({np.float32, np.float64}):
        return False
    return all(math.isfinite(p.relative_intensity) and
               (type(p.relative_intensity) is not int or float(p.relative_intensity) == p.relative_intensity)
               for p in peaks)


class _NumericBlock:
    __slots__ = ("signature", "masses", "fragment_flags", "occurrences", "offsets", "record_flags")

    def __init__(self, block, signature):
        glyceride_classes, loss_types, sphingo_keys = signature
        columns = block.column_buffers
        if columns is not None:
            if columns["masses"] is None:
                raise ValueError("Unsupported numeric block")
            self.signature = signature
            self.masses = columns["masses"]
            flags = [int(t == "FA_Frag" or t in loss_types) | (2 * int(t == "Precursor Ion"))
                     for t in columns["type_names"]]
            self.fragment_flags = array("B", (flags[i] for i in columns["type_ids"]))
            pairs = []
            for compound_class, adduct in columns["class_adduct_names"]:
                key = re.sub(r"[^A-Za-z0-9]+", "", str(compound_class or "").upper())
                glyceride = str(adduct or "").strip().endswith("+") and key in glyceride_classes
                fa = str(compound_class or "").strip().upper() == "FA" and f"{compound_class}_{adduct}" not in sphingo_keys
                pairs.append(int(glyceride) | (2 * int(fa)))
            self.record_flags = array("B", (pairs[i] for i in columns["class_adduct_ids"]))
            self.occurrences = (columns["occurrences"] if columns["occurrences"].typecode == "I"
                                else array("I", columns["occurrences"]))
            self.offsets = columns["offsets"]
            return
        if block.templates is None or any(not _supported_number(f[0]) for f in block.templates):
            raise ValueError("Unsupported numeric block")
        self.signature = signature
        self.masses = array("d", (f[0] for f in block.templates))
        self.fragment_flags = array("B", (
            int(f[2] == "FA_Frag" or f[2] in loss_types) | (2 * int(f[2] == "Precursor Ion"))
            for f in block.templates
        ))
        self.occurrences, self.offsets, self.record_flags = array("I"), array("I", [0]), array("B")
        pair_flags = {}
        for i, value in enumerate(block.values):
            pair = (value[1], value[5])
            flag = pair_flags.get(pair)
            if flag is None:
                class_key = re.sub(r"[^A-Za-z0-9]+", "", str(value[1] or "").upper())
                glyceride = str(value[5] or "").strip().endswith("+") and class_key in glyceride_classes
                fa = (str(value[1] or "").strip().upper() == "FA"
                      and f"{value[1]}_{value[5]}" not in sphingo_keys)
                pair_flags[pair] = flag = int(glyceride) | (2 * int(fa))
            self.record_flags.append(flag)
            self.occurrences.extend(block.occurrence_ids(i))
            self.offsets.append(len(self.occurrences))


def _records_unchanged(block, requested):
    # LibraryRecord is mutable. Do not use a stale packed mass/class view if a
    # caller edited an already decoded candidate in memory.
    for offset in requested:
        record = block.records.get(offset)
        if record is None:
            continue
        value = block.values[offset]
        occurrences = block.occurrence_ids(offset)
        if record.compound_class != value[1] or record.adduct != value[5] or len(record.fragments) != len(occurrences):
            return False
        for fragment, occurrence in zip(record.fragments, occurrences):
            template = block.fragment_template(occurrence)
            if (type(fragment.mz) is not type(template[0]) or fragment.mz != template[0]
                    or fragment.fragment_type != template[2]):
                return False
    return True


def prepare_indexed_matches(library, spectrum, left, right, *, tolerance_da, tolerance_ppm,
                            prefilter, glyceride_classes, loss_types, sphingo_keys):
    """Return ordered visible indexes -> occurrence peak indexes, or Python fallback."""
    if not isinstance(library, IndexedLibrary):
        return None
    # Search and generic scoring have different historical defaults when both
    # tolerances are None; that unusual configuration keeps its Python path.
    tolerance = tolerance_da if tolerance_da is not None else tolerance_ppm
    if tolerance is None or not _supported_number(tolerance) or tolerance < 0:
        return None
    kernel = _kernel()
    if kernel is None:
        return None
    if left >= right:
        return {}
    if any(not _supported_number(p.mz) for p in spectrum.peaks) or not _supported_intensities(spectrum.peaks):
        return None
    mz = array("d", (p.mz for p in spectrum.peaks))
    if any(a > b for a, b in zip(mz, mz[1:])):
        return None
    intensity = array("d", (p.relative_intensity for p in spectrum.peaks))
    signature = (frozenset(glyceride_classes), frozenset(loss_types), frozenset(sphingo_keys))
    result = {}
    for block, indexes, local_indexes in library.candidate_record_blocks(left, right):
        if block.templates is None or not _records_unchanged(block, local_indexes):
            return None
        numeric = block.numeric_view
        if numeric is None or numeric.signature != signature:
            try:
                numeric = _NumericBlock(block, signature)
            except (ValueError, TypeError, OverflowError):
                return None
            block.numeric_view = numeric
        requested = array("I", local_indexes)
        selected = array("B", [0]) * len(requested)
        matches = array("i", [-1]) * len(numeric.occurrences)
        status = kernel.lg_batch_match(
            _pointer(mz, ct.c_double), _pointer(intensity, ct.c_double), len(mz),
            _pointer(numeric.masses, ct.c_double), _pointer(numeric.fragment_flags, ct.c_uint8), len(numeric.masses),
            _pointer(numeric.occurrences, ct.c_uint32), len(numeric.occurrences),
            _pointer(numeric.offsets, ct.c_uint32), _pointer(numeric.record_flags, ct.c_uint8), len(block.values),
            _pointer(requested, ct.c_uint32), len(requested), float(tolerance), tolerance_da is None, bool(prefilter),
            _pointer(selected, ct.c_uint8), _pointer(matches, ct.c_int32),
        )
        if status:
            raise RuntimeError(f"Native fragment matching failed (status {status})")
        for q, local in enumerate(local_indexes):
            if selected[q]:
                result[indexes[q]] = tuple(matches[numeric.offsets[local]:numeric.offsets[local + 1]])
    return result


def materialize_matches(spectrum, record, peak_indexes):
    """Retain original fragment/peak objects and Python's exact error value."""
    bridge = native_bridge.load_bridge()
    query = get_native_query(spectrum)
    if query is not None and bridge is not None and hasattr(bridge, 'materialize_query_matches'):
        native = bridge.materialize_query_matches(query, record.fragments, peak_indexes, FragmentMatch)
        if native is not None:
            return native
    if bridge is not None and hasattr(bridge, 'materialize_matches'):
        native = bridge.materialize_matches(spectrum.peaks, record.fragments, peak_indexes, FragmentMatch)
        if native is not None:
            return native
    return [FragmentMatch(fragment, spectrum.peaks[index], abs(spectrum.peaks[index].mz - fragment.mz))
            for fragment, index in zip(record.fragments, peak_indexes) if index >= 0]
