"""Load bounded numeric search columns and run conservative C++ evidence plans."""
from __future__ import annotations

from array import array
import bisect
import ctypes as ct
from itertools import accumulate, islice
import pickle
import sqlite3
import sys
import zlib

from . import safe_pickle, scoring
from .indexed_library import IndexedLibrary
from .native_engine import _kernel, _pointer, _supported_number, _supported_intensities
from .search_plan import PLAN_VERSION, read_programs
from .record_facts import record_revision, record_unchanged, _current_spectrum_facts, set_native_query


MAX_PLAN_BYTES = 256 * 1024 * 1024


def _column(value, kind):
    if type(value) is not array or value.typecode != kind:
        raise ValueError('Invalid search plan column')
    return value


def _monotone_offsets(value, length):
    _column(value, 'I')
    if not value or value[0] != 0 or value[-1] != length or any(a > b for a, b in zip(value, value[1:])):
        raise ValueError('Invalid search plan offsets')


class NativeSearchPlan:
    def __init__(self, library, kernel):
        self.library, self.kernel = library, kernel
        if kernel.lg_plan_api_version() != PLAN_VERSION:
            raise ValueError('Unsupported native search plan ABI')
        doubles = ct.POINTER(ct.c_double)
        bytes_ = ct.POINTER(ct.c_uint8)
        uints = ct.POINTER(ct.c_uint32)
        kernel.lg_search_plan.argtypes = [
            doubles, doubles, ct.c_uint32, doubles, bytes_, uints, bytes_, uints, ct.c_uint32,
            uints, ct.c_uint32, uints, ct.POINTER(ct.c_uint64), bytes_,
            ct.c_uint32, ct.c_uint32, ct.c_double, ct.c_uint8, ct.c_uint8,
            bytes_, ct.POINTER(ct.c_int32),
        ]
        kernel.lg_search_plan.restype = ct.c_int
        self.plan_offsets, self.arm_offsets, self.masks, self.minimums = read_programs(library._connection)
        _column(self.masks, 'Q'); _column(self.minimums, 'B')
        if len(self.masks) != len(self.minimums) or any(n > 64 for n in self.minimums):
            raise ValueError('Invalid search plan clauses')
        _monotone_offsets(self.arm_offsets, len(self.masks))
        _monotone_offsets(self.plan_offsets, len(self.arm_offsets) - 1)
        self.masses, self.fragment_flags = array('d'), array('B')
        self.offsets, self.record_flags, self.programs = array('I', [0]), array('B'), array('I')
        positions, block_size = library._positions, library.metadata['block_size']
        first = positions[0] // block_size if positions else 0
        last = positions[-1] // block_size if positions else -1
        cursor = library._connection.execute(
            'SELECT block_id,payload FROM search_plans WHERE block_id>=? AND block_id<=? ORDER BY block_id', (first, last))
        visited = 0
        for block_id, payload in cursor:
            left = bisect.bisect_left(positions, block_id * block_size)
            right = bisect.bisect_left(positions, (block_id + 1) * block_size, left)
            if left == right:
                continue
            masses, flags, counts, record_flags, programs = safe_pickle.loads(zlib.decompress(payload))
            for value, kind in ((masses, 'd'), (flags, 'B'), (counts, 'I'), (record_flags, 'B'), (programs, 'I')):
                _column(value, kind)
            if (len(masses) != len(flags) or sum(counts) != len(masses)
                    or len(counts) != len(record_flags) or len(counts) != len(programs)
                    or len(counts) > block_size
                    or any(p >= len(self.plan_offsets) - 1 for p in programs)
                    or any(f > 3 for f in flags) or any(f > 7 for f in record_flags)):
                raise ValueError('Invalid numeric search block')
            base = self.offsets[-1]
            if right - left == len(counts) and positions[left] == block_id * block_size:
                self.masses.extend(masses); self.fragment_flags.extend(flags)
                self.record_flags.extend(record_flags); self.programs.extend(programs)
                self.offsets.extend(islice(accumulate(counts, initial=base), 1, None))
            else:
                boundaries = array('I', accumulate(counts, initial=0))
                for visible in range(left, right):
                    local = positions[visible] - block_id * block_size
                    if local >= len(counts):
                        raise ValueError('Search block is missing a visible record')
                    a, b = boundaries[local], boundaries[local+1]
                    self.masses.extend(masses[a:b]); self.fragment_flags.extend(flags[a:b])
                    self.offsets.append(len(self.masses))
                    self.record_flags.append(record_flags[local]); self.programs.append(programs[local])
            visited += right - left
            if self.bytes > MAX_PLAN_BYTES:
                raise ValueError('Numeric search plan exceeds the memory budget')
        if visited != len(positions) or len(self.offsets) != len(positions) + 1:
            raise ValueError('Incomplete numeric search plan')
        self._masses_ptr = _pointer(self.masses, ct.c_double)
        self._flags_ptr = _pointer(self.fragment_flags, ct.c_uint8)
        self._offsets_ptr = _pointer(self.offsets, ct.c_uint32)
        self._record_flags_ptr = _pointer(self.record_flags, ct.c_uint8)
        self._programs_ptr = _pointer(self.programs, ct.c_uint32)
        self._plan_offsets_ptr = _pointer(self.plan_offsets, ct.c_uint32)
        self._arm_offsets_ptr = _pointer(self.arm_offsets, ct.c_uint32)
        self._masks_ptr = _pointer(self.masks, ct.c_uint64)
        self._minimums_ptr = _pointer(self.minimums, ct.c_uint8)
        from .equivalent_candidates import load_equivalent_candidates
        self.equivalent = load_equivalent_candidates(library)
        if self.bytes > MAX_PLAN_BYTES:
            raise ValueError('Combined numeric search plan exceeds the memory budget')
        from .native_bridge import load_bridge
        self.bridge = load_bridge()
        self.bridge_plan = (None if self.bridge is None else self.bridge.create_plan(
            self.masses, self.fragment_flags, self.offsets, self.record_flags, self.programs,
            self.plan_offsets, self.arm_offsets, self.masks, self.minimums))

    @property
    def bytes(self):
        return sum(len(value) * value.itemsize for value in (
            self.masses, self.fragment_flags, self.offsets, self.record_flags,
            self.programs, self.plan_offsets, self.arm_offsets, self.masks, self.minimums)) + (
                self.equivalent.bytes if getattr(self, 'equivalent', None) is not None else 0)

    def _unchanged(self, left, right):
        """Honor edits to any currently decoded candidate before pruning it."""
        positions, size = self.library._positions, self.library.metadata['block_size']
        if left == right:
            return True
        low, high = positions[left], positions[right-1]
        selected = self.library._selected_positions
        a, b = bisect.bisect_left(selected, low), bisect.bisect_right(selected, high)
        for i in range(a, b):
            entry = self.library._selected_records[selected[i]]
            if not record_unchanged(entry[0], entry[2]):
                return False
        for block_id, block in self.library._blocks.items():
            if block_id * size > high or (block_id + 1) * size <= low:
                continue
            if not block.packed or block.templates is None:
                return False
            for local, record in block.records.items():
                position = block_id * size + local
                if position < low or position > high:
                    continue
                original = block.values[local]
                occurrences = block.occurrence_ids(local)
                if (record.compound_class != original[1] or record.lipid_chain_name != original[3]
                        or record.lipid_name != original[2] or record.adduct != original[5]
                        or len(record.fragments) != len(occurrences)):
                    return False
                for fragment, occurrence in zip(record.fragments, occurrences):
                    template = block.fragment_template(occurrence)
                    if (type(fragment.mz) is not type(template[0]) or fragment.mz != template[0]
                            or fragment.name != template[1] or fragment.fragment_type != template[2]):
                        return False
        return True

    def match(self, spectrum, left, right, *, tolerance_da, tolerance_ppm, prefilter, collapse_equivalent=True):
        tolerance = tolerance_da if tolerance_da is not None else tolerance_ppm
        if (tolerance is None or not _supported_number(tolerance) or tolerance < 0
                or not self._unchanged(left, right)):
            return None
        if self.bridge_plan is not None:
            current = _current_spectrum_facts.get()
            shared = (self.bridge.match_query_plan(self.bridge_plan, spectrum, left, right,
                      float(tolerance), tolerance_da is None, bool(prefilter), scoring.PRECURSOR_CLUSTER_EXCLUSION_DA)
                      if current is not None and current[0] is spectrum
                      and hasattr(self.bridge, 'match_query_plan') else None)
            if shared is not None:
                result, query = shared
                set_native_query(spectrum, query, self.bridge, scoring.PRECURSOR_CLUSTER_EXCLUSION_DA)
            else:
                result = self.bridge.match_plan(self.bridge_plan, spectrum.peaks, left, right,
                                               float(tolerance), tolerance_da is None, bool(prefilter))
            if result is not None:
                return (result if self.equivalent is None or not collapse_equivalent else
                        self.equivalent.collapse(spectrum, result))
        if (any(not _supported_number(p.mz) for p in spectrum.peaks)
                or not _supported_intensities(spectrum.peaks)):
            return None
        mz = array('d', (p.mz for p in spectrum.peaks))
        if any(a > b for a, b in zip(mz, mz[1:])):
            return None
        intensity = array('d', (p.relative_intensity for p in spectrum.peaks))
        selected = array('B', [0]) * (right - left)
        start, stop = self.offsets[left], self.offsets[right]
        matches = array('i', [-1]) * (stop - start)
        status = self.kernel.lg_search_plan(
            _pointer(mz, ct.c_double), _pointer(intensity, ct.c_double), len(mz),
            self._masses_ptr, self._flags_ptr, self._offsets_ptr, self._record_flags_ptr,
            self._programs_ptr, len(self.record_flags), self._plan_offsets_ptr, len(self.plan_offsets)-1,
            self._arm_offsets_ptr, self._masks_ptr, self._minimums_ptr,
            left, right, float(tolerance), tolerance_da is None, bool(prefilter),
            _pointer(selected, ct.c_uint8), _pointer(matches, ct.c_int32))
        if status:
            raise RuntimeError(f'Native evidence plan failed (status {status})')
        result = {}
        for visible in range(left, right):
            flag = selected[visible-left]
            if flag:
                result[visible] = (None if flag == 2 else
                                   tuple(matches[self.offsets[visible]-start:self.offsets[visible+1]-start]))
        return (result if self.equivalent is None or not collapse_equivalent else
                self.equivalent.collapse(spectrum, result))


def load_native_search_plan(library):
    if (not isinstance(library, IndexedLibrary) or sys.byteorder != 'little'
            or library.metadata.get('search_plan_version') != PLAN_VERSION):
        return None
    kernel = _kernel()
    if kernel is None or not hasattr(kernel, 'lg_plan_api_version'):
        return None
    try:
        kernel.lg_plan_api_version.argtypes = []
        kernel.lg_plan_api_version.restype = ct.c_uint32
        return NativeSearchPlan(library, kernel)
    except (ValueError, TypeError, OverflowError, AttributeError, pickle.UnpicklingError,
            EOFError, sqlite3.DatabaseError, zlib.error):
        # A legacy/unsupported cache never changes the scientific search path.
        return None
