"""Search independent mzML files with one library load per process."""

from __future__ import annotations

import ctypes
import gc
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path

_worker_searcher = None


def _available_memory_bytes():
    if os.name == "nt":
        class MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        status = MemoryStatus()
        status.dwLength = ctypes.sizeof(status)
        return int(status.ullAvailPhys) if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)) else None
    try:
        return int(os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return None


def _estimated_worker_memory_bytes(library_path):
    """Conservative estimate based on the parsed cache or source library size."""
    from .indexed_library import prebuilt_index_path
    from .library import _library_cache_path

    library_path = Path(library_path)
    if prebuilt_index_path(library_path) is not None:
        # Shared read-only disk index, small precursor arrays and bounded LRU.
        # Reserve space for native readers, results and compressed raw scans.
        return 768 * 1024**2
    cache_path = _library_cache_path(library_path)
    if cache_path.is_file() and cache_path.stat().st_size >= 100_000_000:
        return max(int(6.5 * 1024**3), cache_path.stat().st_size * 8)
    source_bytes = library_path.stat().st_size
    if library_path.suffix.lower() == ".gz" and source_bytes >= 5_000_000:
        # Slotted library records occupy about 5.9 GiB before scoring;
        # reserve additional working memory even when no cache exists yet.
        return max(int(6.5 * 1024**3), source_bytes * 250)
    if source_bytes >= 100_000_000:
        return max(3 * 1024**3, source_bytes * 9)
    return 0  # Small custom libraries need no large-library memory guard.


def validate_parallel_capacity(workers, file_count, library_path):
    if isinstance(workers, bool) or not isinstance(workers, int) or not 1 <= workers <= 4:
        raise ValueError("MS2 并行进程数必须是 1–4")
    effective = min(workers, file_count)
    if effective <= 1:
        return effective
    per_worker = _estimated_worker_memory_bytes(library_path)
    available = _available_memory_bytes()
    required = per_worker * effective + 3 * 1024**3
    if per_worker and available is not None and available < required:
        raise MemoryError(
            f"{effective} 个 MS2 进程预计需要约 {required / 1024**3:.1f} GB 空闲内存，"
            f"当前约 {available / 1024**3:.1f} GB；请减少进程数"
        )
    return effective


def _make_searcher(options):
    from .search import LipidMS2Searcher

    return LipidMS2Searcher(**options)


def _initialize_worker(options):
    global _worker_searcher
    _worker_searcher = _make_searcher(options)


def _search_file(task):
    path, top_n = task
    if _worker_searcher is None:
        raise RuntimeError("MS2 worker was not initialized")
    return _worker_searcher.search_mzml(path, top_n=top_n)


def search_files(paths, *, search_options, top_n, workers=1, on_file_done=None, on_library_ready=None):
    """Return per-file frames in input order; never share native mzML state."""
    paths = [Path(path) for path in paths]
    effective_workers = validate_parallel_capacity(workers, len(paths), search_options.get("library_path"))
    if effective_workers <= 1:
        searcher = _make_searcher(search_options)
        try:
            if on_library_ready is not None:
                on_library_ready()
            results = []
            for path in paths:
                results.append(searcher.search_mzml(path, top_n=top_n))
                if on_file_done is not None:
                    on_file_done(path, len(results), len(paths))
            return results
        finally:
            searcher.close()
    # A cold cache must be built once before workers start; otherwise every
    # spawned process parses and serializes the same large MSP independently.
    from .library import _library_cache_valid, _install_bundled_prebuilt_cache, load_library
    from .indexed_library import prebuilt_index_path

    library_path = Path(search_options["library_path"])
    if (prebuilt_index_path(library_path) is None and not _library_cache_valid(library_path)
            and not _install_bundled_prebuilt_cache(library_path)):
        records = load_library(library_path)
        del records
        gc.collect()
    if on_library_ready is not None:
        on_library_ready()
    # Spawned workers load the library once each; only file paths and result
    # frames cross process boundaries. Completion order may differ from input
    # order, so place each frame back at its original index.
    with ProcessPoolExecutor(
        max_workers=effective_workers,
        mp_context=get_context("spawn"),
        initializer=_initialize_worker,
        initargs=(search_options,),
    ) as pool:
        pending = {
            pool.submit(_search_file, (path, top_n)): (index, path)
            for index, path in enumerate(paths)
        }
        results = [None] * len(paths)
        for completed, future in enumerate(as_completed(pending), 1):
            index, path = pending[future]
            results[index] = future.result()
            if on_file_done is not None:
                on_file_done(path, completed, len(paths))
        return results
