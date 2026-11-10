"""Greedy, one-observed-peak-per-fragment matching for sorted spectra."""
from __future__ import annotations

import bisect
from typing import Callable, Sequence
from .models import ExperimentalPeak, FragmentRecord, FragmentMatch

def match_fragments(peaks: Sequence[ExperimentalPeak], fragments: Sequence[FragmentRecord],
                    window_for_mz: Callable[[float], float],
                    experimental_mz: Sequence[float] | None = None) -> list[FragmentMatch]:
    if experimental_mz is None:
        experimental_mz = [peak.mz for peak in peaks]
    matches: list[FragmentMatch] = []
    used_peak_indexes = set()
    for fragment in fragments:
        window_da = window_for_mz(fragment.mz)
        left = bisect.bisect_left(experimental_mz, fragment.mz - window_da)
        right = bisect.bisect_right(experimental_mz, fragment.mz + window_da)
        best_index = None
        best_peak = None
        best_error = None
        for peak_index in range(left, right):
            if peak_index in used_peak_indexes:
                continue
            peak = peaks[peak_index]
            error = abs(peak.mz - fragment.mz)
            if best_peak is None or peak.relative_intensity > best_peak.relative_intensity:
                best_peak = peak
                best_index = peak_index
                best_error = error
        if best_peak is not None and best_index is not None and best_error is not None:
            used_peak_indexes.add(best_index)
            matches.append(FragmentMatch(fragment=fragment, experimental_peak=best_peak, mz_error=best_error))
    return matches

