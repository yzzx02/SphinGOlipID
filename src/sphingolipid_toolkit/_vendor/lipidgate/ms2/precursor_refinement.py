"""Conservative, library-independent MS1 confirmation of an MS2 precursor m/z."""
from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from typing import Sequence
import numpy as np


@dataclass(frozen=True)
class MS1Survey:
    native_id: str
    scan_id: str
    rt_seconds: float
    mz: Sequence[float]
    intensity: Sequence[float]
    centroided: bool = True


@dataclass(frozen=True)
class PrecursorRefinement:
    raw_mz: float
    matching_mz: float
    status: str
    source: str = 'MS2_SELECTED_ION'
    ms1_mz: float | None = None
    ms1_scan_id: str | None = None
    ms1_native_id: str | None = None
    ms1_rt_seconds: float | None = None
    confirmation_count: int = 0


class MS1PrecursorRefiner:
    """Use a unique parent-survey centroid confirmed in an adjacent survey.

    No library mass is consulted. Ambiguity, profile data, missing support,
    stale parent references, or multiple nearby centroids retain the raw m/z.
    The association window is NOT the downstream library matching tolerance.
    """
    def __init__(self, surveys: Sequence[MS1Survey], *, association_ppm: float = 15.,
                 confirmation_ppm: float = 3., max_gap_seconds: float = 10.):
        self.surveys = sorted(surveys, key=lambda s:s.rt_seconds)
        self.times = [s.rt_seconds for s in self.surveys]
        self.by_native = {s.native_id:i for i,s in enumerate(self.surveys) if s.native_id}
        self.association_ppm = association_ppm
        self.confirmation_ppm = confirmation_ppm
        self.max_gap_seconds = max_gap_seconds

    @staticmethod
    def _near(survey: MS1Survey, mz: float, ppm: float) -> list[float]:
        delta=abs(mz)*ppm*1e-6
        left=int(np.searchsorted(survey.mz,mz-delta,side='left'))
        right=int(np.searchsorted(survey.mz,mz+delta,side='right'))
        return [float(survey.mz[i]) for i in range(left,right)
                if np.isfinite(survey.mz[i]) and np.isfinite(survey.intensity[i]) and survey.intensity[i]>0]

    def refine(self, raw_mz: float, rt_seconds: float, parent_native_id: str = '') -> PrecursorRefinement:
        def fallback(status):return PrecursorRefinement(raw_mz,raw_mz,status)
        if not np.isfinite(raw_mz) or raw_mz<=0:return fallback('invalid_raw_mz')
        if not self.surveys:return fallback('no_ms1')
        if parent_native_id:
            index=self.by_native.get(parent_native_id)
            if index is None:return fallback('referenced_ms1_missing')
            source='MS1_REFERENCE_CENTROID'
        else:
            index=bisect_left(self.times,rt_seconds)-1
            if index<0:return fallback('no_preceding_ms1')
            source='MS1_PRECEDING_CENTROID'
        parent=self.surveys[index]
        if not 0<=rt_seconds-parent.rt_seconds<=self.max_gap_seconds:return fallback('ms1_time_gap')
        if not parent.centroided:return fallback('ms1_not_centroided')
        candidates=self._near(parent,raw_mz,self.association_ppm)
        if not candidates:return fallback('no_ms1_peak_in_association_window')
        if len(candidates)!=1:return fallback('ambiguous_ms1_peaks')
        chosen=candidates[0]
        confirmations=0
        for neighbor_index in (index-1,index+1):
            if not 0<=neighbor_index<len(self.surveys):continue
            neighbor=self.surveys[neighbor_index]
            if not neighbor.centroided or abs(neighbor.rt_seconds-rt_seconds)>self.max_gap_seconds:continue
            if len(self._near(neighbor,chosen,self.confirmation_ppm))==1:confirmations+=1
        if not confirmations:return fallback('ms1_peak_not_confirmed')
        return PrecursorRefinement(raw_mz,chosen,'confirmed',source,chosen,parent.scan_id,
                                   parent.native_id,parent.rt_seconds,confirmations)

