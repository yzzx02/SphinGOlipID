from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True, slots=True)
class FragmentRecord:
    mz: float
    name: str
    fragment_type: str
    intensity: float = 100.0
    weight: float = 1.0
    required_group: Optional[str] = None

    def __setstate__(self, state):
        _restore_slot_state(self, state)


@dataclass(slots=True)
class LibraryRecord:
    record_id: int
    compound_class: str
    lipid_name: str
    lipid_chain_name: str
    precursor_mz: float
    adduct: str
    formula: str = ""
    polarity: str = "-"
    fragments: List[FragmentRecord] = field(default_factory=list)
    metadata: Dict[str, str] = field(default_factory=dict)

    def __setstate__(self, state):
        _restore_slot_state(self, state)

    @property
    def key(self) -> Tuple[str, str, str, float, str]:
        return (
            self.compound_class,
            self.lipid_name,
            self.lipid_chain_name,
            round(self.precursor_mz, 4),
            self.adduct,
        )


def _restore_slot_state(instance, state):
    """Read both earlier dict-backed caches and new slotted pickles."""
    if isinstance(state, dict):
        values = state.items()
    elif isinstance(state, tuple) and len(state) == 2 and isinstance(state[1], dict):
        values = state[1].items()
    else:
        values = zip((field.name for field in fields(instance)), state)
    for name, value in values:
        object.__setattr__(instance, name, value)


@dataclass(frozen=True, slots=True)
class ExperimentalPeak:
    mz: float
    intensity: float
    relative_intensity: float

    def __setstate__(self, state):
        _restore_slot_state(self, state)


@dataclass
class ExperimentalSpectrum:
    scan_id: str
    precursor_mz: float
    rt_minutes: float
    polarity: str
    peaks: List[ExperimentalPeak]
    metadata: Dict[str, object] = field(default_factory=dict)
    precursor_charge: Optional[int] = None

    @property
    def total_relative_intensity(self) -> float:
        return sum(peak.relative_intensity for peak in self.peaks)

    @property
    def base_peak_intensity(self) -> float:
        return max((peak.intensity for peak in self.peaks), default=0.0)

    @property
    def total_ion_intensity(self) -> float:
        return sum(peak.intensity for peak in self.peaks)


@dataclass(frozen=True, slots=True)
class FragmentMatch:
    fragment: FragmentRecord
    experimental_peak: ExperimentalPeak
    mz_error: float

    def __setstate__(self, state):
        _restore_slot_state(self, state)


# The optional bridge may allocate this exact frozen dataclass directly only
# while both its generated initializer and slot descriptors remain unchanged.
_FRAGMENT_MATCH_INITIALIZER = (FragmentMatch.__init__, FragmentMatch.__init__.__code__)


@dataclass(slots=True)
class PoolScore:
    pool_name: str
    matched_count: int
    total_count: int
    count_ratio: float
    intensity_ratio: float
    weight_ratio: float
    pool_score: float
    scoring_fragment_indexes: tuple[int, ...] = ()
    auxiliary_fragment_indexes: tuple[int, ...] = ()

    def __setstate__(self, state):
        _restore_slot_state(self, state)


@dataclass
class CandidateScore:
    record: LibraryRecord
    total_score: float
    passed_required_gates: bool
    missing_required_groups: List[str]
    ppm_error: float
    resolution_level: str
    matched_fragments: List[FragmentMatch] = field(default_factory=list)
    pool_scores: Dict[str, PoolScore] = field(default_factory=dict)
    matched_intensity_sum: float = 0.0
    matched_relative_intensity_sum: float = 0.0
    downgrade_reason: str = ""


def normalize_peaks(peaks: Sequence[Tuple[float, float]]) -> List[ExperimentalPeak]:
    if not peaks:
        return []
    max_intensity = max(intensity for _, intensity in peaks)
    if max_intensity <= 0:
        return []
    normalized = []
    for mz, intensity in sorted(peaks, key=lambda item: item[0]):
        rel = intensity / max_intensity
        normalized.append(ExperimentalPeak(mz=float(mz), intensity=float(intensity), relative_intensity=rel))
    return normalized
