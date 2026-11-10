"""Product defaults shared by GUI, CLI, API and the 2D runner."""
from dataclasses import dataclass


@dataclass(frozen=True)
class SearchConfig:
    top_n: int = 3
    precursor_tolerance_ppm: float = 5.0
    fragment_tolerance_da: float | None = None
    fragment_tolerance_ppm: float | None = 15.0
    min_relative_intensity: float = 0.002
    min_total_score: float = 50.0
    ms1_precursor_association_ppm: float = 15.0
    ms1_precursor_confirmation_ppm: float = 3.0
    ms1_precursor_max_gap_seconds: float = 10.0


DEFAULT_SEARCH_CONFIG = SearchConfig()
