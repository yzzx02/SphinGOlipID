"""Project-level configuration objects for SphinGOlipID workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class SphinGOlipIDConfig:
    """Unified configuration shared by CLI, GUI, and tests.

    The MS2 workflow can keep raw spectrum text files and precursor Excel files
    in separate folders. When both live together, pass the same path for
    ``raw_ms2_dir`` and ``precursor_dir``.
    """

    raw_ms2_dir: Path
    precursor_dir: Path
    ms1_library_path: Path
    output_dir: Path
    file_indices: Sequence[int] = tuple(range(1, 7))
    text_pattern: str = "hilic-msms-{i}.txt"
    target_pattern: str = "hilic-msms-{i}.xlsx"
    result_pattern: str = "result_msms-{i}.xlsx"
    final_result_name: str = "ms2_annotation_results.xlsx"
    log_file_name: str = "run_log.txt"
    intermediate_dir_name: str = "intermediate"
    encoding: str = "GBK"
    fragment_ppm: float = 20.0
    min_fragment_intensity: float = 20.0
    min_matched_fragments: int = 2
    min_match_score: float = 0.35
    top_n: int = 3
    enable_rt_validation: bool = True
    enable_deduplication: bool = True
    save_intermediate: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "raw_ms2_dir", Path(self.raw_ms2_dir))
        object.__setattr__(self, "precursor_dir", Path(self.precursor_dir))
        object.__setattr__(self, "ms1_library_path", Path(self.ms1_library_path))
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        object.__setattr__(self, "file_indices", tuple(int(i) for i in self.file_indices))

    @property
    def fragment_tolerance_fraction(self) -> float:
        return float(self.fragment_ppm) / 1_000_000

    @classmethod
    def from_single_input_dir(
        cls,
        input_dir: str | Path,
        output_dir: str | Path,
        ms1_library_path: str | Path,
        **kwargs,
    ) -> "SphinGOlipIDConfig":
        """Build a config for the historical one-folder input layout."""

        input_path = Path(input_dir)
        return cls(
            raw_ms2_dir=input_path,
            precursor_dir=input_path,
            ms1_library_path=Path(ms1_library_path),
            output_dir=Path(output_dir),
            **kwargs,
        )
