"""Project-level configuration objects for SphinGOlipID workflows."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from pathlib import Path
from typing import Sequence


# Historical workflow defaults intentionally differ; neither is a manuscript
# final-data assertion. All matching tolerances are supplied through config.
LEGACY_FRAGMENT_PPM = 20.0
TARGETED_FRAGMENT_PPM = 10.0


def normalize_file_indices(values: Sequence[int]) -> tuple[int, ...]:
    """Accept integral identifiers without silently truncating fractional inputs."""

    indices = []
    for value in values:
        try:
            if isinstance(value, bool) or (isinstance(value, Real) and (not isfinite(value) or value != int(value))):
                raise ValueError
            indices.append(int(value))
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"File index must be an integer: {value!r}") from exc
    return tuple(indices)


def validate_matching_parameters(config: object) -> None:
    """Reject settings that silently produce empty or invalid annotations."""

    for name, minimum, maximum in (
        ("fragment_ppm", 0.0, None),
        ("min_fragment_intensity", 0.0, None),
        ("min_match_score", 0.0, 1.0),
    ):
        value = float(getattr(config, name))
        if (
            not isfinite(value)
            or value < minimum
            or (name == "fragment_ppm" and value == 0)
            or (maximum is not None and value > maximum)
        ):
            raise ValueError(f"Invalid {name}: {value}")
    for name in ("top_n", "min_matched_fragments"):
        value = getattr(config, name)
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{name} must be a positive integer")


def validate_batch_settings(config: object) -> None:
    validate_matching_parameters(config)
    if not config.file_indices or any(i < 0 for i in config.file_indices):
        raise ValueError("file_indices must contain at least one non-negative integer")
    if len(set(config.file_indices)) != len(config.file_indices):
        raise ValueError("file_indices must not contain duplicates")
    try:
        result_names = [config.result_pattern.format(i=i) for i in config.file_indices]
        for i in config.file_indices:
            config.text_pattern.format(i=i)
            config.target_pattern.format(i=i)
    except (KeyError, ValueError, IndexError) as exc:
        raise ValueError("File patterns must use only the {i} index placeholder; escape literal braces") from exc
    if len(set(result_names)) != len(result_names):
        raise ValueError("result_pattern must produce a distinct output for each file index")


@dataclass(frozen=True)
class SphinGOlipIDConfig:
    """Unified configuration shared by CLI, GUI, and tests.

    The MS2 workflow can keep raw spectrum text files and precursor Excel files
    in separate folders. When both live together, pass the same path for
    ``raw_ms2_dir`` and ``precursor_dir``.

    ``enable_rt_validation`` and ``enable_deduplication`` are reserved flags;
    ``run_batch`` currently performs MS2 matching only. Call the separate RT
    validation and result-cleaning APIs when those stages are required.
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
    fragment_ppm: float = LEGACY_FRAGMENT_PPM
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
        object.__setattr__(self, "file_indices", normalize_file_indices(self.file_indices))
        validate_batch_settings(self)

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
