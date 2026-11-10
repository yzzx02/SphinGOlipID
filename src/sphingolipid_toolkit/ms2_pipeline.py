"""Text batch adapters for the sole current sphingolipid identification engine."""
from __future__ import annotations
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple
import pandas as pd

from .backend_policy import removed_backend
from .config import TEXT_FRAGMENT_PPM, SphinGOlipIDConfig, normalize_file_indices, validate_batch_settings
from .io_utils import require_file, resolve_ms2_input_paths, validate_output_paths
from .logging_utils import get_logger, attach_file_handler
from .result_writer import write_identifications
from .scoring import add_standard_score_columns

@dataclass
class MS2PipelineConfig:
    """Historical constructor name, always backed by the current engine."""
    input_dir: Path
    output_dir: Path
    ms1_db_path: Path | None = None
    precursor_dir: Path | None = None
    file_indices: Sequence[int] = tuple(range(1, 7))
    text_pattern: str = 'hilic-msms-{i}.txt'
    target_pattern: str = 'hilic-msms-{i}.xlsx'
    result_pattern: str = 'result_msms-{i}.xlsx'
    final_result_name: str = 'ms2_annotation_results.xlsx'
    log_file_name: str = 'run_log.txt'
    intermediate_dir_name: str = 'intermediate'
    encoding: str = 'GBK'
    fragment_ppm: float = TEXT_FRAGMENT_PPM
    min_fragment_intensity: float | None = None
    min_matched_fragments: int | None = None
    min_match_score: float | None = None
    top_n: int = 3
    keep_intermediate: bool = False
    backend: str = 'lipidgate'
    ms2_library_path: Path | None = None
    precursor_ppm: float = 10.0
    fragment_da: float | None = None
    min_total_score: float = 50.0
    feature_rt_tolerance_min: float = 0.5

    def __post_init__(self):
        self.input_dir=Path(self.input_dir);self.output_dir=Path(self.output_dir)
        if self.ms1_db_path is not None:self.ms1_db_path=Path(self.ms1_db_path)
        if self.precursor_dir is not None:self.precursor_dir=Path(self.precursor_dir)
        if self.ms2_library_path is not None:self.ms2_library_path=Path(self.ms2_library_path)
        self.file_indices=normalize_file_indices(self.file_indices)
        validate_batch_settings(self)

    @property
    def raw_ms2_dir(self):return self.input_dir
    @property
    def ms1_library_path(self):return self.ms1_db_path
    @property
    def save_intermediate(self):return self.keep_intermediate
    @property
    def fragment_tolerance_fraction(self):return self.fragment_ppm/1_000_000

ConfigLike=MS2PipelineConfig | SphinGOlipIDConfig

def _coerce_config(config: ConfigLike) -> MS2PipelineConfig:
    validate_batch_settings(config)
    if isinstance(config,MS2PipelineConfig):return config
    return MS2PipelineConfig(input_dir=config.raw_ms2_dir,output_dir=config.output_dir,
        ms1_db_path=config.ms1_library_path,precursor_dir=config.precursor_dir,
        file_indices=config.file_indices,text_pattern=config.text_pattern,target_pattern=config.target_pattern,
        result_pattern=config.result_pattern,final_result_name=config.final_result_name,
        log_file_name=config.log_file_name,intermediate_dir_name=config.intermediate_dir_name,
        encoding=config.encoding,fragment_ppm=config.fragment_ppm,top_n=config.top_n,
        keep_intermediate=config.save_intermediate,backend=config.backend,ms2_library_path=config.ms2_library_path,
        precursor_ppm=config.precursor_ppm,fragment_da=config.fragment_da,
        min_total_score=config.min_total_score,feature_rt_tolerance_min=config.feature_rt_tolerance_min)

def run_one_file(config: ConfigLike,file_index:int,logger=None) -> Path:
    from .reference_pipeline import run_text_file
    config=_coerce_config(config)
    validate_output_paths(config,(file_index,))
    paths=resolve_ms2_input_paths(config,file_index)
    require_file(paths.text_file,'MS/MS txt file');require_file(paths.target_file,'precursor target Excel')
    return run_text_file(config,file_index,logger)

def run_batch(config: ConfigLike,logger=None) -> List[Path]:
    from .engine import library_path
    config=_coerce_config(config)
    validate_output_paths(config,config.file_indices)
    require_file(config.ms2_library_path or library_path('+'),'MS2 sphingolipid library')
    for index in config.file_indices:
        paths=resolve_ms2_input_paths(config,index)
        require_file(paths.text_file,f'MS/MS txt file for index {index}')
        require_file(paths.target_file,f'precursor target Excel for index {index}')
    logger=logger or get_logger('sphingolipid_toolkit.ms2')
    output=Path(config.output_dir);output.mkdir(parents=True,exist_ok=True)
    handler=attach_file_handler(output/config.log_file_name,logger=logger)
    try:
        logger.info('Starting MS2 batch with current sphingolipid engine')
        results=[run_one_file(config,i,logger) for i in config.file_indices]
        _write_batch_result(results,output/config.final_result_name,logger)
        logger.info('MS2 batch completed')
        return results
    finally:
        logger.removeHandler(handler);handler.close()

def _write_batch_result(result_files:Sequence[Path],final_result_file:Path,logger) -> pd.DataFrame:
    frames=[]
    for value in result_files:
        path=Path(value)
        if not path.exists():continue
        frame=pd.read_excel(path)
        if 'evidence_file' in frame:
            frame['evidence_file']=frame['evidence_file'].map(lambda name:
                os.path.relpath(path.parent/str(name),Path(final_result_file).parent) if pd.notna(name) else name)
        frame['source_result_file']=path.name
        frames.append(frame)
    result=add_standard_score_columns(pd.concat(frames,ignore_index=True) if frames else pd.DataFrame())
    logger.info('Exporting results: %s',final_result_file)
    write_identifications(result,final_result_file)
    return result

# Fail explicitly for callers of the retired candidate-generation/coverage route.
_write_db_sheets=removed_backend
_finalize_results=removed_backend

def parse_file_indices(value: str) -> Tuple[int, ...]:
    """Parse values like ``1-6`` or ``1,3,5`` for the CLI."""

    value = value.strip()
    if re.fullmatch(r"\d+-\d+", value):
        a, b = map(int, value.split("-"))
        if b < a:
            raise ValueError("file index range end must be >= start")
        return tuple(range(a, b + 1))
    indices = tuple(int(x.strip()) for x in value.split(",") if x.strip())
    if not indices or any(i < 0 for i in indices):
        raise ValueError("file indices must contain at least one non-negative integer")
    return tuple(dict.fromkeys(indices))
