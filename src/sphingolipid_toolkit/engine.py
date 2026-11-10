"""Sphingolipid-only facade for the verified LipidGate reference engine.

All matching, pool assignment, admission, reporting and ranking are delegated
to the attributed snapshot. This facade only validates input and adds provenance.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


@dataclass(frozen=True)
class IdentificationConfig:
    precursor_ppm: float = 10.0
    fragment_ppm: float | None = 20.0
    fragment_da: float | None = None
    min_total_score: float = 50.0
    min_relative_intensity: float = 0.002
    top_n: int = 3
    library_path: Path | None = None

    def __post_init__(self):
        for name in ('precursor_ppm', 'fragment_ppm', 'fragment_da'):
            value = getattr(self, name)
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(f'{name} must be positive and finite')
        if self.fragment_ppm is None and self.fragment_da is None:
            raise ValueError('Specify fragment_ppm or fragment_da')
        if self.fragment_ppm is not None and self.fragment_da is not None:
            raise ValueError('Use one MS2 tolerance unit at a time')
        if not math.isfinite(self.min_total_score) or not 0 <= self.min_total_score <= 100:
            raise ValueError('min_total_score must be between 0 and 100')
        if not math.isfinite(self.min_relative_intensity) or not 0 <= self.min_relative_intensity <= 1:
            raise ValueError('min_relative_intensity must be between 0 and 1')
        if not isinstance(self.top_n, int) or isinstance(self.top_n, bool) or self.top_n < 1:
            raise ValueError('top_n must be a positive integer')


def library_path(polarity: str = '+') -> Path:
    mode = 'positive' if polarity == '+' else 'negative' if polarity == '-' else None
    if mode is None:
        raise ValueError('Polarity must be + or -')
    result = Path(str(files('sphingolipid_toolkit').joinpath('data', 'ms2', f'sphingolipid_{mode}.sqlite')))
    if not result.is_file():
        raise FileNotFoundError('The packaged sphingolipid library is missing. Install the complete package or run git lfs pull.')
    manifest=result.parent.parent/'library_revision.json'
    stat=result.stat()
    _verify_bundled_library(str(result),stat.st_size,stat.st_mtime_ns,manifest.stat().st_mtime_ns,mode)
    return result


@lru_cache(maxsize=4)
def _verify_bundled_library(path,size,mtime,manifest_mtime,mode):
    path=Path(path)
    manifest=json.loads((path.parent.parent/'library_revision.json').read_text(encoding='utf-8'))
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):digest.update(chunk)
    if digest.hexdigest()!=manifest[mode]['index_sha256']:
        raise ValueError('Bundled sphingolipid library fingerprint mismatch; restore the complete library package.')


class SphingolipidIdentifier:
    """Reuse one bounded library reader per polarity; no LipidGate install needed."""
    def __init__(self, config: IdentificationConfig | None = None):
        self.config = config or IdentificationConfig()
        self._searchers = {}

    def _searcher(self, polarity):
        if polarity not in ('+', '-'):
            raise ValueError('Polarity must be + or -')
        if polarity not in self._searchers:
            from ._vendor.lipidgate.ms2.search import LipidMS2Searcher
            from ._vendor.lipidgate.ms2.positive_sphingolipid import CLASS_KEYS
            selected = self.config.library_path or library_path(polarity)
            self._searchers[polarity] = LipidMS2Searcher(
                selected, precursor_tolerance_ppm=self.config.precursor_ppm,
                fragment_tolerance_ppm=self.config.fragment_ppm,
                fragment_tolerance_da=self.config.fragment_da,
                min_relative_intensity=self.config.min_relative_intensity,
                min_total_score=self.config.min_total_score,
                allowed_classes=tuple(CLASS_KEYS) + ('Hex-IPC',),
            )
        return self._searchers[polarity]

    def identify(self, precursor_mz: float, peaks: Sequence[tuple[float, float]], *,
                 scan_id: str = '', rt_minutes: float = 0., polarity: str = '+',
                 precursor_charge: int | None = None, metadata: dict | None = None) -> list[dict]:
        from ._vendor.lipidgate.ms2.models import ExperimentalSpectrum, normalize_peaks
        if not math.isfinite(precursor_mz) or precursor_mz <= 0 or not math.isfinite(rt_minutes):
            raise ValueError('Precursor m/z must be positive and RT finite')
        values = [(float(mz), float(intensity)) for mz, intensity in peaks]
        if any(not math.isfinite(mz) or mz <= 0 or not math.isfinite(intensity) or intensity < 0 for mz, intensity in values):
            raise ValueError('Peak masses must be positive and intensities nonnegative and finite')
        normalized = [peak for peak in normalize_peaks(values)
                      if peak.relative_intensity >= self.config.min_relative_intensity]
        spectrum = ExperimentalSpectrum(str(scan_id), float(precursor_mz), float(rt_minutes),
            polarity, normalized, metadata=dict(metadata or {}), precursor_charge=precursor_charge)
        return self.score_spectrum(spectrum)

    def score_spectrum(self, spectrum) -> list[dict]:
        # A prepared reference spectrum is used as-is; never renormalize a
        # second time or reapply historical absolute-intensity/coverage gates.
        rows = self._searcher(spectrum.polarity).score_spectrum(spectrum, top_n=self.config.top_n)
        for row in rows:
            from .library_catalog import reference_aliases_for_class
            row['reference_class_aliases']='; '.join(reference_aliases_for_class(row['compound_class']))
            row['annotation_engine'] = 'SphinGOlipID'
            row['engine_upstream'] = 'LipidGate'
            row['engine_revision'] = 'sphingolipid-migration-20261110'
            row['ms1_tolerance_ppm'] = self.config.precursor_ppm
            row['ms2_tolerance_ppm'] = self.config.fragment_ppm
            row['ms2_tolerance_da'] = self.config.fragment_da
        return rows

    def close(self):
        for searcher in self._searchers.values():
            searcher.close()
        self._searchers.clear()

    def identify_mzml(self,path):
        from .reference_pipeline import iter_prepared_mzml
        path=Path(path)
        rows=[]
        for spectrum in iter_prepared_mzml(path,self.config.min_relative_intensity):
            identified=self.score_spectrum(spectrum)
            for row in identified:
                row['source_file']=path.name
                row['native_scan_id']=spectrum.metadata.get('native_scan_id','')
                row['collision_energy']=spectrum.metadata.get('collision_energy')
            rows.extend(identified)
        return compatible_frame(rows)

    def __enter__(self): return self
    def __exit__(self, *exc): self.close()


def compatible_frame(rows: Iterable[dict], *, abundance=None) -> pd.DataFrame:
    """Keep exact reference columns and add aliases for existing RT/export tools."""
    data = pd.DataFrame(rows)
    if data.empty:
        return pd.DataFrame(columns=['lipid_name','observed_mz','observed_rt','rank','total_score',
            'ms2_match_score','matched_fragment_count','ms2_evidence_json'])
    data['lipid_name'] = data['matched_name']
    data['注释'] = data['matched_name']
    data['observed_mz'] = data['precursor_mz']
    data['target'] = data['precursor_mz']
    data['observed_rt'] = data['rt_minutes']
    data['RT'] = data['rt_minutes']
    data['rank'] = data['result_rank']
    data['ms2_match_score'] = data['total_score'] / 100.
    data['匹配度分数'] = data['ms2_match_score']
    data['总分数'] = data['total_score']
    evidence=data['ms2_evidence_json'].map(json.loads)
    matched=evidence.map(lambda item:[f for f in item['fragments'] if f['observed_mz'] is not None])
    data['实际mz']=matched.map(lambda items:[f['observed_mz'] for f in items])
    data['强度']=matched.map(lambda items:[f['intensity'] for f in items])
    data['碎片']=matched.map(lambda items:[f['name'] for f in items])
    data['母离子']=data['library_species_name']
    data['Abund'] = abundance
    data['signal_intensity'] = abundance
    data['强度总和'] = data['matched_intensity_sum']
    return data


def config_from_workflow(config) -> IdentificationConfig:
    return IdentificationConfig(precursor_ppm=getattr(config, 'precursor_ppm', getattr(config, 'ms1_candidate_ppm',10.)),
        fragment_ppm=None if getattr(config,'fragment_da',None) is not None else config.fragment_ppm,
        fragment_da=getattr(config,'fragment_da',None),
        min_total_score=getattr(config,'min_total_score',50.),top_n=config.top_n,
        library_path=getattr(config,'ms2_library_path',None))
