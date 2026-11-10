"""Historical MS1 reference names, separate from MS2 identification evidence."""
import csv
from functools import lru_cache
import gzip
from importlib.resources import files
import json
import math


def iter_reference_names():
    source=files('sphingolipid_toolkit').joinpath('data','reference_names.csv.gz')
    with source.open('rb') as binary,gzip.open(binary,'rt',encoding='utf-8',newline='') as handle:
        yield from csv.DictReader(handle)


def reference_name_candidates(precursor_mz,ppm=10.,*,include_reference_only=False):
    """Return MS1-only candidates with explicit MS2-coverage status.

    This function never admits, confirms chains or ranks an MS2 identification.
    """
    if not math.isfinite(precursor_mz) or precursor_mz<=0 or not math.isfinite(ppm) or ppm<=0:
        raise ValueError('Precursor m/z and ppm must be positive and finite')
    result=[]
    for row in iter_reference_names():
        if not row['precursor_mz']:continue
        if not include_reference_only and row['coverage_status']!='ms2_composition_covered':continue
        value=float(row['precursor_mz'])
        error=(precursor_mz-value)/value*1e6
        if abs(error)<=ppm:
            result.append(dict(row,precursor_error_ppm=error,evidence_level='MS1 reference candidate'))
    return result


@lru_cache(maxsize=1)
def _aliases():
    path=files('sphingolipid_toolkit').joinpath('data','name_aliases.json')
    try:return json.loads(path.read_text(encoding='utf-8'))
    except FileNotFoundError:return {}


def reference_aliases_for_class(compound_class):
    return tuple(_aliases().get(compound_class,()))
