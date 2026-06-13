"""Formula-library module placeholder.

Planned responsibility:
- generate the theoretical sphingolipid formula library from class templates,
  LCB/FA chain ranges, unsaturation, hydroxylation, headgroups, and glycan
  composition rules;
- export the library to the MS1 database format required by ms2_pipeline.py.

This file is intentionally lightweight now so later formula-enumeration code can
be added without changing the MS2 matching API.
"""

from __future__ import annotations


def build_formula_library(*args, **kwargs):
    raise NotImplementedError("Formula-library construction will be added in the next module integration step.")
