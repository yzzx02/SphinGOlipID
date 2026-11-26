"""Positive glycosphingolipid headgroup channels, by precursor and product charge."""
from dataclasses import replace
from functools import lru_cache

from .library_fragment_policy import WATER_MONOISOTOPIC_MASS

GM2_B = frozenset({'B1alpha', 'B1alpha-H2O', 'B2Y2alpha', 'B2', 'B3'})
GM2_CORE = GM2_B | {'Y0', 'Z0', 'Y2beta'}
GANGLIOSIDES = frozenset({'GM1', 'GD1', 'GD3', 'GT1', 'Gal-Fuc-GM1', 'Gal-Fuc-GD1'})
CURATED = GANGLIOSIDES | {'GM2', 'GM3', 'Gb4'}
REVISION = 'sphingolipid-hg8-20261124-r1'


@lru_cache(maxsize=32)
def yz_core(cls, double_charge):
    n = {'GM1': 1, 'GD1': 2, 'GD3': 2, 'GT1': 3, 'Gal-Fuc-GM1': 1, 'Gal-Fuc-GD1': 2}[cls]
    if cls.startswith('Gal-Fuc'):
        ending = '^2+' if double_charge else ''
        products = (f'[Hex3;HexNAc1;NeuAc{n}]', f'[Hex4;HexNAc1;NeuAc{n}]'+ending,
                    f'[Hex3;HexNAc1;NeuAc{n};Fuc1]'+ending)
    elif cls == 'GD3' and double_charge:
        products = ('[Hex2]', '[Hex2;NeuAc1]', '[Hex2;NeuAc1]^2+')
    else:
        base = 'Hex2' if cls == 'GD3' else 'Hex3;HexNAc1'
        residual = '[' + base + (f';NeuAc{n-1}' if n>1 else '') + ']'
        products = ('1', residual, residual+'^2+' if double_charge and n>1 else '[Hex2]')
    return frozenset({'Y0','Z0'} | {p+s for p in ('Y','Z') for s in products})


def curate_headgroup(cls, adduct, fragments):
    if cls not in CURATED or not adduct.endswith('+'):
        return fragments
    result = []
    seen = set()
    for f in fragments:
        name = f.name.removeprefix('Aux ')
        if cls == 'GM2':
            auxiliary = (f.fragment_type in {'Diagnostic_HG','Diagnostic_B'} or name.startswith(('Y','Z','B'))) and name not in GM2_CORE
        elif cls in GANGLIOSIDES:
            auxiliary = name.startswith(('Y','Z')) and name != 'Z0-H2O' and name not in yz_core(cls, adduct.endswith('2+'))
        elif cls == 'Gb4':
            auxiliary = name == 'M+H-HexNAc-3Hex-2H2O'
        else:
            auxiliary = name in {'[Neu5Ac+Hex+H-H2O]+','[Neu5Ac+2Hex+H-H2O]+'}
        if auxiliary:
            f = replace(f, name='Aux '+name, fragment_type='Common', required_group=None)
        elif cls=='GM2' and name in GM2_B and f.fragment_type!='Diagnostic_B':
            f = replace(f, fragment_type='Diagnostic_B')
        # Imported older family MSPs can contain a regenerated counterpart.
        key = (f.name,round(f.mz,4))
        if key not in seen:
            result.append(f)
            seen.add(key)
    if cls=='GM2' and not any(f.name=='B1alpha-H2O' for f in result):
        b1=next((f for f in result if f.name=='B1alpha'),None)
        if b1:
            result.append(replace(b1,name='B1alpha-H2O',mz=round(b1.mz-WATER_MONOISOTOPIC_MASS,4)))
    return sorted(result,key=lambda f:f.mz)
