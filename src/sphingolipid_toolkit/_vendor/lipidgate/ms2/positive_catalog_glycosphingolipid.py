"""Conservative composition-level coverage for historical glycosphingolipid names.

Only the glycan composition is represented. Y0/Z0 are protonated retained
ceramide products; B ions are protonated glycan oxonium products. No linkage,
hexose epimer or ceramide chain assignment is inferred from these channels.
"""
from __future__ import annotations

from functools import lru_cache

from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C, HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N, OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON, WATER_MONOISOTOPIC_MASS as WATER,
)
from .models import FragmentRecord

RESIDUE_NAMES = ('Hex', 'HexNAc', 'Fuc', 'NeuAc', 'NeuGc', 'KDN')
RESIDUE_FORMULAS = ((6,10,0,5),(8,13,1,5),(6,10,0,4),(11,17,1,8),(11,17,1,9),(9,14,0,8))
RESIDUE_MASSES = tuple(c*C+h*H+n*N+o*O for c,h,n,o in RESIDUE_FORMULAS)

# Existing HexCer/LacCer/Hex3Cer/Gb4/GM1/GM2/GM3 policies are unchanged.
# Each new family has only [M+H]+ total-composition entries.
COMPOSITIONS = (
    (2,2,0,0,0,0),(2,2,1,0,0,0),(2,3,1,0,0,0),(2,1,0,0,0,0),
    (2,1,0,0,0,1),(2,0,0,0,0,1),(2,0,0,0,1,0),(3,0,1,0,0,0),
    (3,2,0,0,0,0),(3,2,2,0,0,0),(3,2,1,0,0,0),(3,1,2,0,0,0),
    (3,1,1,0,0,0),(3,1,0,0,0,1),(4,0,0,0,0,0),(4,0,1,0,0,0),
    (4,2,0,0,0,0),(4,2,2,0,0,0),(4,2,3,0,0,0),(4,2,1,0,0,0),
    (4,3,0,0,0,0),(4,3,2,0,0,0),(4,3,1,0,0,0),(4,1,0,0,0,0),
    (4,1,2,0,0,0),(4,1,1,0,0,0),(5,0,0,0,0,0),(5,2,1,0,0,0),
    (5,3,0,0,0,0),(5,3,1,0,0,0),(5,1,1,0,0,0),(1,0,0,1,0,0),
)


def class_name(composition):
    return ''.join(kind+(str(n) if n != 1 else '')
                   for kind,n in zip(RESIDUE_NAMES,composition) if n)+'Cer'


PROFILES = {class_name(composition): composition for composition in COMPOSITIONS}
CLASSES = frozenset(PROFILES)


def is_positive_catalog_glycosphingolipid(compound_class, adduct):
    return compound_class in CLASSES and adduct == '[M+H]+'


@lru_cache(maxsize=64)
def signature_names(compound_class):
    composition = PROFILES[compound_class]
    for kind,n in zip(RESIDUE_NAMES[3:],composition[3:]):
        if n:
            return frozenset({'B-'+kind,'B-'+kind+'-H2O'})
    return frozenset({'B-glycan','B-glycan-H2O'})


def fragments(compound_class, precursor_mz):
    composition = PROFILES[compound_class]
    glycan = sum(n*mass for n,mass in zip(composition,RESIDUE_MASSES))
    core = float(precursor_mz)-glycan
    values = [
        FragmentRecord(round(core,4),'Y0','Diagnostic_HG'),
        FragmentRecord(round(core-WATER,4),'Z0','Diagnostic_HG'),
        FragmentRecord(round(glycan+PROTON,4),'B-glycan','Diagnostic_B'),
        FragmentRecord(round(glycan+PROTON-WATER,4),'B-glycan-H2O','Diagnostic_B'),
    ]
    for i in (3,4,5):
        if composition[i]:
            name = 'B-'+RESIDUE_NAMES[i]
            values.extend((FragmentRecord(round(RESIDUE_MASSES[i]+PROTON,4),name,'Diagnostic_B'),
                FragmentRecord(round(RESIDUE_MASSES[i]+PROTON-WATER,4),name+'-H2O','Diagnostic_B')))
    values.extend((FragmentRecord(round(precursor_mz-WATER,4),'M+H-H2O','Common'),
                   FragmentRecord(round(precursor_mz,4),'[M+H]+','Precursor Ion')))
    # Retain one occurrence per physical product, never two gate votes for
    # alternate names of a product with the same mass and role.
    unique = {}
    for fragment in values:
        unique.setdefault((fragment.mz,fragment.fragment_type),fragment)
    return list(unique.values())
