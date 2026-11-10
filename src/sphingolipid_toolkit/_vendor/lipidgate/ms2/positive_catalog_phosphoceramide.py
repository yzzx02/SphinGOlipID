"""Composition-level positive phosphoceramide coverage retained from the catalog."""
from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C, HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N, OXYGEN_MONOISOTOPIC_MASS as O,
    PHOSPHORUS_MONOISOTOPIC_MASS as P, PROTON_MONOISOTOPIC_MASS as PROTON,
    WATER_MONOISOTOPIC_MASS as WATER,
)
from .models import FragmentRecord

# Neutral head formulas (the lipid contains head minus water).
HEADS = {'CAEP': (2,8,1,3,1), 'N-methylCAEP': (3,10,1,3,1), 'PG-Cer': (3,9,0,6,1)}
ALIASES = {'CAEP':'CAEP', 'N-methylCAEP':'N-methylCAEP', 'PGCer':'PG-Cer'}
CLASSES = frozenset(HEADS)


def is_positive_catalog_phosphoceramide(compound_class, adduct):
    return compound_class in CLASSES and adduct == '[M+H]+'


def gate_names(compound_class):
    return {'M+H-head'} if compound_class == 'PG-Cer' else {'M+H-head','HG-phosphonate'}


def fragments(compound_class, precursor_mz):
    head = sum(n*mass for n,mass in zip(HEADS[compound_class],(C,H,N,O,P)))
    values = [FragmentRecord(round(precursor_mz-head,4),'M+H-head','Diagnostic_HG')]
    if compound_class != 'PG-Cer':
        values.append(FragmentRecord(round(head+PROTON,4),'HG-phosphonate','Diagnostic_HG'))
    values.extend((FragmentRecord(round(precursor_mz-head-WATER,4),'M+H-head-H2O','Common'),
                   FragmentRecord(round(precursor_mz-WATER,4),'M+H-H2O','Common'),
                   FragmentRecord(round(precursor_mz,4),'[M+H]+','Precursor Ion')))
    return values
