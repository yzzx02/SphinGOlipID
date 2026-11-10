"""Retain well-defined legacy subtype families at composition resolution.

FMC aliases describe acetylation/hydroxylation families rather than unique
modification sites. All masses are elemental, including acetic-acid loss.
"""
from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,HYDROGEN_MONOISOTOPIC_MASS as H,
    OXYGEN_MONOISOTOPIC_MASS as O,PHOSPHORUS_MONOISOTOPIC_MASS as P,
    NITROGEN_MONOISOTOPIC_MASS as N,PROTON_MONOISOTOPIC_MASS as PROTON,
    WATER_MONOISOTOPIC_MASS as WATER,
)
from .models import FragmentRecord

S=31.9720711744
HEX=6*C+10*H+5*O
ACETYL=2*C+2*H+O
ACETIC_ACID=2*C+4*H+2*O
INOSITOL_PHOSPHATE=6*C+13*H+9*O+P
SO3=S+3*O

# Attached residues in CHNOPS order, excluding the ceramide/base backbone.
HEADS={
    'Ac1HexCer':(8,12,0,6,0,0),'Ac2HexCer':(10,14,0,7,0,0),
    'Ac5HexCer':(16,20,0,10,0,0),'Ac6HexCer':(18,22,0,11,0,0),
    'MIP2C':(18,32,0,21,2,0),
    'HexSPB':(6,10,0,5,0,0),'Hex3SPB':(18,30,0,15,0,0),'SHexSPB':(6,10,0,8,0,1),
    'KetoSPB':(0,-2,0,0,0,0),
    'N,N-dimethylSo':(2,4,0,0,0,0),
}
ALIASES={'FMC-1':'Ac1HexCer','FMC-2':'Ac1HexCer','FMC-3':'Ac2HexCer','FMC-4':'Ac2HexCer',
    'FMC-5':'Ac5HexCer','FMC-6':'Ac5HexCer','FMC-7':'Ac6HexCer','M(IP)2C':'MIP2C',
    'GlcSo':'HexSPB','GalSo':'HexSPB','Gb3So':'Hex3SPB','Lyso-sulfo':'SHexSPB'}
ALIASES['ketoSa']='KetoSPB'
ALIASES['N,N-dimethylSo']='N,N-dimethylSo'
CLASSES=frozenset(HEADS)
FREE_BASE_CLASSES=frozenset({'HexSPB','Hex3SPB','SHexSPB','KetoSPB','N,N-dimethylSo'})
ACETYL_CLASSES=frozenset({'Ac1HexCer','Ac2HexCer','Ac5HexCer','Ac6HexCer'})


def is_positive_catalog_legacy_sphingolipid(compound_class,adduct):
    return compound_class in CLASSES and adduct=='[M+H]+'


def required_names(compound_class):
    if compound_class=='N,N-dimethylSo':return {'DMS-C7H12N+'}
    if compound_class=='KetoSPB':return {'M+H-CH2O'}
    if compound_class in ACETYL_CLASSES:return {'M+H-HOAc','M+H-glycan-HOAc'}
    return {'Y0','Z0'}


def alternative_names(compound_class):
    if compound_class in ACETYL_CLASSES:return {'M+H-glycan','M+H-glycan-H2O'}
    if compound_class=='MIP2C':return {'M+H-InositolP','M+H-InositolP-H2O'}
    return set()


def hg_minimum(compound_class):
    if compound_class=='KetoSPB':return 1
    return 3 if compound_class in ACETYL_CLASSES or compound_class=='MIP2C' else 2


def fragments(compound_class,precursor_mz,core_oxygens=3):
    products=[]
    def add(name,mass,kind='Diagnostic_HG'):
        if mass>0:products.append(FragmentRecord(round(mass,4),name,kind))
    if compound_class=='N,N-dimethylSo':
        # Zhu et al., Nat Protoc 2013, Fig. 7b: native positive CID, 20/40 V.
        # Exact formulas below are elemental assignments to published 3-decimal
        # masses, not claimed published exact formula annotations. Only the
        # d18:1 reference is curated; never extrapolate a methylated-base grid.
        for name,carbon,hydrogen in [('DMS-C3H8N+',3,8),('DMS-C5H10N+',5,10),
                                     ('DMS-C7H12N+',7,12)]:
            add(name,carbon*C+(hydrogen-1)*H+N+PROTON)
        add('DMS-C6H10N+',6*C+9*H+N+PROTON,'Common')
        add('M+H-2H2O',precursor_mz-2*WATER,'Common')
        add('M+H-CH2O-H2O',precursor_mz-(C+2*H+O)-WATER,'Common')
    elif compound_class=='KetoSPB':
        formaldehyde=C+2*H+O
        add('M+H-CH2O',precursor_mz-formaldehyde)
        add('M+H-CH2O-H2O',precursor_mz-formaldehyde-WATER,'Common')
    elif compound_class in ACETYL_CLASSES:
        count={'Ac1HexCer':0,'Ac2HexCer':1,'Ac5HexCer':4,'Ac6HexCer':4}[compound_class]
        loss=HEX+count*ACETYL
        add('M+H-glycan',precursor_mz-loss)
        add('M+H-glycan-H2O',precursor_mz-loss-WATER)
        add('M+H-HOAc',precursor_mz-ACETIC_ACID)
        add('M+H-glycan-HOAc',precursor_mz-loss-ACETIC_ACID)
        add('M+H-glycan-H2O-HOAc',precursor_mz-loss-WATER-ACETIC_ACID)
    elif compound_class=='MIP2C':
        head=18*C+32*H+21*O+2*P
        add('Y0',precursor_mz-head)
        add('Z0',precursor_mz-head-WATER)
        add('M+H-InositolP',precursor_mz-INOSITOL_PHOSPHATE)
        add('M+H-InositolP-H2O',precursor_mz-INOSITOL_PHOSPHATE-WATER)
        add('M+H-Hex',precursor_mz-HEX,'Common')
    elif compound_class=='SHexSPB':
        head=HEX+SO3
        add('Y0',precursor_mz-head)
        add('Z0',precursor_mz-SO3)
        add('LCB-H2O',precursor_mz-head-WATER,'Common')
        if core_oxygens>=2:add('LCB-2H2O',precursor_mz-head-2*WATER,'Common')
    else:
        head=HEX*(3 if compound_class=='Hex3SPB' else 1)
        add('Y0',precursor_mz-head)
        add('Z0',precursor_mz-head-WATER)
        if core_oxygens>=2:add('LCB-2H2O',precursor_mz-head-2*WATER,'Common')
        if compound_class=='Hex3SPB':
            add('M+H-Hex',precursor_mz-HEX,'Common')
            add('M+H-2Hex',precursor_mz-2*HEX,'Common')
    add('M+H-H2O',precursor_mz-WATER,'Common')
    add('[M+H]+',precursor_mz,'Precursor Ion')
    unique={}
    for product in products:unique.setdefault((product.mz,product.fragment_type),product)
    return list(unique.values())
