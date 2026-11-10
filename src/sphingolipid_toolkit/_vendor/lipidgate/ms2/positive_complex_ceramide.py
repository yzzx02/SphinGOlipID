"""Positive IPC and acylhexosylceramide products, calculated from explicit chains."""
from __future__ import annotations

from dataclasses import dataclass,replace
import re

from .models import FragmentRecord
from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,WATER_MONOISOTOPIC_MASS as WATER,
)

P=30.97376199842
HEX=6*C+10*H+5*O
IP_RESIDUE=6*C+11*H+8*O+P
IP_ACID=IP_RESIDUE+WATER
IPC_GATE='M+H-C6H13O9P'
AHEX_CORE=frozenset({'Cer+H','Cer+H-H2O','Cer+H-2H2O'})
AHEX_ACYL=frozenset({'AcylHex+H-H2O','M+H-RCOOH','M+H-(RCOOH-H2O)'})
LCB_NAMES=frozenset({'LCB-H2O','LCB-2H2O'})
LCB_BASE=re.compile(r'(?:\(|\s)([dt])(\d+):(\d+)(?:/|\(O-)',re.IGNORECASE)
AHEX_NAME=re.compile(r'^AHexCer\s+([dt])(\d+):(\d+)\(O-(\d+):(\d+)\)/(\d+):(\d+)(\(OH\))?$',re.IGNORECASE)
IPC_NAME=re.compile(r'^PI-Cer\(([dt])(\d+):(\d+)/(h?)(\d+):(\d+)\)$',re.IGNORECASE)


def is_positive_complex_ceramide(compound_class,adduct):
    return (compound_class=='PI-Cer' and adduct=='[M+H]+'
            or compound_class=='AHexCer' and adduct in {'[M+H]+','[M+H-H2O]+'})


@dataclass(frozen=True,order=True)
class CeramideCore:
    series:str
    lcb_carbon:int
    lcb_db:int
    hydroxy_acyl:bool
    acyl_carbon:int
    acyl_db:int

    @property
    def counts(self):
        carbon=self.lcb_carbon+self.acyl_carbon
        return carbon,2*carbon+1-2*(self.lcb_db+self.acyl_db),1,(3 if self.series=='t' else 2)+1+int(self.hydroxy_acyl)

    @property
    def mass(self):return sum(n*m for n,m in zip(self.counts,(C,H,N,O)))

    @property
    def lcb_mass(self):
        return self.lcb_carbon*C+(2*self.lcb_carbon+3-2*self.lcb_db)*H+N+(3 if self.series=='t' else 2)*O+PROTON

    @property
    def lcb_ions(self):
        return [FragmentRecord(round(self.lcb_mass-loss*WATER,4),name,'LCB碎片') for name,loss in [('LCB-H2O',1),('LCB-2H2O',2)]]


@dataclass(frozen=True,order=True)
class ComplexCeramide:
    compound_class:str
    core:CeramideCore
    outer_carbon:int=0
    outer_db:int=0

    @property
    def name(self):
        b=self.core
        if self.compound_class=='PI-Cer':return f'PI-Cer({b.series}{b.lcb_carbon}:{b.lcb_db}/{"h" if b.hydroxy_acyl else ""}{b.acyl_carbon}:{b.acyl_db})'
        return f'AHexCer {b.series}{b.lcb_carbon}:{b.lcb_db}(O-{self.outer_carbon}:{self.outer_db})/{b.acyl_carbon}:{b.acyl_db}{"(OH)" if b.hydroxy_acyl else ""}'

    @property
    def species_name(self):
        b=self.core;oxygen=(3 if b.series=='t' else 2)+int(b.hydroxy_acyl)
        return f'{self.compound_class} {b.lcb_carbon+b.acyl_carbon+self.outer_carbon}:{b.lcb_db+b.acyl_db+self.outer_db};O{oxygen}'

    @property
    def formula(self):
        c,h,n,o=self.core.counts
        if self.compound_class=='PI-Cer':return f'C{c+6}H{h+11}NO{o+8}P'
        return f'C{c+6+self.outer_carbon}H{h+10+2*self.outer_carbon-2*self.outer_db-2}NO{o+6}'

    @property
    def outer_fa(self):return self.outer_carbon*C+(2*self.outer_carbon-2*self.outer_db)*H+2*O

    def precursor_mz(self,adduct):
        mass=(self.core.mass+IP_RESIDUE if self.compound_class=='PI-Cer' else self.core.mass+HEX+self.outer_fa-WATER)
        return mass+PROTON-(WATER if adduct=='[M+H-H2O]+' else 0.)

    def fragments(self,adduct):
        precursor=self.precursor_mz(adduct);core=self.core.mass+PROTON
        def f(mz,name,kind='Common'):return FragmentRecord(round(mz,4),name,kind)
        if self.compound_class=='PI-Cer':
            # Loss of neutral inositol phosphate yields dehydrated Cer, not intact Cer.
            peaks=[f(precursor-IP_ACID,IPC_GATE,'Diagnostic_HG'),
                f(precursor-IP_ACID-WATER,IPC_GATE+'-H2O'),
                f(precursor-WATER,'M+H-H2O')]
        else:
            peaks=[f(self.outer_fa+HEX-WATER+PROTON,'AcylHex+H-H2O','Diagnostic_HG'),
                f(precursor-self.outer_fa,'M+H-RCOOH','Diagnostic_HG'),
                f(precursor-self.outer_fa+WATER,'M+H-(RCOOH-H2O)','Diagnostic_HG'),
                f(core-WATER,'Cer+H-H2O','Diagnostic_HG'),
                f(core-2*WATER,'Cer+H-2H2O','Diagnostic_HG')]
            if adduct=='[M+H]+':peaks.extend([f(core,'Cer+H','Diagnostic_HG'),f(precursor-WATER,'M+H-H2O')])
        peaks.extend(self.core.lcb_ions)
        if self.compound_class=='AHexCer' and self.core.series=='d':
            lcb=self.core.lcb_mass
            # Retain the legacy d-series support channels as optional Other.
            # These do not establish a chain or add an admission requirement.
            peaks.extend([f(lcb-(C+4*H+2*O),'LCB-CH4O2'),
                          f(lcb-(2*C+7*H+N+O),'LCB-C2H7NO')])
        peaks.append(f(precursor,adduct,'Precursor Ion'))
        return sorted(peaks,key=lambda p:p.mz)


def from_name(name):
    match=AHEX_NAME.fullmatch(name)
    if match:
        series,lc,ld,oc,od,ac,ad,oh=match.groups()
        return ComplexCeramide('AHexCer',CeramideCore(series.lower(),int(lc),int(ld),bool(oh),int(ac),int(ad)),int(oc),int(od))
    match=IPC_NAME.fullmatch(name)
    if match:
        series,lc,ld,oh,ac,ad=match.groups()
        return ComplexCeramide('PI-Cer',CeramideCore(series.lower(),int(lc),int(ld),bool(oh),int(ac),int(ad)))
    return None


def species_name(record):
    model=from_name(record.lipid_chain_name)
    return model.species_name if model else record.lipid_name


def normalize_fragments(compound_class,adduct,chain_name,fragments):
    """Migrate known legacy AHexCer labels without inventing observed products."""
    if compound_class!='AHexCer' or not is_positive_complex_ceramide(compound_class,adduct):return fragments
    model=from_name(chain_name)
    if model is None:return fragments
    result=[]
    for fragment in fragments:
        name=fragment.name
        if name.startswith('M+H-Acyl(') and '-C6H10O5' in name:
            name='Cer+H'+('-2H2O' if name.endswith('-2H2O') else '-H2O' if name.endswith('-H2O') else '')
        elif name.startswith('O-') and name.endswith('-Hex+'):name='AcylHex+H-H2O'
        kind=('Diagnostic_HG' if name in AHEX_CORE|AHEX_ACYL else 'LCB碎片' if name in LCB_NAMES else
              'Precursor Ion' if fragment.fragment_type=='Precursor Ion' else 'Common')
        result.append(replace(fragment,name=name,fragment_type=kind,required_group='hg' if kind=='Diagnostic_HG' else None))
    return result
