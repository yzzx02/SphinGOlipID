"""Composition-based positive ganglioside family spectra, without linkage claims.

Products follow glycosidic cleavages of the conventional family glycan trees.
Single and internal two-cut products keep their glycan composition; changing
the ceramide shifts only the Y/Z products. Ammoniated parents lose ammonia,
while glycan and ceramide products are protonated rather than ammoniated.
"""
from __future__ import annotations

from dataclasses import dataclass,replace
from functools import lru_cache
import re

from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,WATER_MONOISOTOPIC_MASS as WATER,
)
from .models import FragmentRecord

LEGACY_CLASSES=frozenset({'GM1','GD1','GT1','GD3'})
MODIFIED_CLASSES=frozenset({'Gal-Fuc-GM1','Gal-Fuc-GD1'})
CLASSES=LEGACY_CLASSES | MODIFIED_CLASSES
NH3=N+3*H
# charge, ammonium count, pre-existing neutral water loss
ADDUCTS={
    '[M+H]+':(1,0,0),'[M+NH4]+':(1,1,0),
    '[M+2H]2+':(2,0,0),'[M+H+NH4]2+':(2,1,0),'[M+2NH4]2+':(2,2,0),
    '[M+H-H2O]+':(1,0,1),'[M+2H-H2O]2+':(2,0,1),
}
HEX=(6,10,0,5);HEXNAC=(8,13,1,5);NEUAC=(11,17,1,8);FUC=(6,10,0,4)
RESIDUES=(HEX,HEXNAC,NEUAC,FUC)
RESIDUE_MASSES=tuple(c*C+h*H+n*N+o*O for c,h,n,o in RESIDUES)
# Each entry gives residue types (Hex=0, HexNAc=1, NeuAc=2) and tree edges.
# The root is the glucose attached to Cer. Family profiles pool linkage isomers.
TREES={
 'GM1':(
   ((0,0,2,1,0),((0,1),(1,2),(1,3),(3,4))),
   ((0,0,1,0,2),((0,1),(1,2),(2,3),(3,4)))),
 'GD1':(
   ((0,0,2,1,0,2),((0,1),(1,2),(1,3),(3,4),(4,5))),
   ((0,0,2,1,0,2),((0,1),(1,2),(1,3),(3,4),(2,5)))),
 'GT1':(
   ((0,0,2,1,0,2,2),((0,1),(1,2),(1,3),(3,4),(4,5),(5,6))),
   ((0,0,2,1,0,2,2),((0,1),(1,2),(1,3),(3,4),(2,5),(4,6))),
   ((0,0,2,1,0,2,2),((0,1),(1,2),(1,3),(3,4),(2,5),(5,6)))),
 'GD3':(((0,0,2,2),((0,1),(1,2),(2,3))),),
}
# Append a Hex/Fuc cap to the existing conventional glycan trees. This pools
# linkage alternatives; a family name never establishes the positions of the cap.
for _modified,_base in (('Gal-Fuc-GM1','GM1'),('Gal-Fuc-GD1','GD1')):
    _profiles=[]
    for _residues,_edges in TREES[_base]:
        _terminal_hex=max(i for i,kind in enumerate(_residues) if kind==0)
        _n=len(_residues)
        _profiles.append((_residues+(0,3),_edges+((_terminal_hex,_n),(_terminal_hex,_n+1))))
    TREES[_modified]=tuple(_profiles)
CHAIN=re.compile(r'^(GM1|GD1|GT1|GD3|Gal-Fuc-GM1|Gal-Fuc-GD1)\(([dt])(\d+):(\d+)/(h?)(\d+):(\d+)\)$')
NEU_NAMES=frozenset({'B-NeuAc','B-NeuAc-H2O'})
REFERENCE_B={'Reference-HG-366.1400':'B-HexHexNAc','Reference-HG-657.2354':'B[Hex1;HexNAc1;NeuAc1]'}

def is_positive_ganglioside(compound_class,adduct):
    return compound_class in CLASSES and adduct in ADDUCTS

def _mass(composition):
    return sum(n*m for n,m in zip(composition,RESIDUE_MASSES))

def _label(composition):
    return ';'.join(name+str(n) for name,n in zip(('Hex','HexNAc','NeuAc','Fuc'),composition) if n)

def _dimension(compound_class):
    return 4 if compound_class in MODIFIED_CLASSES else 3

def _b_name(composition):
    if not any(composition[3:]):
        if composition[:3]==(0,0,1):return 'B-NeuAc'
        if composition[:3]==(1,1,0):return 'B-HexHexNAc'
    return 'B['+_label(composition)+']'

def _connected(nodes,edges):
    reached={next(iter(nodes))}
    while True:
        previous=len(reached)
        for a,b in edges:
            if a in nodes and b in nodes and (a in reached or b in reached):reached.update((a,b))
        if len(reached)==previous:return reached==nodes

@lru_cache(maxsize=6)
def glycan_compositions(compound_class):
    dimension=_dimension(compound_class)
    backbones={(0,)*dimension};glycans=set()
    for residues,edges in TREES[compound_class]:
        full=(1<<len(residues))-1
        for mask in range(1,full+1):
            nodes={i for i in range(len(residues)) if mask>>i&1}
            if not _connected(nodes,edges):continue
            counts=tuple(sum(residues[i]==kind for i in nodes) for kind in range(dimension))
            cuts=sum((a in nodes)!=(b in nodes) for a,b in edges)
            if 0 in nodes and mask!=full and cuts<=2:backbones.add(counts)
            if cuts+int(0 in nodes)<=2:glycans.add(counts)
    if compound_class in MODIFIED_CLASSES:
        # The labile Hex/Fuc cap may be lost before further core fragmentation.
        # Retain the established core products as sequential cleavage products;
        # these shared peaks still cannot satisfy the separate cap-loss gate.
        core=compound_class.removeprefix('Gal-Fuc-')
        core_backbones,core_glycans=glycan_compositions(core)
        backbones.update(c+(0,) for c in core_backbones)
        glycans.update(c+(0,) for c in core_glycans)
    return tuple(sorted(backbones)),tuple(sorted(glycans))

@lru_cache(maxsize=12)
def fragment_specs(compound_class,charge):
    """Return names, Cer-relative/absolute mass, product charge and pool type."""
    backbones,glycans=_selected_compositions(compound_class)
    full=tuple(TREES[compound_class][0][0].count(i) for i in range(_dimension(compound_class)))
    specs=[]
    for composition in backbones:
        suffix='0' if not any(composition) else '1' if composition[0]==1 and not any(composition[1:]) else '['+_label(composition)+']'
        retains_two_charges=(charge==2 and composition[:2]==full[:2] and 0<composition[2]<full[2])
        if compound_class in MODIFIED_CLASSES:
            # Cap loss can retain both protons even with all NeuAc residues.
            retains_two_charges=charge==2 and composition[2]>0
        for product_charge in ((1,2) if retains_two_charges else (1,)):
            ending='^2+' if product_charge==2 else ''
            specs.extend((prefix+suffix+ending,_mass(composition)-loss,product_charge,'Diagnostic_HG',True)
                         for prefix,loss in (('Y',0.),('Z',WATER)))
    for composition in glycans:
        name=_b_name(composition)
        # Small HexNAc ions are useful support but cannot admit a ganglioside.
        if composition[:3]==(0,1,0) and not any(composition[3:]):
            specs.extend((label,_mass(composition)+PROTON-loss,1,'Common',False)
                         for label,loss in (('HexNAc',0.),('HexNAc-H2O',WATER)))
        elif composition[2] or composition[1] and sum(composition)>=2 or any(composition[3:]) and sum(composition)>=2:
            specs.extend((name+ending,_mass(composition)+PROTON-loss,1,'Diagnostic_B',False)
                         for ending,loss in (('',0.),('-H2O',WATER)))
    return tuple(specs)

@lru_cache(maxsize=6)
def _selected_compositions(compound_class):
    """Keep the canonical Cer/lactosyl and desialylation ladders, not every cut.

    The tree is still used to reject disconnected or unsupported compositions.
    HexNAc-Hex and sialyl/disialyl/trisialyl glycan products supply independent
    head evidence. Low-mass HexNAc monomers remain optional Other support.
    """
    backbones,glycans=glycan_compositions(compound_class)
    full=tuple(TREES[compound_class][0][0].count(i) for i in range(_dimension(compound_class)))
    yz={(0,0,0),(1,0,0),(2,0,0),(2,1,0),(3,1,0)}
    yz.update((full[0],full[1],n) for n in range(full[2]))
    yz.update((2,0,n) for n in range(1,full[2]+1))
    b={(0,0,1),(0,1,0),(1,1,0),(1,0,1),(1,1,1),(2,1,1),(0,0,2),(0,0,3),
       full,(full[0]-1,full[1],full[2])}
    if compound_class in MODIFIED_CLASSES:
        yz={c+(0,) for c in yz}
        yz.update((hexose,1,neuac,fucose)
                  for hexose in (3,4) for neuac in range(full[2]+1) for fucose in (0,1))
        b={c+(0,) for c in b if len(c)==3} | {full,(full[0]-1,full[1],full[2],full[3])}
        b.update({(1,0,0,1),(2,0,0,1),(1,1,0,1),(2,1,0,1)})
    return tuple(c for c in backbones if c in yz),tuple(c for c in glycans if c in b)

@lru_cache(maxsize=6)
def extended_glycan_names(compound_class):
    _,glycans=_selected_compositions(compound_class)
    names=set()
    for composition in glycans:
        if not composition[1] or sum(composition)<2:continue
        name=_b_name(composition)
        names.update((name,name+'-H2O'))
    names.update(REFERENCE_B)
    return frozenset(names)

@lru_cache(maxsize=2)
def modified_core_hg_names(compound_class):
    """Require loss of the Hex+Fuc cap; shared NeuAc/Cer peaks alone cannot admit."""
    if compound_class not in MODIFIED_CLASSES:return frozenset()
    neuac=1 if compound_class=='Gal-Fuc-GM1' else 2
    suffix='['+_label((3,1,neuac,0))+']'
    return frozenset(prefix+suffix+ending for prefix in ('Y','Z') for ending in ('','^2+'))

@lru_cache(maxsize=2)
def modified_fucose_evidence_names(compound_class):
    """Distinguish the cap from near-isobaric ammoniated sialylated parents.

    Shared Cer/NeuAc and the underlying core can also occur after loss of
    NeuAc from another family. Require a Fuc-containing HG/B composition or
    a separate, single-Fuc-loss Y/Z product; Other cannot satisfy this group.
    """
    if compound_class not in MODIFIED_CLASSES:return frozenset()
    neuac=1 if compound_class=='Gal-Fuc-GM1' else 2
    lost_fuc='['+_label((4,1,neuac,0))+']'
    single_loss={prefix+lost_fuc+ending for prefix in ('Y','Z') for ending in ('','^2+')}
    return frozenset(name for charge in (1,2) for name,_,_,kind,_ in fragment_specs(compound_class,charge)
                     if kind in {'Diagnostic_HG','Diagnostic_B'} and ('Fuc1' in name or name in single_loss))

@dataclass(frozen=True,order=True)
class GangliosideComposition:
    compound_class:str
    series:str
    lcb_carbon:int
    lcb_double_bonds:int
    hydroxy_acyl:bool
    acyl_carbon:int
    acyl_double_bonds:int

    @property
    def name(self):
        return f'{self.compound_class}({self.series}{self.lcb_carbon}:{self.lcb_double_bonds}/{"h" if self.hydroxy_acyl else ""}{self.acyl_carbon}:{self.acyl_double_bonds})'

    @property
    def species_name(self):
        oxygen=(2 if self.series=='d' else 3)+self.hydroxy_acyl
        prefix={2:'d',3:'t',4:'t'}[oxygen]
        extra='(OH)' if oxygen==4 else ''
        return f'{self.compound_class}({prefix}{self.lcb_carbon+self.acyl_carbon}:{self.lcb_double_bonds+self.acyl_double_bonds}){extra}'

    @property
    def ceramide_counts(self):
        carbon=self.lcb_carbon+self.acyl_carbon
        return (carbon,2*carbon+1-2*(self.lcb_double_bonds+self.acyl_double_bonds),1,
                (2 if self.series=='d' else 3)+1+self.hydroxy_acyl)

    @property
    def ceramide_mass(self):
        return sum(n*m for n,m in zip(self.ceramide_counts,(C,H,N,O)))

    @property
    def formula_counts(self):
        counts=list(self.ceramide_counts)
        residues=TREES[self.compound_class][0][0]
        for kind in residues:
            for i,n in enumerate(RESIDUES[kind]):counts[i]+=n
        return tuple(counts)

    @property
    def formula(self):
        return ''.join(atom+str(n) for atom,n in zip(('C','H','N','O'),self.formula_counts))

    @property
    def neutral_mass(self):
        return sum(n*m for n,m in zip(self.formula_counts,(C,H,N,O)))

    def precursor_mz(self,adduct):
        charge,ammonium,waters=ADDUCTS[adduct]
        return (self.neutral_mass+charge*PROTON+ammonium*NH3-waters*WATER)/charge

    def fragments(self,adduct):
        charge,ammonium,_=ADDUCTS[adduct]
        result=[]
        for name,mass,product_charge,kind,contains_cer in fragment_specs(self.compound_class,charge):
            mz=(self.ceramide_mass+mass+product_charge*PROTON)/product_charge if contains_cer else mass
            result.append(FragmentRecord(round(mz,4),name,kind))
        lcb_mass=self.lcb_carbon*C+(2*self.lcb_carbon+3-2*self.lcb_double_bonds)*H+N+(2 if self.series=='d' else 3)*O+PROTON
        result.extend(FragmentRecord(round(lcb_mass-n*WATER,4),name,'Common')
                      for name,n in (('LCB-H2O',1),('LCB-2H2O',2)))
        parent=self.precursor_mz(adduct)
        result.extend((FragmentRecord(round(parent,4),adduct,'Precursor Ion'),
                       FragmentRecord(round(parent-WATER/charge,4),'Precursor-H2O','Common'),
                       FragmentRecord(round(self.ceramide_mass+PROTON-2*WATER,4),'Z0-H2O','Common')))
        result.extend(FragmentRecord(round(parent-n*NH3/charge,4),f'Precursor-{n}NH3','Common') for n in range(1,ammonium+1))
        return sorted(result,key=lambda f:f.mz)

def from_name(name,compound_class=None):
    matched=CHAIN.fullmatch(name)
    if matched is None or compound_class is not None and matched[1]!=compound_class:return None
    cls,series,lc,ld,hydroxy,fc,fd=matched.groups()
    return GangliosideComposition(cls,series,int(lc),int(ld),bool(hydroxy),int(fc),int(fd))

def fragment_type(fragment):
    if fragment.name in REFERENCE_B or fragment.name.startswith('B-') or fragment.name.startswith('B['):return 'Diagnostic_B'
    if fragment.name.startswith(('Y','Z')) and fragment.name!='Z0-H2O':return 'Diagnostic_HG'
    return 'Precursor Ion' if fragment.fragment_type=='Precursor Ion' else 'Common'

def normalize_fragments(compound_class,name,adduct,fragments):
    model=from_name(name,compound_class)
    if model is None:return list(fragments)
    result=[replace(f,fragment_type=fragment_type(f),required_group=None) for f in fragments]
    names={f.name for f in result}
    aliases={REFERENCE_B[f.name] for f in result if f.name in REFERENCE_B}
    for fragment in model.fragments(adduct):
        if fragment.name not in names and fragment.name not in aliases:
            result.append(fragment);names.add(fragment.name)
    return result
