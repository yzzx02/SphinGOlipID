"""Central author-specified scoring weights and explicit subclass policies."""
from dataclasses import dataclass
from math import isfinite
import re
from types import MappingProxyType


PRIMARY_POOL_SATURATION_HALF_INTENSITY = 0.10
SECONDARY_POOL_SATURATION_HALF_INTENSITY = 0.05
PRECURSOR_CLUSTER_EXCLUSION_DA = 4.1
SUPPORT_SATURATION_COUNT = 3
POOL_WEIGHTS = MappingProxyType({"primary": 60., "secondary": 20., "support": 20.})


@dataclass(frozen=True)
class EvidenceScoringConfig:
    primary_half_intensity: float = PRIMARY_POOL_SATURATION_HALF_INTENSITY
    secondary_half_intensity: float = SECONDARY_POOL_SATURATION_HALF_INTENSITY
    precursor_cluster_exclusion_da: float = PRECURSOR_CLUSTER_EXCLUSION_DA
    minimum_normalized_intensity: float = 0.
    minimum_group_fraction: float = .50

    def __post_init__(self):
        for key in ("primary_half_intensity", "secondary_half_intensity"):
            value = getattr(self, key)
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be finite and positive")
        if not isfinite(self.precursor_cluster_exclusion_da) or self.precursor_cluster_exclusion_da < 0:
            raise ValueError("Precursor exclusion window must be finite and nonnegative")
        if not isfinite(self.minimum_normalized_intensity) or not 0 <= self.minimum_normalized_intensity <= 1:
            raise ValueError("Minimum normalized intensity must be in [0, 1]")
        if not isfinite(self.minimum_group_fraction) or not 0 < self.minimum_group_fraction <= 1:
            raise ValueError("Minimum group fraction must be in (0, 1]")

@dataclass(frozen=True)
class ClassEvidenceRule:
    lipid_class: str
    generator_functions: tuple[str, ...]
    has_HG: bool
    species_groups: tuple[str, ...]
    molecular_groups: tuple[str, ...]
    primary_evidence: str | None
    secondary_evidence: str | None
    chain_count: int | None = 2
    status: str = "provisional_author_policy"
    notes: str = ""

    @property
    def scoring_policy(self):
        return "structural_60_20_20" if self.primary_evidence and self.secondary_evidence else "UNSPECIFIED"

    def __post_init__(self):
        allowed = {"HG", "LCB", "diagnostic_nl", None}
        if self.primary_evidence not in allowed or self.secondary_evidence not in allowed:
            raise ValueError("Structural pools accept only HG, LCB or diagnostic_nl")
        if self.primary_evidence is not None and self.primary_evidence == self.secondary_evidence:
            raise ValueError("Primary and secondary must be different evidence groups")
        if self.chain_count != 2 and (self.primary_evidence is not None or self.secondary_evidence is not None):
            raise ValueError("Single/multichain scoring policies remain UNSPECIFIED")
        if any(group not in {"HG","LCB","NL"} for group in self.species_groups + self.molecular_groups):
            raise ValueError("Only HG, LCB and diagnostic NL can be structural gates")
        if self.molecular_groups and not set(self.species_groups) <= set(self.molecular_groups):
            raise ValueError("Molecular gate must retain every species prerequisite")


_rules = {}
def _add(name, functions, hg, primary, secondary, *, chains=2, notes=""):
    species = ("HG",) if hg else ("LCB",) if name == "So" else ("NL",)
    molecular = tuple(dict.fromkeys((*species,"LCB")))
    _rules[name] = ClassEvidenceRule(name, tuple(functions), hg, species, molecular, primary, secondary, chains, notes=notes)


for _name, _function in {
    "SM":"SM", "PE-Cer":"PE_cer", "PI-Cer":"PI_cer", "PG-Cer":"PG_cer",
    "CAEP":"CAEP", "N-CAEP":"N_CAEP", "Cer1P":"CerP",
}.items():
    _add(_name, [_function], True, "HG", "LCB")
for _name in ("GM1","GM1a","GM1b","GM2","GM3","GD1","GD1a","GD1b","GD2","GD3"):
    _add(_name, ["GSL_fragments"], True, "HG", "LCB", notes="NeuAc HG only when actually generated; branch isomers require explicit structure")
for _name in ("HexCer","GlcCer","GalCer","LacCer","Hex2Cer","Gb3","Gb4","GA1","GA2","type I B antigen"):
    _add(_name, ["GSL_fragments"], False, "diagnostic_nl", "LCB",
        notes="First/complete unhydrated glycan losses are provisional scoring anchors; topology/observability need validation")
_add("Cer", [], False, "LCB", "diagnostic_nl", notes="Author-designated M+H-H2O and M+H-2H2O diagnostic NL; not universal subclass specificity")
for _name in ("FMC_1","FMC_3","FMC_5"):
    _add(_name, [_name], False, "diagnostic_nl", "LCB")
_add("EO-Cer", ["EO_cer"], False, "LCB", "diagnostic_nl", notes="Existing fixed E 18:2 loss is NL, not direct FA-ion confirmation")
_add("EO-GlcCer", ["EO_Glccer"], False, "diagnostic_nl", "LCB")
# Single-chain policy is deliberately unspecified; gates remain independent.
for _name, _source, _hg in (("So","So",False),("S1P","S1P",True),
    ("Lyso-SM","Lyso_SM",True),("Lyso-sulfo","Lyso_sulfo",True),
    ("Glu-So","Glu_So",False),("Gb3-So","Gb3_So",False)):
    _add(_name, [_source], _hg, None, None, chains=1,
        notes="UNSPECIFIED single-chain scoring policy; current so_fuc does not append LCB dictionary fragments")
_rules["1-O-acetyl-Cer"] = ClassEvidenceRule("1-O-acetyl-Cer",("O_FA_cer","OCer"),False,(),(),None,None,3,
    "UNSPECIFIED","Legacy fixed three-chain enumeration needs author review; no complementary two-chain inference")
STRUCTURAL_EVIDENCE_REGISTRY = MappingProxyType(_rules)
CLASS_RULES = STRUCTURAL_EVIDENCE_REGISTRY  # compatibility name for registry callers


def _key(value):
    return re.sub(r"[-_\s]", "", str(value)).lower()


_aliases = {_key(name):name for name in CLASS_RULES}
_aliases.update({_key(alias):name for alias,name in {
    "PE_cer":"PE-Cer","PECer":"PE-Cer","PI_cer":"PI-Cer","PICer":"PI-Cer",
    "PG_cer":"PG-Cer","PGCer":"PG-Cer","CerP":"Cer1P","Cer-1-P":"Cer1P",
    "Gb3Cer":"Gb3","Gb4Cer":"Gb4","GlcSo":"Glu-So","Glu_So":"Glu-So",
    "Gb3So":"Gb3-So","EO_cer":"EO-Cer","EO_Glccer":"EO-GlcCer","O_FA_cer":"1-O-acetyl-Cer","type-I-B":"type I B antigen",
}.items()})


def get_class_rule(lipid_class, registry=None):
    """Unknown subclasses fail closed; custom registries can explicitly override."""
    if registry is not None and lipid_class in registry:
        return registry[lipid_class]
    canonical = _aliases.get(_key(lipid_class))
    if canonical:
        return CLASS_RULES[canonical]
    return ClassEvidenceRule(str(lipid_class),(),False,(),(),None,None,None,"UNSPECIFIED","No reviewed class rule")


def class_mapping_records():
    from .evidence import HG_LABELS, LOSS_TOKENS
    rows = []
    for rule in CLASS_RULES.values():
        def gate(groups):
            return "Precursor AND " + " AND ".join(groups) if groups else "UNSPECIFIED"
        rows.append({"lipid_class":rule.lipid_class,"has_HG":rule.has_HG,
            "species_gate":gate(rule.species_groups),"molecular_species_gate":gate(rule.molecular_groups),
            "primary_evidence":rule.primary_evidence or "UNSPECIFIED",
            "secondary_evidence":rule.secondary_evidence or "UNSPECIFIED",
            "scoring_policy":rule.scoring_policy,
            "scoring_weights":"60 / 20 / 20" if rule.scoring_policy != "UNSPECIFIED" else "UNSPECIFIED",
            "HG examples":"; ".join(label for fn in rule.generator_functions for label in HG_LABELS.get(fn,())),
            "LCB evidence":"LCB_fragment / LCB_mz" if rule.chain_count == 2 else "UNSPECIFIED: no LCB append in current single/three-chain path",
            "NL examples":"M+H-H2O; M+H-2H2O (author policy)" if rule.lipid_class == "Cer" else
                "first/complete glycan loss" if "GSL_fragments" in rule.generator_functions else
                "; ".join(loss for fn in rule.generator_functions for loss in LOSS_TOKENS.get(fn,())),
            "common examples":"[C2H5NO+H]+; other precursor dehydration",
            "generator_functions":"; ".join(rule.generator_functions),"status":rule.status,"notes":rule.notes})
    return rows
