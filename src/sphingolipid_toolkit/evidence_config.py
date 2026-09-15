"""Central author-specified scoring weights and explicit subclass policies."""
from dataclasses import dataclass
from math import isfinite
import re
from types import MappingProxyType


@dataclass(frozen=True)
class ScoringTemplate:
    name: str
    weights: tuple[tuple[str, float], ...]

    def __post_init__(self):
        if len(dict(self.weights)) != len(self.weights) or any(
            category not in {"HG", "LCB", "NL", "common"} or not isfinite(weight) or weight <= 0
            for category, weight in self.weights
        ) or sum(weight for _, weight in self.weights) != 100:
            raise ValueError("Unique MS/MS category weights must be positive and sum to 100")


SCORING_TEMPLATES = MappingProxyType({
    "HG-dominant": ScoringTemplate("HG-dominant", (("HG",60.),("LCB",20.),("NL",15.),("common",5.))),
    "NL-dominant glycosphingolipid": ScoringTemplate("NL-dominant glycosphingolipid", (("LCB",35.),("NL",55.),("common",10.))),
    "LCB-dominant": ScoringTemplate("LCB-dominant", (("LCB",60.),("NL",30.),("common",10.))),
})


@dataclass(frozen=True)
class EvidenceScoringConfig:
    k_hg: float = .10
    k_lcb: float = .10
    k_nl_diagnostic: float = .10
    k_nl: float = .05
    k_common: float = .05
    minimum_normalized_intensity: float = 0.
    minimum_group_fraction: float = .50

    def __post_init__(self):
        for key in ("k_hg","k_lcb","k_nl_diagnostic","k_nl","k_common"):
            value = getattr(self, key)
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be finite and positive")
        if not isfinite(self.minimum_normalized_intensity) or not 0 <= self.minimum_normalized_intensity <= 1:
            raise ValueError("Minimum normalized intensity must be in [0, 1]")
        if not isfinite(self.minimum_group_fraction) or not 0 < self.minimum_group_fraction <= 1:
            raise ValueError("Minimum group fraction must be in (0, 1]")

    def k_for(self, category, strength):
        return {"HG":self.k_hg,"LCB":self.k_lcb,"NL":self.k_nl_diagnostic if strength == "diagnostic" else self.k_nl,
            "common":self.k_common}[category]


@dataclass(frozen=True)
class ClassEvidenceRule:
    lipid_class: str
    generator_functions: tuple[str, ...]
    has_HG: bool
    species_groups: tuple[str, ...]
    molecular_groups: tuple[str, ...]
    scoring_template: str | None
    chain_count: int | None = 2
    status: str = "provisional_author_policy"
    notes: str = ""

    def __post_init__(self):
        if any(group not in {"HG","LCB","NL"} for group in self.species_groups + self.molecular_groups):
            raise ValueError("Only HG, LCB and diagnostic NL can be structural gates")
        if self.molecular_groups and not set(self.species_groups) <= set(self.molecular_groups):
            raise ValueError("Molecular gate must retain every species prerequisite")


_rules = {}
def _add(name, functions, hg, template, *, chains=2, notes=""):
    species = ("HG",) if hg else ("LCB",) if name == "So" else ("NL",)
    molecular = tuple(dict.fromkeys((*species,"LCB")))
    _rules[name] = ClassEvidenceRule(name, tuple(functions), hg, species, molecular, template, chains, notes=notes)


for _name, _function in {
    "SM":"SM", "PE-Cer":"PE_cer", "PI-Cer":"PI_cer", "PG-Cer":"PG_cer",
    "CAEP":"CAEP", "N-CAEP":"N_CAEP", "Cer1P":"CerP",
}.items():
    _add(_name, [_function], True, "HG-dominant")
for _name in ("GM1","GM1a","GM1b","GM2","GM3","GD1","GD1a","GD1b","GD2","GD3"):
    _add(_name, ["GSL_fragments"], True, "HG-dominant", notes="NeuAc HG only when actually generated; branch isomers require explicit structure")
for _name in ("HexCer","GlcCer","GalCer","LacCer","Gb3","Gb4","GA1","GA2","type I B antigen"):
    _add(_name, ["GSL_fragments"], False, "NL-dominant glycosphingolipid",
        notes="First/complete unhydrated glycan losses are provisional scoring anchors; topology/observability need validation")
_add("Cer", [], False, "LCB-dominant", notes="Author-designated M+H-H2O and M+H-2H2O diagnostic NL; not universal subclass specificity")
for _name in ("FMC_1","FMC_3","FMC_5"):
    _add(_name, [_name], False, "NL-dominant glycosphingolipid")
_add("EO-Cer", ["EO_cer"], False, "LCB-dominant", notes="Existing fixed E 18:2 loss is NL, not direct FA-ion confirmation")
_add("EO-GlcCer", ["EO_Glccer"], False, "NL-dominant glycosphingolipid")
# No invented single-chain weights: all three author templates contain category
# assumptions not met by these production functions. Gates are separate and
# still report observed support; confidence remains unavailable until specified.
for _name, _source, _hg in (("So","So",False),("S1P","S1P",True),
    ("Lyso-SM","Lyso_SM",True),("Lyso-sulfo","Lyso_sulfo",True),
    ("Glu-So","Glu_So",False),("Gb3-So","Gb3_So",False)):
    _add(_name, [_source], _hg, None, chains=1,
        notes="UNSPECIFIED single-chain scoring template; current so_fuc does not append LCB dictionary fragments")
_rules["1-O-acetyl-Cer"] = ClassEvidenceRule("1-O-acetyl-Cer",("O_FA_cer","OCer"),False,(),(),None,3,
    "UNSPECIFIED","Legacy fixed three-chain enumeration needs author review; no complementary two-chain inference")
CLASS_RULES = MappingProxyType(_rules)


def _key(value):
    return re.sub(r"[-_\s]", "", str(value)).lower()


_aliases = {_key(name):name for name in CLASS_RULES}
_aliases.update({_key(alias):name for alias,name in {
    "PE_cer":"PE-Cer","PECer":"PE-Cer","PI_cer":"PI-Cer","PICer":"PI-Cer",
    "PG_cer":"PG-Cer","PGCer":"PG-Cer","CerP":"Cer1P","Cer-1-P":"Cer1P",
    "Gb3Cer":"Gb3","Gb4Cer":"Gb4","GlcSo":"Glu-So","Glu_So":"Glu-So",
    "Gb3So":"Gb3-So","EO_cer":"EO-Cer","EO_Glccer":"EO-GlcCer","O_FA_cer":"1-O-acetyl-Cer",
}.items()})


def get_class_rule(lipid_class, registry=None):
    """Unknown subclasses fail closed; custom registries can explicitly override."""
    if registry is not None and lipid_class in registry:
        return registry[lipid_class]
    canonical = _aliases.get(_key(lipid_class))
    if canonical:
        return CLASS_RULES[canonical]
    return ClassEvidenceRule(str(lipid_class),(),False,(),(),None,None,"UNSPECIFIED","No reviewed class rule")


def class_mapping_records():
    from .evidence import HG_LABELS, LOSS_TOKENS
    rows = []
    for rule in CLASS_RULES.values():
        def gate(groups):
            return "Precursor AND " + " AND ".join(groups) if groups else "UNSPECIFIED"
        rows.append({"lipid_class":rule.lipid_class,"has_HG":rule.has_HG,
            "species_gate":gate(rule.species_groups),"molecular_species_gate":gate(rule.molecular_groups),
            "scoring_template":rule.scoring_template or "UNSPECIFIED",
            "HG examples":"; ".join(label for fn in rule.generator_functions for label in HG_LABELS.get(fn,())),
            "LCB evidence":"LCB_fragment / LCB_mz" if rule.chain_count == 2 else "UNSPECIFIED: no LCB append in current single/three-chain path",
            "NL examples":"M+H-H2O; M+H-2H2O (author policy)" if rule.lipid_class == "Cer" else
                "first/complete glycan loss" if "GSL_fragments" in rule.generator_functions else
                "; ".join(loss for fn in rule.generator_functions for loss in LOSS_TOKENS.get(fn,())),
            "common examples":"[C2H5NO+H]+; other precursor dehydration",
            "generator_functions":"; ".join(rule.generator_functions),"status":rule.status,"notes":rule.notes})
    return rows
