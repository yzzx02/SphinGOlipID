from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Optional


POOL_FRAGMENT_TYPES = {
    "fah": {"Diagnostic_FA"},
    "hg": {"Diagnostic_HG"},
}


@dataclass(frozen=True)
class ScoreProfile:
    pool_weights: Dict[str, float] = field(
        default_factory=lambda: {"fah": 60.0, "hg": 20.0, "other": 20.0}
    )
    metric_weights: Dict[str, float] = field(
        default_factory=lambda: {"coverage": 0.35, "intensity": 0.65}
    )
    single_group_fallback_min_matches: int = 2
    single_group_fallback_min_relative_intensity_sum: float = 0.5


@dataclass(frozen=True)
class ClassRule:
    lipid_class: str
    required_groups: Dict[str, int] = field(default_factory=dict)
    chain_level_min_fah: int = 2
    allow_hg_only_if_no_fah: bool = False
    positive_hg_min_matches_if_no_fah: int = 1
    positive_hg_min_fraction: float = 0.0
    positive_hg_min_matches: int = 0
    positive_loss_min_fraction: float = 0.0
    positive_loss_min_matches: int = 0
    positive_hg_complete_can_resolve_chain: bool = False
    positive_signature_min_matches_if_no_fah: int = 0
    allow_loss_only_if_no_fah: bool = False
    require_loss_with_fah_only: bool = False
    negative_hg_min_matches: int = 0
    negative_require_all_fah: bool = False
    score_profile: ScoreProfile = field(default_factory=ScoreProfile)


@dataclass
class RuleSet:
    by_class: Dict[str, ClassRule]
    default_rule: ClassRule

    def get(self, lipid_class: str) -> ClassRule:
        return self.by_class.get(lipid_class, self.default_rule)


def build_default_rules(classes: Optional[Iterable[str]] = None) -> RuleSet:
    classes = list(classes or [
        "PC", "PC-P", "PE", "PG", "PI", "PS",
        "LPC", "LPE", "LPG", "LPI", "LPS",
        "PC-O", "PE-O", "PI-O", "PS-O",
        "LPC-O", "LPE-O", "LPE-P", "LPG-O", "LPI-O",
        "PA", "PIP", "PIP2", "PIP3", "CL", "DLCL", "MLCL",
    ])
    default_profile = ScoreProfile()
    lps_profile = ScoreProfile(pool_weights={"fah": 20.0, "hg": 60.0, "other": 20.0})
    default_rule = ClassRule(lipid_class="DEFAULT", score_profile=default_profile)
    by_class = {}
    for lipid_class in classes:
        min_fah = 1 if lipid_class.startswith("L") else 2
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=min_fah,
            allow_hg_only_if_no_fah=lipid_class.startswith("L"),
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=False,
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=lps_profile if lipid_class == "LPS" else default_profile,
        )
    by_class["LNAPE"] = ClassRule(
        lipid_class="LNAPE",
        negative_hg_min_matches=1,
        negative_require_all_fah=True,
        score_profile=ScoreProfile(pool_weights={"fah": 20.0, "hg": 60.0, "other": 20.0}),
    )
    by_class["LPE-P"] = ClassRule(
        lipid_class="LPE-P",
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_min_matches=2,
        positive_hg_min_fraction=1.0,
        positive_hg_complete_can_resolve_chain=True,
        score_profile=default_profile,
    )
    # BMP positive-mode MAG fragments encode both chains, so both MAG peaks
    # are required before we promote the match to chain level.
    by_class["BMP"] = ClassRule(
        lipid_class="BMP",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    # Positive-mode classes such as GPNAE/N-acyl amino lipids often only expose
    # one headgroup-diagnostic fragment in the library, so they still use the
    # original HG-only fallback path.
    for lipid_class in {"GPNAE", "NAGly", "NAGlySer", "NAOrn"}:
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=1,
            allow_hg_only_if_no_fah=True,
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=False,
            positive_signature_min_matches_if_no_fah=3,
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=default_profile,
        )
    by_class["NAE"] = ClassRule(
        lipid_class="NAE",
        chain_level_min_fah=1,
        score_profile=ScoreProfile(pool_weights={"fah": 100.0, "hg": 0.0, "other": 0.0}),
    )
    mg_profile = ScoreProfile(pool_weights={"fah": 0.0, "hg": 75.0, "other": 25.0})
    for lipid_class in {"BA", "BASulfate", "SSulfate", "MG", "NAT"}:
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=1,
            allow_hg_only_if_no_fah=True,
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=lipid_class == "MG",
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=mg_profile if lipid_class == "MG" else default_profile,
        )
    by_class["VD"] = ClassRule(
        lipid_class="VD",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["VE"] = ClassRule(
        lipid_class="VE",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["HBMP"] = ClassRule(
        lipid_class="HBMP",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=True,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    # Positive NAPE/NAPS records expose two complementary structural ions:
    # the DAG ion and the N-acyl headgroup ion.  Both must be present.  NAPS
    # negative mode instead requires every glycerol-chain RCOO ion plus the
    # PA-H ion produced by loss of the N-acyl chain.
    by_class["NAPE"] = ClassRule(
        lipid_class="NAPE",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["NAPS"] = ClassRule(
        lipid_class="NAPS",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        negative_hg_min_matches=1,
        negative_require_all_fah=True,
        score_profile=default_profile,
    )
    # Am-PS positive mode is identified by the 347.0617 Da headgroup loss.
    # RCO+ ions are optional for the class gate but are required to promote a
    # passing result from molecular-species to chain level.  The negative-mode
    # rule retains its strict all-RCOO and two-headgroup gate.
    by_class["Am-PS"] = ClassRule(
        lipid_class="Am-PS",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        negative_hg_min_matches=2,
        negative_require_all_fah=True,
        score_profile=default_profile,
    )
    # PS negative mode uses a two-member HG pool: the glycerophosphate ion
    # (m/z 153) and the serine neutral loss (87.0320 Da).  Either HG is
    # sufficient; phosphate m/z 78 remains supporting evidence only.
    # CE-PE negative mode requires every RCOO chain ion and
    # two of its three headgroup ions; positive mode can fall back to species
    # level from the single class-specific headgroup loss when RCO+ is absent.
    by_class["PS"] = ClassRule(
        lipid_class="PS",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["CE-PE"] = ClassRule(
        lipid_class="CE-PE",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        negative_hg_min_matches=2,
        negative_require_all_fah=True,
        score_profile=default_profile,
    )
    by_class["NAGPS"] = ClassRule(
        lipid_class="NAGPS",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["APCS"] = ClassRule(
        lipid_class="APCS",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_min_fraction=0.5,
        positive_hg_min_matches=1,
        positive_loss_min_fraction=0.5,
        positive_loss_min_matches=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["PHEG"] = ClassRule(
        lipid_class="PHEG",
        required_groups={},
        chain_level_min_fah=2,
        negative_hg_min_matches=2,
        negative_require_all_fah=True,
        score_profile=default_profile,
    )
    # At least one HG and every distinct RCOO ion are mandatory. In the
    # acetylated class, measured glycerophosphate can supply HG without the
    # theoretical glycosyl-loss ion. Neutral losses only support Other.
    for lipid_class in ("Glc-GP", "6-Ac-Glc-GP", "PtdGlc"):
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            chain_level_min_fah=2,
            negative_hg_min_matches=1,
            negative_require_all_fah=True,
            score_profile=default_profile,
        )
    by_class["ADGGA"] = ClassRule(
        lipid_class="ADGGA",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_min_fraction=0.5,
        positive_hg_min_matches=1,
        positive_loss_min_fraction=1.0,
        positive_loss_min_matches=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        negative_require_all_fah=True,
        score_profile=default_profile,
    )
    return RuleSet(by_class=by_class, default_rule=default_rule)


DEFAULT_RULES = build_default_rules()

# Backward-compatible aliases for older callers. The default rules are now
# named by role rather than polarity because scoring is adduct-aware.
build_negative_rules = build_default_rules
DEFAULT_NEGATIVE_RULES = DEFAULT_RULES
