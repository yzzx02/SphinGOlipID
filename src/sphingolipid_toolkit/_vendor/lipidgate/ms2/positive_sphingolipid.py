"""Scope and explicit acyl identities for positive sphingolipid curation."""
from __future__ import annotations

from functools import lru_cache
import re

from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK


CLASS_KEYS = frozenset(rule.lipid_class.upper().replace("-", "")
                       for rule in SPHINGOLIPID_RULEBOOK.values()) | {
    "SPH", "DHSPH", "PHYTOSPH", "CERP",
}
OTHER_MIN_RELATIVE_INTENSITY = 0.005
ACYL_MAX_DOUBLE_BONDS = 3
LCB_MAX_DOUBLE_BONDS = 2
LCB_THREE_DOUBLE_BOND_CARBONS = frozenset({18})
_CHAIN = re.compile(r"(?<![A-Za-z0-9])(?P<prefix>[mdth]|O-|P-)?(?P<carbon>\d+):(?P<db>\d+)", re.IGNORECASE)
_LCB = re.compile(r"(?P<prefix>[mdt])(?P<carbon>\d+):(?P<db>\d+)", re.IGNORECASE)
_SINGLE_LCB_CLASSES = frozenset({"SPB", "SPH", "DHSPH", "PHYTOSPH", "LSM", "S1P", "HEXSPB", "HEX3SPB", "SHEXSPB", "KETOSPB"})
_TWO_CHAIN_SUM = re.compile(r"^(?:Cer|Cer1P|SM)\([mdt]?(?P<carbon>\d+):(?P<db>\d+)\)(?:\(OH\))?$", re.IGNORECASE)


@lru_cache(maxsize=512)
def is_positive_sphingolipid(compound_class: str, adduct: str) -> bool:
    return (str(adduct).strip().endswith("+")
            and str(compound_class).strip().upper().replace("-", "") in CLASS_KEYS)


def explicit_acyl_carbons(chain_name: str) -> tuple[int, ...]:
    """Exclude marked LCBs; include N-acyl, hydroxy-acyl and attached O-acyl chains.

    A sum composition has no explicit acyl chain and is never decomposed here.
    """
    return tuple(int(match["carbon"]) for match in _CHAIN.finditer(chain_name)
                 if (match["prefix"] or "").lower() not in {"m", "d", "t"}
                 and int(match["carbon"]) > 0)


def explicit_chain_compositions(chain_name: str, compound_class: str) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Return actual LCBs and acyls; do not interpret molecular/core sums as LCBs.

    Cer's esterified names may start ``Cerd18:1`` without a delimiter. ASM's
    marked token describes the sphingomyelin core total, followed by an outer
    acyl; it is not an explicitly annotated LCB.
    """
    key = str(compound_class).strip().upper().replace("-", "")
    lcbs = tuple((int(m['carbon']), int(m['db'])) for m in _LCB.finditer(chain_name))
    acyls = tuple((int(m['carbon']), int(m['db'])) for m in _CHAIN.finditer(chain_name)
                  if (m['prefix'] or '').lower() not in {'m', 'd', 't'} and int(m['carbon']) > 0)
    if key == 'ASM':
        return (), acyls
    if len(lcbs) + len(acyls) == 1 and key not in _SINGLE_LCB_CLASSES:
        return (), ()
    return lcbs, acyls


def chain_double_bond_violations(chain_name: str, compound_class: str) -> tuple[str, ...]:
    lcbs, acyls = explicit_chain_compositions(chain_name, compound_class)
    reasons = []
    if any(db > ACYL_MAX_DOUBLE_BONDS for _, db in acyls):
        reasons.append('acyl_double_bonds_above_3')
    if any(db > LCB_MAX_DOUBLE_BONDS + (carbon in LCB_THREE_DOUBLE_BOND_CARBONS) for carbon, db in lcbs):
        reasons.append('lcb_double_bonds_outside_c18_exception')
    return tuple(reasons)


def impossible_two_chain_sum_double_bonds(chain_name: str) -> bool:
    """A necessary total bound only; never assign an acyl or LCB from a sum."""
    match = _TWO_CHAIN_SUM.fullmatch(chain_name)
    return match is not None and int(match['db']) > ACYL_MAX_DOUBLE_BONDS + LCB_MAX_DOUBLE_BONDS + 1
