"""Positive hexosylceramide fragment roles and chain confirmation."""
from __future__ import annotations

from dataclasses import replace
import re

from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,
    HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,
    OXYGEN_MONOISOTOPIC_MASS as O,
    PROTON_MONOISOTOPIC_MASS as PROTON,
    WATER_MONOISOTOPIC_MASS as WATER,
    fragment_with_role,
    synthetic_fragment,
)
from .models import FragmentRecord
from .positive_catalog_glycosphingolipid import CLASSES as CATALOG_GLYCOSPHINGOLIPIDS
from .positive_catalog_phosphoceramide import CLASSES as CATALOG_PHOSPHOCERAMIDES
from .positive_catalog_legacy_sphingolipid import CLASSES as CATALOG_LEGACY_SPHINGOLIPIDS
from .positive_complex_ceramide import LCB_BASE as COMPLEX_LCB_BASE
from .positive_ganglioside import (
    CLASSES as EXTENDED_GANGLIOSIDE_CLASSES,
    is_positive_ganglioside as is_positive_extended_ganglioside,
    fragment_type as extended_ganglioside_fragment_type,
)
from .positive_gb4 import B2_NAME as GB4_B2_NAME, GB4_HG
from .positive_gd2 import B_IONS as GD2_B_IONS, YZ_LOSSES as GD2_YZ
from .positive_gm2 import B_IONS as GM2_B_IONS, HG_LOSSES as GM2_HG
from .positive_gm3 import B_FRAGMENT_TYPE, B_IONS as GM3_B_IONS, GM3_HG, LCB_NAME, upgrade_legacy_spectrum


_GANGLIOSIDE_IONS = {"GD2": (GD2_YZ, GD2_B_IONS), "GM2": (GM2_HG, GM2_B_IONS), "GM3": (GM3_HG, GM3_B_IONS)}
GB_CLASSES = frozenset({"Hex3Cer", "Gb4"})
SHEXCER_CLASSES = frozenset({"SHexCer", "SHexCer+O"})
CLASSES = frozenset({"HexCer", "Hex2Cer", "LacCer", *SHEXCER_CLASSES, *GB_CLASSES, *_GANGLIOSIDE_IONS, *EXTENDED_GANGLIOSIDE_CLASSES, *CATALOG_GLYCOSPHINGOLIPIDS, *CATALOG_PHOSPHOCERAMIDES, *CATALOG_LEGACY_SPHINGOLIPIDS})
CHAIN_CONFIRMATION_CLASSES = CLASSES | {"SM", "Cer", "Cer1P", "PI-Cer", "AHexCer"}
DIHEXOSYL_CLASSES = frozenset({"Hex2Cer", "LacCer"})
LCB_POOL_CLASSES = frozenset({"HexCer", *DIHEXOSYL_CLASSES})
HG_TYPES = frozenset({"Diagnostic_HG", "C类碎片"})
LCB_TYPE = "LCB碎片"
LCB_MIN_RELATIVE_INTENSITY = 0.01
T_LCB_WATER_LOSSES = {"LCB": 0, "LCB-H2O": 1, "LCB-2H2O": 2, "LCB-3H2O": 3}
T_OTHER_NAME = "LCB-CH6O3"
HEXCER_DOUBLE_WATER_LOSS = "M+H-C6H10O5-2H2O"
DIHEXOSYL_REMOVED_LOSS = "M+H-C6H12O6-H2O"
HEXCER_HG = frozenset({
    "M+H-C6H10O5-H2O", HEXCER_DOUBLE_WATER_LOSS,
})
HEXCER_LCB_POOLS = {"d": frozenset({"LCB-H2O", "LCB-2H2O"}), "t": frozenset({"LCB-H2O", "LCB-2H2O"})}
DIHEXOSYL_HG = frozenset({
    "M+H-C12H20O10-H2O", "M+H-C12H20O10-2H2O",
})
DIHEXOSYL_CORE_HG = DIHEXOSYL_HG
HEX3CER_HG = frozenset({
    "M+H-1Hex", "M+H-1Hex-H2O", "M+H-2Hex", "M+H-2Hex-H2O",
    "M+H-3Hex", "M+H-3Hex-H2O",
})
HEX3CER_REMOVED_LOSSES = frozenset({"M+H-1Hex-2H2O", "M+H-2Hex-2H2O"})
SHEXCER_HG = frozenset({"M+H-H2SO4", "M+H-H2SO4-C6H10O5"})
SHEXCER_HYDROXY_HG_ALIAS = "Cer+H-H2O"
SHEXCER_OTHER_LOSS = "M+H-H2SO4-C6H10O5-H2O"
HEXOSE_RESIDUE = 6 * C + 10 * H + 5 * O
T_CHAIN = re.compile(r"(?:\(|\s)t(\d+):(\d+)/", re.IGNORECASE)
HEXCER_BASE = re.compile(r"(?:\(|\s)([dt])(\d+):(\d+)/", re.IGNORECASE)
CHAIN_COMPOSITION = re.compile(r"(?:\(|\s)([mdt])(\d+):(\d+)/(h?)(\d+):(\d+)\)?$", re.IGNORECASE)


def is_positive_hexosylceramide(compound_class: str, adduct: str) -> bool:
    return (adduct == "[M+H]+" and compound_class in CLASSES
            or is_positive_extended_ganglioside(compound_class,adduct))


def is_positive_ganglioside(compound_class: str, adduct: str) -> bool:
    return (adduct == "[M+H]+" and compound_class in _GANGLIOSIDE_IONS
            or is_positive_extended_ganglioside(compound_class,adduct))


def is_positive_globoside(compound_class: str, adduct: str) -> bool:
    return adduct == "[M+H]+" and compound_class in GB_CLASSES


def uses_t_lcb_pool(compound_class: str, chain_name: str) -> bool:
    return compound_class in {"HexCer", "Hex2Cer", "LacCer"} and T_CHAIN.search(chain_name) is not None


def fragment_type(compound_class: str, chain_name: str, fragment: FragmentRecord) -> str:
    name = fragment.name
    if compound_class in CATALOG_GLYCOSPHINGOLIPIDS or compound_class in CATALOG_PHOSPHOCERAMIDES or compound_class in CATALOG_LEGACY_SPHINGOLIPIDS:
        return fragment.fragment_type
    if compound_class in EXTENDED_GANGLIOSIDE_CLASSES:
        return extended_ganglioside_fragment_type(fragment)
    if compound_class == "Cer1P":
        if name == "M+H-H3PO4":
            return "Diagnostic_HG"
        return "Precursor Ion" if fragment.fragment_type == "Precursor Ion" else "Common"
    if compound_class == "Cer":
        base = HEXCER_BASE.search(chain_name)
        if base:
            if name in HEXCER_LCB_POOLS[base.group(1).lower()]:
                return LCB_TYPE
            return "Precursor Ion" if fragment.fragment_type == "Precursor Ion" else "Common"
        return fragment.fragment_type if fragment.fragment_type in {LCB_TYPE, "Precursor Ion"} else "Common"
    if compound_class == "SM":
        return "Common" if fragment.fragment_type == LCB_TYPE or name == "LCB" or name.startswith("LCB-") else fragment.fragment_type
    if compound_class in LCB_POOL_CLASSES:
        if name in (HEXCER_HG if compound_class == "HexCer" else DIHEXOSYL_HG):
            return "Diagnostic_HG"
        if name == "M+H-H2O":
            return "Common"
        base = HEXCER_BASE.search(chain_name)
        if base and name in T_LCB_WATER_LOSSES:
            return LCB_TYPE if name in HEXCER_LCB_POOLS[base.group(1).lower()] else "Common"
    if compound_class == "Hex3Cer":
        if name in HEX3CER_HG:
            return "Diagnostic_HG"
        if name == "M+H-3Hex-2H2O":
            return "Common"
    if compound_class == "Gb4":
        if name in GB4_HG:
            return "Diagnostic_HG"
        if name in (GB4_B2_NAME, LCB_NAME):
            return "Common"
    if compound_class in _GANGLIOSIDE_IONS:
        hg_names, b_names = _GANGLIOSIDE_IONS[compound_class]
        if name in hg_names:
            return "Diagnostic_HG"
        if name in b_names:
            return B_FRAGMENT_TYPE
        if name == LCB_NAME:
            return "Common"
    if compound_class in SHEXCER_CLASSES:
        if name in SHEXCER_HG or (compound_class == "SHexCer+O" and name == SHEXCER_HYDROXY_HG_ALIAS):
            return "Diagnostic_HG"
        base = HEXCER_BASE.search(chain_name)
        if base and name in HEXCER_LCB_POOLS[base.group(1).lower()]:
            return LCB_TYPE
        # Preserve existing support labels without making them structural evidence.
        return fragment.fragment_type if fragment.fragment_type in {"Precursor Ion", "NL", "Neutral_Loss"} else "Common"
    if name == "LCB-2H2O" or (name in T_LCB_WATER_LOSSES and uses_t_lcb_pool(compound_class, chain_name)):
        return LCB_TYPE
    if fragment.fragment_type == LCB_TYPE or name == "LCB" or name.startswith("LCB-"):
        return "Common"
    if fragment.fragment_type in HG_TYPES:
        return "Common"
    return fragment.fragment_type


def normalize_fragments(compound_class, chain_name, precursor_mz, fragments):
    if compound_class in CATALOG_LEGACY_SPHINGOLIPIDS:
        # Total-only extensions must not acquire synthetic chain channels.
        # Preserve the branch-wide removal of the old t-specific support.
        return [f for f in fragments if not (T_CHAIN.search(chain_name) and f.name == T_OTHER_NAME)]
    if compound_class == "GM3":
        legacy = upgrade_legacy_spectrum(chain_name, precursor_mz, fragments)
        if legacy is not None:
            fragments = [fragment_with_role(fragment, fragment.fragment_type) for fragment in legacy]
    t_chain = T_CHAIN.search(chain_name) if uses_t_lcb_pool(compound_class, chain_name) else None
    lcb_base = HEXCER_BASE.search(chain_name) if compound_class in CHAIN_CONFIRMATION_CLASSES else None
    lcb_masses = {}
    if lcb_base:
        series, carbon, double_bonds = lcb_base.groups()
        carbon, double_bonds = int(carbon), int(double_bonds)
        oxygens = 3 if series.lower() == "t" else 2
        intact = carbon * C + (2 * carbon + 3 - 2 * double_bonds) * H + N + oxygens * O + PROTON
        lcb_masses = {name: round(intact - losses * WATER, 4) for name, losses in T_LCB_WATER_LOSSES.items()}
    result = []
    for fragment in fragments:
        name = fragment.name
        if name == T_OTHER_NAME and lcb_base and lcb_base.group(1).lower() == "t":
            continue
        if compound_class in DIHEXOSYL_CLASSES and name == DIHEXOSYL_REMOVED_LOSS:
            continue
        if compound_class == "Hex3Cer" and name in HEX3CER_REMOVED_LOSSES:
            continue
        allowed_lcb = HEXCER_LCB_POOLS[lcb_base.group(1).lower()] if lcb_base else frozenset()
        if (compound_class in _GANGLIOSIDE_IONS or compound_class in GB_CLASSES) and name not in allowed_lcb | {LCB_NAME} and (
            fragment.fragment_type == LCB_TYPE or name == "LCB" or name.startswith("LCB-")
        ):
            continue
        if compound_class in LCB_POOL_CLASSES and name == "LCB" and t_chain is None:
            continue
        role = fragment_type(compound_class, chain_name, fragment)
        if role != fragment.fragment_type:
            fragment = fragment_with_role(fragment, role)
        if name in lcb_masses and (compound_class in LCB_POOL_CLASSES or abs(fragment.mz - lcb_masses[name]) > 0.002):
            fragment = replace(fragment, mz=lcb_masses[name])
        if compound_class == "HexCer" and name == HEXCER_DOUBLE_WATER_LOSS:
            expected_mz = precursor_mz - HEXOSE_RESIDUE - 2 * WATER
            # Allow the original four-decimal rounding of precursor and product.
            if abs(fragment.mz - expected_mz) > 0.0001:
                fragment = replace(fragment, mz=round(expected_mz, 4))
        result.append(fragment)

    if lcb_base:
        names = (T_LCB_WATER_LOSSES if compound_class in LCB_POOL_CLASSES and lcb_base.group(1).lower() == "t"
                 else HEXCER_LCB_POOLS[lcb_base.group(1).lower()])
        existing = {fragment.name for fragment in result}
        for name in sorted(names):
            if name not in existing:
                fragment = synthetic_fragment(lcb_masses[name], name, "Common")
                result.append(fragment_with_role(fragment, fragment_type(compound_class, chain_name, fragment)))
    return result


def has_chain_evidence(record, matches) -> bool:
    """Shared positive sphingolipid chain confirmation, independent of score pools."""
    base = HEXCER_BASE.search(record.lipid_chain_name)
    if base is None and record.compound_class=='AHexCer':base=COMPLEX_LCB_BASE.search(record.lipid_chain_name)
    names = HEXCER_LCB_POOLS.get(base.group(1).lower(), frozenset()) if base else frozenset()
    return any(match.fragment.name in names
               and match.experimental_peak.relative_intensity >= LCB_MIN_RELATIVE_INTENSITY
               for match in matches)


def species_name(record) -> str:
    if "/" not in record.lipid_name:
        return record.lipid_name
    matched = CHAIN_COMPOSITION.search(record.lipid_chain_name)
    if matched is None:
        return record.lipid_name
    series, lcb_c, lcb_db, hydroxy, fa_c, fa_db = matched.groups()
    series = series.lower()
    extra_oh = ""
    if hydroxy:
        if series == "t":
            extra_oh = "(OH)"
        else:
            series = {"m": "d", "d": "t"}[series]
    carbon, double_bonds = int(lcb_c) + int(fa_c), int(lcb_db) + int(fa_db)
    return f"{record.compound_class}({series}{carbon}:{double_bonds}){extra_oh}"
