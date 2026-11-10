from __future__ import annotations

import re
from dataclasses import replace
from typing import Iterable

from .models import FragmentRecord


CARBON_MONOISOTOPIC_MASS = 12.0
HYDROGEN_MONOISOTOPIC_MASS = 1.00782503223
OXYGEN_MONOISOTOPIC_MASS = 15.99491462
NITROGEN_MONOISOTOPIC_MASS = 14.00307400443
PHOSPHORUS_MONOISOTOPIC_MASS = 30.97376199842
PROTON_MONOISOTOPIC_MASS = 1.007276466621
WATER_MONOISOTOPIC_MASS = 18.01056468
AMMONIA_MONOISOTOPIC_MASS = (
    NITROGEN_MONOISOTOPIC_MASS + 3 * HYDROGEN_MONOISOTOPIC_MASS
)
FORMALDEHYDE_MONOISOTOPIC_MASS = (
    CARBON_MONOISOTOPIC_MASS
    + 2 * HYDROGEN_MONOISOTOPIC_MASS
    + OXYGEN_MONOISOTOPIC_MASS
)
AM_PS_HEADGROUP_NEUTRAL_LOSS_MASS = (
    9 * CARBON_MONOISOTOPIC_MASS
    + 18 * HYDROGEN_MONOISOTOPIC_MASS
    + NITROGEN_MONOISOTOPIC_MASS
    + 11 * OXYGEN_MONOISOTOPIC_MASS
    + PHOSPHORUS_MONOISOTOPIC_MASS
)
PS_HEADGROUP_NEUTRAL_LOSS_MASS = (
    3 * CARBON_MONOISOTOPIC_MASS
    + 8 * HYDROGEN_MONOISOTOPIC_MASS
    + NITROGEN_MONOISOTOPIC_MASS
    + 6 * OXYGEN_MONOISOTOPIC_MASS
    + PHOSPHORUS_MONOISOTOPIC_MASS
)
CE_PE_POSITIVE_HEADGROUP_NEUTRAL_LOSS_MASS = (
    5 * CARBON_MONOISOTOPIC_MASS
    + 12 * HYDROGEN_MONOISOTOPIC_MASS
    + NITROGEN_MONOISOTOPIC_MASS
    + 6 * OXYGEN_MONOISOTOPIC_MASS
    + PHOSPHORUS_MONOISOTOPIC_MASS
)
PE_HEADGROUP_NEUTRAL_LOSS_MASS = (
    2 * CARBON_MONOISOTOPIC_MASS + 8 * HYDROGEN_MONOISOTOPIC_MASS
    + NITROGEN_MONOISOTOPIC_MASS + 4 * OXYGEN_MONOISOTOPIC_MASS
    + PHOSPHORUS_MONOISOTOPIC_MASS
)


def curate_positive_lyso_pe_fragments(
    compound_class: str, precursor_mz: float, adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Keep the distinct LPE-O and LPE-P positive-mode diagnostic patterns."""
    source = list(fragments)
    if adduct != "[M+H]+" or compound_class not in {"LPE-O", "LPE-P"}:
        return source
    precursor = float(precursor_mz)

    def ion(loss: float, name: str, role: str) -> FragmentRecord:
        return synthetic_fragment(round(precursor - loss, 4), name, role)

    if compound_class == "LPE-O":
        return [
            ion(PE_HEADGROUP_NEUTRAL_LOSS_MASS, "[M+H-141]+", "Diagnostic_HG"),
            ion(WATER_MONOISOTOPIC_MASS, "[M+H-H2O]+", "Common"),
            ion(0, "[M+H]+", "Common"),
        ]
    # Morel et al. (2021): loss of the dehydrated glycerol backbone
    # C3H6O2 plus H3PO4 or HPO3. Nitrogen stays in the product ion.
    # "Glycerol backbone" here is 74.0368 Da, not free glycerol (92.0473).
    glycerol_backbone = (
        3 * CARBON_MONOISOTOPIC_MASS + 6 * HYDROGEN_MONOISOTOPIC_MASS
        + 2 * OXYGEN_MONOISOTOPIC_MASS
    )
    phosphoric_acid = (
        3 * HYDROGEN_MONOISOTOPIC_MASS + PHOSPHORUS_MONOISOTOPIC_MASS
        + 4 * OXYGEN_MONOISOTOPIC_MASS
    )
    metaphosphoric_acid = (
        HYDROGEN_MONOISOTOPIC_MASS + PHOSPHORUS_MONOISOTOPIC_MASS
        + 3 * OXYGEN_MONOISOTOPIC_MASS
    )
    return [
        ion(glycerol_backbone + phosphoric_acid, "[M+H-C3H9O6P]+", "Diagnostic_HG"),
        ion(glycerol_backbone + metaphosphoric_acid, "[M+H-C3H7O5P]+", "Diagnostic_HG"),
        ion(glycerol_backbone, "[M+H-C3H6O2]+", "Common"),
        ion(WATER_MONOISOTOPIC_MASS, "[M+H-H2O]+", "Common"),
        ion(0, "[M+H]+", "Common"),
    ]


def curate_positive_pe_p_fragments(
    compound_class: str, lipid_chain_name: str, adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Keep PE-P support as its precursor and P-chain/ethanolamine Common ion."""
    source = list(fragments)
    if compound_class != "PE-P" or adduct != "[M+H]+":
        return source
    source = [fragment for fragment in source if not re.match(
        r"^(?:\(R=O\)\+|RCO\(|\[RCO\]\+)", str(fragment.name or "").strip(),
    )]
    chain = re.search(r"(?:^|[ (])P-(\d+):(\d+)[/_]\d+:\d+(?:\)|$)", lipid_chain_name)
    if chain is None:
        return source
    carbon, double_bonds = map(int, chain.groups())
    if carbon <= 0 or double_bonds < 0 or double_bonds >= carbon:
        return source
    # C(n+2)H(2n+4-2d)N+; use the proton mass with one fewer neutral H.
    mz = round((carbon + 2) * CARBON_MONOISOTOPIC_MASS
               + (2 * carbon + 3 - 2 * double_bonds) * HYDROGEN_MONOISOTOPIC_MASS
               + NITROGEN_MONOISOTOPIC_MASS + PROTON_MONOISOTOPIC_MASS, 4)
    name = f"[P-chain+PEtn-H3PO4]+(P-{carbon}:{double_bonds})"
    retained = [fragment for fragment in source if abs(fragment.mz - mz) > 0.0001]
    existing = next((fragment for fragment in source if abs(fragment.mz - mz) <= 0.0001), None)
    retained.append(replace(existing, name=name, fragment_type="Common", required_group=None)
                    if existing is not None else synthetic_fragment(mz, name, "Common"))
    return retained


def without_pi_ammonium_acyl_losses(
    compound_class: str, adduct: str, fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    source = list(fragments)
    if compound_class != "PI" or adduct != "[M+NH4]+":
        return source
    return [fragment for fragment in source if re.match(
        r'^\[M-\((?:ROOH|R=O)\)\+H\]\+', str(fragment.name or "")
    ) is None]


def required_group_for_fragment(fragment_type: str) -> str | None:
    if fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}:
        return "fah"
    if fragment_type == "Diagnostic_HG":
        return "hg"
    return None


def fragment_with_role(fragment: FragmentRecord, fragment_type: str) -> FragmentRecord:
    return replace(
        fragment,
        fragment_type=fragment_type,
        required_group=required_group_for_fragment(fragment_type),
    )


def synthetic_fragment(mz: float, name: str, fragment_type: str) -> FragmentRecord:
    return FragmentRecord(
        mz=float(mz),
        intensity=100.0,
        name=name,
        fragment_type=fragment_type,
        weight=1.0,
        required_group=required_group_for_fragment(fragment_type),
    )


def fatty_acid_neutral_mass(
    carbons: int,
    double_bonds: int,
    *,
    oxygens: int = 2,
) -> float:
    hydrogens = 2 * int(carbons) - 2 * int(double_bonds)
    return (
        int(carbons) * CARBON_MONOISOTOPIC_MASS
        + hydrogens * HYDROGEN_MONOISOTOPIC_MASS
        + int(oxygens) * OXYGEN_MONOISOTOPIC_MASS
    )


def positive_fatty_acyl_ion_mz(
    carbons: int,
    double_bonds: int,
    *,
    oxygens: int = 2,
) -> float:
    return fatty_acid_neutral_mass(carbons, double_bonds, oxygens=oxygens) - (
        HYDROGEN_MONOISOTOPIC_MASS + OXYGEN_MONOISOTOPIC_MASS
    )


def _chain_parts(token: str) -> tuple[int, int]:
    carbons, double_bonds = str(token).split(":", 1)
    return int(carbons), int(double_bonds)


def _parse_tg_est_chains(name: str) -> tuple[str, str, str, str] | None:
    matched = re.fullmatch(
        r"TG-EST(?:\s+|\()"
        r"(?P<sn1>\d+:\d+)_(?P<sn2>\d+:\d+)_(?P<hfa>\d+:\d+);O"
        r"\(FA\s+(?P<attached>\d+:\d+)\)\)?",
        str(name or "").strip(),
        flags=re.IGNORECASE,
    )
    if matched is None:
        return None
    return (
        matched.group("sn1"),
        matched.group("sn2"),
        matched.group("hfa"),
        matched.group("attached"),
    )


def normalize_tg_o_fragments(
    compound_class: str,
    lipid_chain_name: str,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Assign TG-O ether/acyl losses to its three-substituent evidence pool."""

    source = list(fragments)
    if str(compound_class or "").strip().upper().replace("-", "") != "TGO":
        return source
    if str(adduct or "").strip() not in {"[M+NH4]+", "[M+H]+"}:
        return source
    matched = re.fullmatch(
        r"TG-O\((?P<ether>O-\d+:\d+)_(?P<fa1>\d+:\d+)_(?P<fa2>\d+:\d+)\)",
        str(lipid_chain_name or "").strip(),
        flags=re.IGNORECASE,
    )
    if matched is None:
        return source

    ether_token = matched.group("ether")
    curated: list[FragmentRecord] = []
    for fragment in source:
        compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
        if compact_name == "[M-R1-OH+H]+" or compact_name == f"[M-R1-OH+H]+({ether_token})".upper():
            curated.append(
                fragment_with_role(
                    replace(fragment, name=f"[M-R1-OH+H]+({ether_token})"),
                    "Diagnostic_FA_Loss",
                )
            )
        elif (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and "M-NH3-(ROOH)+NH4" in compact_name
        ):
            curated.append(fragment_with_role(fragment, "Diagnostic_FA_Loss"))
        else:
            curated.append(fragment)
    return sorted(curated, key=lambda fragment: fragment.mz)


def normalize_tg_est_fragments(
    compound_class: str,
    lipid_chain_name: str,
    precursor_mz: float,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Build the TG-EST FA1/FA2/FAHFA gate and ordinary support pool."""

    source = list(fragments)
    if str(compound_class or "").strip().upper().replace("-", "") != "TGEST":
        return source
    adduct_text = str(adduct or "").strip()
    if adduct_text not in {"[M+NH4]+", "[M+H]+"}:
        return source
    parsed = _parse_tg_est_chains(lipid_chain_name)
    if parsed is None:
        return source

    sn1, sn2, hydroxy, attached = parsed
    protonated_mz = (
        float(precursor_mz) - AMMONIA_MONOISOTOPIC_MASS
        if adduct_text == "[M+NH4]+"
        else float(precursor_mz)
    )
    curated = [
        fragment_with_role(fragment, "Precursor Ion")
        for fragment in source
        if fragment.fragment_type == "Precursor Ion"
    ]
    for token in dict.fromkeys((sn1, sn2)):
        carbons, double_bonds = _chain_parts(token)
        name = f"[M+H-FA]+({token})"
        fragment = next(
            (fragment for fragment in source if str(fragment.name or "").strip() == name),
            None,
        ) or synthetic_fragment(
            protonated_mz - fatty_acid_neutral_mass(carbons, double_bonds),
            name,
            "Diagnostic_FA_Loss",
        )
        curated.append(fragment_with_role(fragment, "Diagnostic_FA_Loss"))

    hydroxy_c, hydroxy_db = _chain_parts(hydroxy)
    attached_c, attached_db = _chain_parts(attached)
    fahfa_name = f"[M+H-FAHFA]+({hydroxy};O/FA {attached})"
    fahfa = next(
        (
            fragment
            for fragment in source
            if "O-acyl hydroxy FA" in str(fragment.name or "")
            or str(fragment.name or "").strip().startswith("[M+H-FAHFA]+")
        ),
        None,
    )
    if fahfa is None:
        fahfa_mass = (
            fatty_acid_neutral_mass(hydroxy_c, hydroxy_db, oxygens=3)
            + fatty_acid_neutral_mass(attached_c, attached_db)
            - WATER_MONOISOTOPIC_MASS
        )
        fahfa = synthetic_fragment(
            protonated_mz - fahfa_mass,
            fahfa_name,
            "Diagnostic_FA_Loss",
        )
    curated.append(
        fragment_with_role(replace(fahfa, name=fahfa_name), "Diagnostic_FA_Loss")
    )

    curated.append(synthetic_fragment(protonated_mz, "[M+H]+", "Common"))
    acyl_specs = (
        (sn1, 2, sn1),
        (sn2, 2, sn2),
        (hydroxy, 3, f"{hydroxy};O"),
        (attached, 2, attached),
    )
    seen_names: set[str] = set()
    for token, oxygens, label in acyl_specs:
        name = f"(R=O)+({label})"
        if name in seen_names:
            continue
        seen_names.add(name)
        carbons, double_bonds = _chain_parts(token)
        curated.append(
            synthetic_fragment(
                positive_fatty_acyl_ion_mz(carbons, double_bonds, oxygens=oxygens),
                name,
                "Common",
            )
        )
    unique = {
        (str(fragment.name), round(float(fragment.mz), 6)): fragment
        for fragment in curated
    }
    return sorted(unique.values(), key=lambda fragment: fragment.mz)


def _positive_ps(
    lipid_chain_name: str,
    precursor_mz: float,
    fragments: list[FragmentRecord],
) -> list[FragmentRecord]:
    def compact_name(fragment: FragmentRecord) -> str:
        return re.sub(r"\s+", "", str(fragment.name or "")).upper()

    def is_headgroup(fragment: FragmentRecord) -> bool:
        return compact_name(fragment) in {
            "[M-C3H8O6NP+H]+",
            "[M-C3H8NO6P+H]+",
        }

    def is_chain_ketene(fragment: FragmentRecord) -> bool:
        name = compact_name(fragment)
        return name.startswith("[M-R=O-C3H8O6NP+H]+(") or name.startswith(
            "[M-R=O-C3H8NO6P+H]+("
        )

    curated: list[FragmentRecord] = []
    for fragment in fragments:
        if compact_name(fragment) == "[M-C3H5O2N+H]+":
            continue  # Positive PS has no valid 87 Da loss.
        if is_headgroup(fragment):
            role = "Diagnostic_HG"
        elif is_chain_ketene(fragment):
            role = "Diagnostic_FA_Loss"
        elif fragment.fragment_type == "Precursor Ion":
            role = "Precursor Ion"
        else:
            role = "Common"
        curated.append(fragment_with_role(fragment, role))
    if not any(is_headgroup(fragment) for fragment in curated):
        curated.append(
            synthetic_fragment(
                float(precursor_mz) - PS_HEADGROUP_NEUTRAL_LOSS_MASS,
                "[M-C3H8O6NP+H]+",
                "Diagnostic_HG",
            )
        )

    chain_match = re.fullmatch(
        r"PS\((?P<sn1>\d+:\d+)[_/](?P<sn2>\d+:\d+)\)",
        str(lipid_chain_name or "").strip(),
        flags=re.IGNORECASE,
    )
    existing_tokens = {
        matched.group(1)
        for fragment in curated
        if is_chain_ketene(fragment)
        for matched in [re.search(r"\((\d+:\d+)\)\s*$", str(fragment.name or ""))]
        if matched is not None
    }
    if chain_match is not None:
        for token in dict.fromkeys((chain_match.group("sn1"), chain_match.group("sn2"))):
            if token in existing_tokens:
                continue
            carbons, double_bonds = _chain_parts(token)
            ketene_mass = fatty_acid_neutral_mass(carbons, double_bonds) - WATER_MONOISOTOPIC_MASS
            curated.append(
                synthetic_fragment(
                    float(precursor_mz) - PS_HEADGROUP_NEUTRAL_LOSS_MASS - ketene_mass,
                    f"[M-R=O-C3H8O6NP+H]+({token})",
                    "Diagnostic_FA_Loss",
                )
            )
    return sorted(curated, key=lambda fragment: fragment.mz)


def _positive_am_ps(
    precursor_mz: float,
    fragments: list[FragmentRecord],
) -> list[FragmentRecord]:
    rco_fragments: list[FragmentRecord] = []
    seen_names: set[str] = set()
    for fragment in fragments:
        name = str(fragment.name or "").strip()
        upper_name = name.upper()
        if not (
            upper_name.startswith("(R=O)+(")
            or upper_name.startswith("RCO(")
            or upper_name.startswith("[RCO]+")
        ) or name in seen_names:
            continue
        seen_names.add(name)
        rco_fragments.append(fragment_with_role(fragment, "Diagnostic_FA"))
    rco_fragments.extend(
        [
            synthetic_fragment(
                float(precursor_mz) - AM_PS_HEADGROUP_NEUTRAL_LOSS_MASS,
                "[M-C9H18NO11P+H]+",
                "Diagnostic_HG",
            ),
            synthetic_fragment(
                float(precursor_mz) - 3 * WATER_MONOISOTOPIC_MASS - FORMALDEHYDE_MONOISOTOPIC_MASS,
                "[M-3H2O-HCHO+H]+",
                "Common",
            ),
            synthetic_fragment(float(precursor_mz) - 3 * WATER_MONOISOTOPIC_MASS, "[M-3H2O+H]+", "Common"),
            synthetic_fragment(float(precursor_mz) - 2 * WATER_MONOISOTOPIC_MASS, "[M-2H2O+H]+", "Common"),
            synthetic_fragment(float(precursor_mz) - WATER_MONOISOTOPIC_MASS, "[M-H2O+H]+", "Common"),
            synthetic_fragment(float(precursor_mz), "[M+H]+", "Common"),
        ]
    )
    return sorted(rco_fragments, key=lambda fragment: fragment.mz)


def _naps(fragments: list[FragmentRecord], *, positive: bool) -> list[FragmentRecord]:
    positive_diagnostics = {
        "[DAG-H2O+H]+",
        "[N-acylserine-H2O+H]+",
        "[PA+Na]+",
        "[M+Na-DAG]+",
    }
    curated = []
    for fragment in fragments:
        name = str(fragment.name or "").strip()
        if positive:
            role = "Diagnostic_HG" if name in positive_diagnostics else "Common"
        elif "RCOO" in name.upper():
            role = "Diagnostic_FA"
        elif name == "[PA-H]-":
            role = "Diagnostic_HG"
        else:
            role = "Common"
        curated.append(fragment_with_role(fragment, role))
    return sorted(curated, key=lambda fragment: fragment.mz)


def _negative_ps(fragments: list[FragmentRecord]) -> list[FragmentRecord]:
    curated = []
    for fragment in fragments:
        name = str(fragment.name or "").strip()
        if name in {"[C3H6O5P]-", "[M-C3H5O2N-H]-"}:
            role = "Diagnostic_HG"
        elif name in {"[PO3]-", "[H2PO4]-", "PO3-", "H2PO4-"}:
            role = "Common"
        else:
            role = str(fragment.fragment_type or "Common")
        curated.append(fragment_with_role(fragment, role))
    return sorted(curated, key=lambda fragment: fragment.mz)


def _negative_ce_pe(fragments: list[FragmentRecord]) -> list[FragmentRecord]:
    curated: list[FragmentRecord] = []
    seen_names: set[str] = set()
    for fragment in fragments:
        name = str(fragment.name or "").strip()
        if "RCOO" in name.upper():
            role = "Diagnostic_FA"
        elif name.startswith("[M-(ROOH)-H]-") or name.startswith("[M-(R=O)-H]-"):
            role = "Common"
        elif fragment.fragment_type == "Precursor Ion" or name == "[M-H]-":
            role = "Precursor Ion"
        else:
            continue
        if name in seen_names:
            continue
        seen_names.add(name)
        curated.append(fragment_with_role(fragment, role))
    curated.extend(
        [
            synthetic_fragment(78.9591, "[PO3]-", "Common"),
            synthetic_fragment(152.9953, "[C3H6O5P]-", "Diagnostic_HG"),
            synthetic_fragment(212.0334, "[C5H11NO6P]-", "Diagnostic_HG"),
            synthetic_fragment(268.0601, "[C8H15NO7P]-", "Diagnostic_HG"),
        ]
    )
    return sorted(curated, key=lambda fragment: fragment.mz)


def _positive_ce_pe(
    precursor_mz: float,
    fragments: list[FragmentRecord],
) -> list[FragmentRecord]:
    curated: list[FragmentRecord] = []
    seen_names: set[str] = set()
    for fragment in fragments:
        name = str(fragment.name or "").strip()
        upper_name = name.upper()
        if not (
            upper_name.startswith("(R=O)+(")
            or upper_name.startswith("RCO(")
            or upper_name.startswith("[RCO]+")
        ) or name in seen_names:
            continue
        seen_names.add(name)
        curated.append(fragment_with_role(fragment, "Diagnostic_FA"))
    curated.extend(
        [
            synthetic_fragment(
                float(precursor_mz) - CE_PE_POSITIVE_HEADGROUP_NEUTRAL_LOSS_MASS,
                "[M-C5H12NO6P+H]+",
                "Diagnostic_HG",
            ),
            synthetic_fragment(float(precursor_mz), "[M+H]+", "Common"),
        ]
    )
    return sorted(curated, key=lambda fragment: fragment.mz)


def is_lnape_n_acyl_headgroup(name: str) -> bool:
    return re.fullmatch(r"\[C\d+H\d+O5NP\(\d+:\d+\)\]-", str(name).strip()) is not None


def normalize_special_aminophospholipid_fragments(
    compound_class: str,
    lipid_chain_name: str,
    precursor_mz: float,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Apply class-specific library roles without embedding them in the parser."""

    class_name = str(compound_class or "").strip()
    adduct_text = str(adduct or "").strip()
    source = list(fragments)
    if class_name in {"PtdGlc", "Glc-GP", "6-Ac-Glc-GP"} and adduct_text == "[M-H]-":
        from .negative_ptdglc import chain_pair, negative_glc_gp_fragments

        if chain_pair(lipid_chain_name) is not None:
            return negative_glc_gp_fragments(lipid_chain_name, precursor_mz=precursor_mz)[2]
    if class_name == "LNAPE" and adduct_text.endswith("-"):
        return [
            replace(fragment, name=fragment.name.strip(), fragment_type="Diagnostic_HG", required_group="hg")
            if is_lnape_n_acyl_headgroup(fragment.name) else fragment
            for fragment in source
        ]
    if class_name == "PS" and adduct_text == "[M+H]+":
        return _positive_ps(lipid_chain_name, precursor_mz, source)
    if class_name == "Am-PS" and adduct_text == "[M+H]+":
        return _positive_am_ps(precursor_mz, source)
    if class_name == "NAPS" and adduct_text.endswith("-"):
        return _naps(source, positive=False)
    if class_name == "NAPS" and adduct_text.endswith("+"):
        return _naps(source, positive=True)
    if class_name == "PAsp" and adduct_text == "[M-H]-":
        if not any(str(fragment.name or "").strip() == "[M-H2O-H]-" for fragment in source):
            source.append(
                synthetic_fragment(
                    float(precursor_mz) - WATER_MONOISOTOPIC_MASS,
                    "[M-H2O-H]-",
                    "Common",
                )
            )
        return sorted(source, key=lambda fragment: fragment.mz)
    if class_name == "PS" and adduct_text == "[M-H]-":
        return _negative_ps(source)
    if class_name == "CE-PE" and adduct_text == "[M-H]-":
        return _negative_ce_pe(source)
    if class_name == "CE-PE" and adduct_text == "[M+H]+":
        return _positive_ce_pe(precursor_mz, source)
    return source
