from __future__ import annotations

import re
from dataclasses import dataclass
from dataclasses import replace
from typing import Iterable

from .models import FragmentRecord, LibraryRecord
from .sphingolipid_naming import canonicalize_multichain_sphingolipid_name


FREE_SPHINGOID_BASE_CLASSES = {"SPB", "Sph", "DHSph", "PhytoSph"}
SL_SUPPORTED_POSITIVE_ADDUCTS = {"[M+H]+", "[M+NH4]+"}
SL_HEADGROUP_MZ = 124.0063
SL_NH3_MASS = 17.0265
SL_H2O_MASS = 18.0106
SL_FRAGMENT_MZ_TOLERANCE = 0.02


def _sl_chain_fragment_name(lipid_chain_name: str, fragment_mz: float) -> str:
    match = re.search(r"SL\(m(\d+):(\d+)/(\d+):(\d+)\)", str(lipid_chain_name or ""))
    if match is None:
        return f"SL chain fragment (m/z {float(fragment_mz):.4f})"

    base_c, base_db, acyl_c, acyl_db = map(int, match.groups())
    carbon_step = 14.015650064
    double_bond_step = 2.015650064
    targets = (
        (
            carbon_step * base_c - double_bond_step * base_db + 113.98549936,
            f"NL acyl({acyl_c}:{acyl_db})",
        ),
        (
            carbon_step * base_c
            - double_bond_step * base_db
            + 113.98549936
            - SL_H2O_MASS,
            f"NL acyl({acyl_c}:{acyl_db})-H2O",
        ),
        (
            carbon_step * acyl_c - double_bond_step * acyl_db + 142.01679936,
            f"NL (SPB(m{base_c}:{base_db})-C2H8N)",
        ),
        (
            carbon_step * acyl_c
            - double_bond_step * acyl_db
            + 142.01679936
            - SL_H2O_MASS,
            f"NL (SPB(m{base_c}:{base_db})-C2H8N)-H2O",
        ),
    )
    target_mz, name = min(targets, key=lambda item: abs(float(fragment_mz) - item[0]))
    if abs(float(fragment_mz) - target_mz) <= SL_FRAGMENT_MZ_TOLERANCE:
        return name
    return f"SL chain fragment (m/z {float(fragment_mz):.4f})"


def annotate_positive_sl_fragments(
    lipid_chain_name: str,
    precursor_mz: float,
    adduct: str,
    fragments: Iterable[FragmentRecord],
) -> list[FragmentRecord]:
    """Assign the curated SL positive-mode fragments to HG, FAH, and common pools."""

    result = list(fragments)
    adduct = str(adduct or "").strip()
    if adduct not in SL_SUPPORTED_POSITIVE_ADDUCTS:
        return result

    precursor_mz = float(precursor_mz)
    common_targets = [(precursor_mz - SL_H2O_MASS, "[M+H-H2O]+")]
    if adduct == "[M+NH4]+":
        common_targets = [
            (precursor_mz - SL_NH3_MASS, "[M+H]+"),
            (precursor_mz - SL_NH3_MASS - SL_H2O_MASS, "[M+H-H2O]+"),
        ]

    unassigned_indices: list[int] = []
    annotated: list[FragmentRecord] = []
    for index, fragment in enumerate(result):
        mz = float(fragment.mz)
        if abs(mz - precursor_mz) <= SL_FRAGMENT_MZ_TOLERANCE:
            annotated.append(
                replace(
                    fragment,
                    name=adduct,
                    fragment_type="Precursor Ion",
                    required_group=None,
                )
            )
            continue
        if abs(mz - SL_HEADGROUP_MZ) <= SL_FRAGMENT_MZ_TOLERANCE:
            annotated.append(
                replace(
                    fragment,
                    name="[C2H6NO3S]+",
                    fragment_type="Diagnostic_HG",
                    required_group="hg",
                )
            )
            continue
        common_name = next(
            (
                name
                for target_mz, name in common_targets
                if abs(mz - target_mz) <= SL_FRAGMENT_MZ_TOLERANCE
            ),
            None,
        )
        if common_name is not None:
            annotated.append(
                replace(
                    fragment,
                    name=common_name,
                    fragment_type="Common",
                    required_group=None,
                )
            )
            continue
        unassigned_indices.append(index)
        annotated.append(fragment)

    # Curated MS-DIAL SL entries contain at most four remaining chain fragments.
    # If a future source carries extra experimental peaks, only water-loss pairs
    # are promoted automatically instead of treating every unexplained peak as FAH.
    fah_indices = set(unassigned_indices)
    if len(unassigned_indices) > 4:
        fah_indices = {
            index
            for index in unassigned_indices
            if any(
                other_index != index
                and abs(
                    abs(float(result[index].mz) - float(result[other_index].mz))
                    - SL_H2O_MASS
                )
                <= SL_FRAGMENT_MZ_TOLERANCE
                for other_index in unassigned_indices
            )
        }
    for index in fah_indices:
        fragment = annotated[index]
        annotated[index] = replace(
            fragment,
            name=_sl_chain_fragment_name(lipid_chain_name, fragment.mz),
            fragment_type="Diagnostic_FA_Loss",
            required_group="fah",
        )
    return annotated


@dataclass(frozen=True)
class SourceMspRecord:
    name: str
    precursor_mz: float
    adduct: str
    compound_class: str
    formula: str = ""
    peaks: list[tuple[float, float]] | None = None


def _sum_chains(tokens: list[str]) -> tuple[int, int]:
    carbons = 0
    double_bonds = 0
    for token in tokens:
        m = re.search(r"(\d+):(\d+)", token)
        if not m:
            continue
        carbons += int(m.group(1))
        double_bonds += int(m.group(2))
    return carbons, double_bonds


def transform_source_name(compound_class: str, source_name: str) -> tuple[str, str]:
    """Convert selected MS-DIAL sphingolipid names into the local naming style."""

    cls = str(compound_class or "").strip()
    name = str(source_name or "").strip()

    if cls in FREE_SPHINGOID_BASE_CLASSES:
        if "(" in name:
            canonical_name = f"SPB({name.split('(', 1)[1]}"
        else:
            canonical_name = re.sub(
                r"^(?:SPB|Sph|DHSph|PhytoSph)\s*",
                "SPB",
                name,
                count=1,
            )
        return canonical_name, canonical_name

    if cls == "SL":
        chains = re.findall(r"\d+:\d+(?:;O)?", name)
        if len(chains) >= 2:
            local_chains = []
            has_oh = False
            for index, chain in enumerate(chains[:2]):
                if chain.endswith(";O"):
                    has_oh = True
                    chain = chain[:-2]
                local_chains.append(("m" if index == 0 else "") + chain)
            total_c, total_db = _sum_chains(local_chains)
            suffix = "(OH)" if has_oh and chains[1].endswith(";O") else ""
            return (
                f"SL({local_chains[0]}/{local_chains[1]}){suffix}",
                f"SL(m{total_c}:{total_db}){suffix}",
            )

    if cls == "ASM":
        detail = name
        detail = re.sub(r"^SM\s+", "ASM ", detail)
        canonical = canonicalize_multichain_sphingolipid_name(detail, "ASM")
        return canonical, canonical

    return name, name


def build_library_records(records: Iterable[SourceMspRecord]) -> list[LibraryRecord]:
    out: list[LibraryRecord] = []
    for idx, source in enumerate(records):
        chain_name, species_name = transform_source_name(source.compound_class, source.name)
        compound_class = (
            "SPB"
            if str(source.compound_class or "").strip() in FREE_SPHINGOID_BASE_CLASSES
            else str(source.compound_class)
        )
        fragments = [
            FragmentRecord(
                mz=float(mz),
                intensity=float(intensity),
                name=f"m/z {float(mz):.4f}",
                fragment_type="Common",
                required_group=None,
            )
            for mz, intensity in sorted(source.peaks or [], key=lambda item: item[0])
        ]
        if compound_class == "SL":
            fragments = annotate_positive_sl_fragments(
                chain_name,
                float(source.precursor_mz),
                str(source.adduct),
                fragments,
            )
        out.append(
            LibraryRecord(
                record_id=idx,
                compound_class=compound_class,
                lipid_name=species_name,
                lipid_chain_name=chain_name,
                precursor_mz=float(source.precursor_mz),
                adduct=str(source.adduct),
                formula=str(source.formula or ""),
                polarity="+" if str(source.adduct).strip().endswith("+") else "-",
                fragments=fragments,
            )
        )
    return out
