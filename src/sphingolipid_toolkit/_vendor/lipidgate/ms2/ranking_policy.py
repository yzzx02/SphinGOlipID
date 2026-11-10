from __future__ import annotations

from collections.abc import Callable, Sequence
import re

from .chain_utils import chain_token_multiplicity, extract_fragment_chain_token
from .gate_policy import (
    glyceride_substituent_evidence_group,
    glyceride_substituent_group_multiplicity,
)
from .models import CandidateScore
from .positive_hexosylceramide import CHAIN_CONFIRMATION_CLASSES
from .positive_ganglioside import is_positive_ganglioside as is_positive_extended_ganglioside


def confirmed_sphingolipid_isomer_pair(left: CandidateScore, right: CandidateScore) -> bool:
    """Identify confirmed positive d/h and t alternatives independently of score."""
    a, b = left.record, right.record
    if (a.compound_class not in CHAIN_CONFIRMATION_CLASSES or a.compound_class != b.compound_class
            or (a.adduct != "[M+H]+" and not is_positive_extended_ganglioside(a.compound_class,a.adduct))
            or a.adduct != b.adduct or a.lipid_name != b.lipid_name
            or round(a.precursor_mz, 4) != round(b.precursor_mz, 4)
            or not left.passed_required_gates or not right.passed_required_gates
            or left.resolution_level != "chain_level" or right.resolution_level != "chain_level"):
        return False
    for hydroxy, t in ((a, b), (b, a)):
        match = re.fullmatch(r"([^()]+)\(d(\d+):(\d+)/h(\d+):(\d+)\)", hydroxy.lipid_chain_name, re.I)
        if match and int(match[3]) >= 1:
            expected = f"{match[1]}(t{int(match[2])}:{int(match[3])-1}/{int(match[4])}:{int(match[5])+1})"
            if t.lipid_chain_name.lower() == expected.lower():
                return True
    return False


def confirmed_sphingolipid_isomer_score_tie(left: CandidateScore, right: CandidateScore) -> bool:
    """Compatibility predicate for equal-scoring confirmed isomer pairs."""
    return (confirmed_sphingolipid_isomer_pair(left, right)
            and abs(float(left.total_score)-float(right.total_score)) <= 1e-9)


def multiplicity_adjusted_fragment_count(result: CandidateScore) -> int:
    """Count matched fragments while honoring repeated structural chains."""

    multiplicity = chain_token_multiplicity(result.record)
    substituent_multiplicity = glyceride_substituent_group_multiplicity(result.record)
    if len(substituent_multiplicity) == 1 and sum(substituent_multiplicity.values()) == 3:
        # A single physical peak may support three identical TG chains for the
        # gate, but it remains one observed fragment in a Top1 tie-break.
        return len(result.matched_fragments)
    effective_count = 0
    for match in result.matched_fragments:
        substituent_group = glyceride_substituent_evidence_group(
            result.record,
            match.fragment,
        )
        if substituent_group is not None:
            effective_count += substituent_multiplicity.get(substituent_group, 1)
            continue
        token = extract_fragment_chain_token(match.fragment)
        effective_count += multiplicity.get(token, 1) if token is not None else 1
    return effective_count


def build_original_rank_tiers(
    ordered_results: Sequence[CandidateScore],
    *,
    scores_tied: Callable[[float, float], bool],
    class_key: Callable[[object], str],
) -> list[list[CandidateScore]]:
    """Build immutable original-score tiers before shared-peak rescoring.

    Only the original Top1 score tier is split by fragment count, and that
    split is performed independently inside each lipid subclass.  Lower score
    tiers are never subdivided by fragment count.
    """

    ordered = list(ordered_results)
    if not ordered:
        return []

    highest_score = ordered[0].total_score
    raw_top_tier = [
        result
        for result in ordered
        if scores_tied(result.total_score, highest_score)
    ]
    maximum_fragments_by_class: dict[str, int] = {}
    fragment_counts = {}
    for result in raw_top_tier:
        result_class = class_key(result.record.compound_class)
        count = multiplicity_adjusted_fragment_count(result)
        fragment_counts[id(result)] = count
        maximum_fragments_by_class[result_class] = max(
            maximum_fragments_by_class.get(result_class, 0),
            count,
        )

    top1_tier = [
        result
        for result in raw_top_tier
        if fragment_counts[id(result)]
        == maximum_fragments_by_class[class_key(result.record.compound_class)]
    ]
    top1_ids = {id(result) for result in top1_tier}
    demoted_top_tier = [result for result in raw_top_tier if id(result) not in top1_ids]
    raw_top_ids = {id(result) for result in raw_top_tier}
    demoted_ids = {id(result) for result in demoted_top_tier}
    for result in raw_top_tier:
        setattr(
            result,
            "_demoted_from_top1_by_fragment_count",
            id(result) in demoted_ids,
        )

    tiers: list[list[CandidateScore]] = [top1_tier]
    if demoted_top_tier:
        tiers.append(demoted_top_tier)
    for result in ordered:
        if id(result) in raw_top_ids:
            continue
        if scores_tied(result.total_score, tiers[-1][-1].total_score):
            tiers[-1].append(result)
        else:
            tiers.append([result])
    return tiers
