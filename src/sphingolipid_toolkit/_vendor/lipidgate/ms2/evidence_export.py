"""Read-only spectrum evidence snapshot for the result browser.

This adapter does not match, rank or rescore anything. It serializes the exact
CandidateScore matches and the measured peaks used by the search engine.
"""

import json
import re

from .scoring import (
    _fragment_counts_as_hg,
    _fragment_counts_as_fa_loss_gate,
    _pool_for_fragment,
    _is_positive_glyceride_chain_gate_fragment,
    _is_cl_chain_info_fragment,
    _uses_lcb_other_support,
)
from .resolution_policy import fragment_is_chain_evidence
from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK
from .negative_hexcer import is_negative_hexcer, logical_fragment_types
from .positive_hexosylceramide import LCB_POOL_CLASSES, is_positive_hexosylceramide
from .positive_gm3 import LCB_NAME
from .positive_spb import is_positive_spb
from .esterified_ceramide import is_esterified_ceramide
from .positive_cer_sodium import is_positive_cer_sodium, GATE_NAME as CER_SODIUM_GATE
from .positive_cer1p import is_positive_cer1p, GATE_NAMES as CER1P_GATE_NAMES


def fragment_role(record, fragment):
    """Use the same class predicates as the engine, never label by FA text alone."""
    if is_positive_cer1p(record):
        return "headgroup" if fragment.name == CER1P_GATE_NAMES[record.adduct] else "other"
    if (record.compound_class == "Cer" and record.adduct == "[M+H]+"
            and not is_esterified_ceramide(record.compound_class, record.lipid_chain_name)):
        return "chain" if _pool_for_fragment(record, fragment) == "lcb" else "other"
    if is_positive_spb(record):
        return "chain" if _pool_for_fragment(record, fragment) == "lcb" else "other"
    if fragment.name == LCB_NAME and record.compound_class not in LCB_POOL_CLASSES and is_positive_hexosylceramide(record.compound_class, record.adduct):
        return "chain"
    if _fragment_counts_as_hg(record, fragment):
        return "headgroup"
    if is_positive_hexosylceramide(record.compound_class, record.adduct):
        pool = _pool_for_fragment(record, fragment)
        if pool == "b":
            return "headgroup"
        return "chain" if pool == "lcb" else "other"
    # Auxiliary neutral losses are intentionally not a separate display role.
    if fragment.fragment_type == "Neutral_Loss":
        return "other"
    key = f"{record.compound_class}_{record.adduct}"
    if is_esterified_ceramide(record.compound_class, record.lipid_chain_name):
        key = f"Cer-esterified_{record.adduct}"
    rule = SPHINGOLIPID_RULEBOOK.get(key)
    if rule:
        types = (
            logical_fragment_types(fragment)
            if is_negative_hexcer(record.compound_class, record.adduct)
            else {fragment.fragment_type}
        )
        required_types = (
            set().union(*rule.required_type_any_groups)
            if rule.required_type_any_groups
            else set()
        )
        for group, _ in [
            *rule.required_type_count_groups,
            *rule.required_type_fraction_groups,
        ]:
            required_types.update(group)
        for alternatives in rule.required_type_count_any_groups:
            for group, _ in alternatives:
                required_types.update(group)
        names = set(rule.required_all)
        series_match = re.search(
            r"(?:\(|\s)([mdt])\d", record.lipid_chain_name, flags=re.I
        )
        series = series_match.group(1).lower() if series_match else "other"
        for group in [
            *rule.required_any_groups,
            rule.required_any_by_series.get(series, set()),
        ]:
            names.update(group)
        explicit = (
            bool(types & required_types)
            or fragment.name in names
            or any(
                re.search(pattern, fragment.name, flags=re.I)
                for pattern, _ in rule.required_name_pattern_fraction_groups
            )
        )
        if explicit and _pool_for_fragment(record, fragment) in {"fah", "lcb"}:
            return "chain"
        return "other"
    if record.adduct.endswith("-"):
        if fragment.fragment_type == "Diagnostic_FA" or _is_cl_chain_info_fragment(
            record, fragment
        ):
            return "chain"
        return "other"
    if _is_positive_glyceride_chain_gate_fragment(record, fragment):
        return "chain"
    if (
        _fragment_counts_as_fa_loss_gate(record, fragment)
        or fragment_is_chain_evidence(record, fragment)
    ) and _pool_for_fragment(record, fragment) == "fah":
        return "chain"
    return "other"


def evidence_json(spectrum, candidate):
    matches = {match.fragment: match for match in candidate.matched_fragments}
    record = candidate.record
    lcb_support = (_uses_lcb_other_support(record)
                   and "other" in candidate.pool_scores and "lcb" in candidate.pool_scores
                   and candidate.pool_scores["lcb"].total_count > 0)
    cer_sodium = is_positive_cer_sodium(record)
    cer1p = is_positive_cer1p(record)
    include_pool = cer1p or cer_sodium or lcb_support or "b" in candidate.pool_scores or is_positive_hexosylceramide(record.compound_class, record.adduct)
    auxiliary = set(candidate.pool_scores["other"].auxiliary_fragment_indexes) if lcb_support else set()
    credited = (set(candidate.pool_scores["other"].scoring_fragment_indexes)
                | set(candidate.pool_scores["lcb"].scoring_fragment_indexes)) if lcb_support and "lcb" in candidate.pool_scores else set()
    fragments = []
    for index, fragment in enumerate(record.fragments):
        match = matches.get(fragment)
        observed = match.experimental_peak if match else None
        fragments.append(
            dict(
                theoretical_mz=fragment.mz,
                name=fragment.name,
                role=fragment_role(record, fragment),
                fragment_type=fragment.fragment_type,
                observed_mz=observed.mz if observed else None,
                intensity=observed.intensity if observed else None,
                error_ppm=(observed.mz - fragment.mz) / fragment.mz * 1e6
                if observed
                else None,
            )
        )
        if include_pool:
            fragments[-1]["scoring_pool"] = _pool_for_fragment(record, fragment)
        if cer_sodium:
            fragments[-1]["admission_gate"] = fragment.name == CER_SODIUM_GATE
        if cer1p:
            fragments[-1]["admission_gate"] = fragment.name == CER1P_GATE_NAMES[record.adduct]
        if lcb_support:
            pool = _pool_for_fragment(record, fragment)
            fragments[-1]["library_pool"] = pool
            if pool in {"lcb", "other"}:
                fragments[-1]["contributes_to_score"] = index in credited
            if index in auxiliary:
                fragments[-1].update(scoring_pool="other", role="other", auxiliary_lcb_support=True)
    return json.dumps(
        dict(
            schema=1,
            experimental_scope="matching_input",
            spectrum=[[p.mz, p.intensity] for p in spectrum.peaks],
            fragments=fragments,
        ),
        ensure_ascii=False,
        separators=(",", ":"),
    )
