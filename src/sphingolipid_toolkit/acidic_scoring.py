"""Configurable positive sphingolipid scores over reference matches."""
from contextvars import ContextVar
from dataclasses import replace
import json
import math

from ._vendor.lipidgate.ms2.models import FragmentRecord, PoolScore
from ._vendor.lipidgate.ms2.esterified_ceramide import is_esterified_ceramide
from ._vendor.lipidgate.ms2.positive_catalog_glycosphingolipid import PROFILES as GLYCAN_COMPOSITIONS
from ._vendor.lipidgate.ms2.positive_catalog_legacy_sphingolipid import CLASSES as LEGACY_CLASSES
from ._vendor.lipidgate.ms2.positive_ganglioside import CLASSES as GANGLIOSIDE_FAMILIES
from ._vendor.lipidgate.ms2.headgroup_core import GM2_CORE, CURATED
from ._vendor.lipidgate.ms2.positive_gm3 import B1_NAME, GM3_HG, WATER
from ._vendor.lipidgate.ms2.positive_pc_sodium import PHOSPHOCHOLINE
from ._vendor.lipidgate.ms2.matching import match_fragments
from ._vendor.lipidgate.ms2.scoring import (
    _library_fragments_by_pool, _scoring_pools_for_fragment, _uses_lcb_other_support,
    _non_precursor_quality_overrides, _saturation_fragment_quality, _total_score_from_pool_scores,
)
from ._vendor.lipidgate.ms2.search import LipidMS2Searcher
from ._vendor.lipidgate.ms2.sphingolipid_rules import SPHINGOLIPID_RULEBOOK


_TYPED_GLYCAN_HEADGROUPS = GANGLIOSIDE_FAMILIES | {'GD2'} | GLYCAN_COMPOSITIONS.keys()
_NAMED_HEADGROUPS = {'GM2': GM2_CORE, 'GM3': GM3_HG | {B1_NAME, '[Neu5Ac+H-2H2O]+'}}
_HG_LCB_GLYCANS = _TYPED_GLYCAN_HEADGROUPS | _NAMED_HEADGROUPS.keys()
_LINKAGE_GLYCANS = {'GD1a', 'GD1b', 'GT1b', 'GQ1b'}
_PHOSPHOCERAMIDES = {'CAEP', 'N-methylCAEP', 'PG-Cer'}
_CHAIN_PHOSPHOCERAMIDES = {'CAEP', 'PG-Cer'}
_GLYCOSYLATED_BASES = {'HexSPB', 'Hex3SPB', 'SHexSPB'}
_REASSIGNED = {'LSM', 'Hex-IPC', 'Gb4', 'SM', 'Cer1P'} | _LINKAGE_GLYCANS | _GLYCOSYLATED_BASES
_STRUCTURAL_CLASSES = _HG_LCB_GLYCANS | _REASSIGNED
_ORDERED_MATCHING = CURATED | _REASSIGNED | _CHAIN_PHOSPHOCERAMIDES
_THREE_POOLS = ('hg', 'lcb', 'other')

_PROFILE_GROUPS = {
    ('hg', 'other'): (LEGACY_CLASSES - _GLYCOSYLATED_BASES) | _PHOSPHOCERAMIDES | {
        'SM [M+Na]+', 'Cer [M+Na]+', 'Cer1P [M+Na]+', 'S1P', 'MIP2C',
    },
    ('lcb', 'other'): {'Cer', 'SPB'},
    _THREE_POOLS: _STRUCTURAL_CLASSES | {f'{cls} (chain)' for cls in _CHAIN_PHOSPHOCERAMIDES} | {
        'HexCer', 'Hex2Cer', 'Hex3Cer', 'SHexCer', 'SHexCer+O', 'AHexCer', 'PI-Cer', 'PE-Cer',
    },
    ('hg', 'fah', 'other'): {'ASM'},
    ('fah', 'hg', 'other'): {'SL'},
    ('fah', 'lcb', 'other'): {'Esterified Cer'},
}
SCORING_PROFILES = {
    name: (pools, (60, 20, 20) if len(pools) == 3 else (60, 40))
    for pools, names in _PROFILE_GROUPS.items() for name in names
}

_LCB_DEHYDRATIONS = {'LCB-H2O', 'LCB-2H2O'}
_SM_PHOSPHOCHOLINE_LOSS = 'M+H-183'
_SCORED = ContextVar('acidic_candidates', default=None)


def _is_sm_phosphocholine_loss(record, fragment):
    return (record.compound_class == 'SM' and record.adduct == '[M+H]+'
            and fragment.name == _SM_PHOSPHOCHOLINE_LOSS and fragment.fragment_type == 'NL')


def normalize_acidic_weights(overrides):
    normalized = {}
    for cls, values in overrides.items():
        if cls in {'Cer1P', 'SM', 'Gb4'} | _GLYCOSYLATED_BASES and len(values) == 2:
            # Older saved parameters contain explicit HG and Other weights.
            hg, other = values
            if cls in {'SM', 'Cer1P'} and f'{cls} [M+Na]+' not in overrides:
                normalized.update(normalize_acidic_weights({f'{cls} [M+Na]+': (hg, other)}))
            values = (hg, 0, other)
        if cls not in SCORING_PROFILES or len(values) != len(SCORING_PROFILES[cls][0]):
            raise ValueError(f'鞘脂权重与类别 {cls} 的证据池不匹配。')
        values = tuple(float(value) for value in values)
        if any(not math.isfinite(value) or value < 0 for value in values) or not any(values):
            raise ValueError(f'{cls} 权重须为非负有限数字，且至少一项大于零。')
        scale = max(values)
        total = sum(value / scale for value in values)
        weights = tuple(100 * (value / scale) / total for value in values)
        if not all(math.isclose(a, b, rel_tol=1e-12) for a, b in zip(weights, SCORING_PROFILES[cls][1])):
            normalized[cls] = weights
    return normalized


def _structural_pool(fragment, compound_class):
    if compound_class in _TYPED_GLYCAN_HEADGROUPS:
        is_headgroup = fragment.fragment_type in {'Diagnostic_HG', 'Diagnostic_B'}
    elif compound_class in _NAMED_HEADGROUPS:
        is_headgroup = fragment.name in _NAMED_HEADGROUPS[compound_class]
    else:
        is_headgroup = (fragment.fragment_type == 'Diagnostic_HG'
                        or (compound_class in _LINKAGE_GLYCANS and fragment.name.startswith('Reference-HG-'))
                        or (compound_class in {'SM', 'LSM'} and fragment.name == '[C2H6O4P]+')
                        or (compound_class == 'Gb4' and fragment.name == '[HexNAc+Hex+H]+'))
    if is_headgroup:
        return 'hg'
    if compound_class == 'Hex-IPC' and fragment.name in {'Ceramide fragment B', 'Ceramide fragment D', 'Ceramide fragment P'}:
        return 'lcb'
    return 'lcb' if fragment.name in _LCB_DEHYDRATIONS else 'other'


def _structural_score(result, overrides, weights):
    compound_class = result.record.compound_class
    fragments = result.record.fragments
    pool_names = tuple(_structural_pool(fragment, compound_class) for fragment in fragments)
    indexes = {fragment: index for index, fragment in enumerate(fragments)}
    unique = {'hg': {}, 'lcb': {}, 'other': {}}
    claimed = set()
    for pool in unique:
        for match in result.matched_fragments:
            if pool_names[indexes[match.fragment]] != pool:
                continue
            peak = match.experimental_peak
            if peak.mz in claimed:
                continue
            raw = peak.relative_intensity
            if (pool == 'lcb' and raw < .01) or (pool == 'other' and raw <= .005):
                continue
            value = raw if pool == 'other' else overrides.get(id(match), raw)
            best = unique[pool].get(peak.mz)
            option = (value, -abs(match.mz_error), indexes[match.fragment])
            if best is None or option[:2] > best[:2]:
                unique[pool][peak.mz] = option
        claimed.update(unique[pool])

    hg_limit = 2 if compound_class in {'GD2', 'GM2', 'GM3'} else 1
    selected = {name: sorted(unique[name].values(), reverse=True)[:limit]
                for name, limit in (('hg', hg_limit), ('lcb', 1), ('other', 3))}
    for name, weight in zip(('hg', 'lcb', 'other'), weights):
        picked = selected[name]
        n = len(unique[name])
        possible = pool_names.count(name)
        if name == 'other':
            quality = min(n / min(3, possible), 1.) if possible else 0.
        else:
            mean = sum(item[0] for item in picked) / len(picked) if picked else 0.
            quality = _saturation_fragment_quality(mean, .10 if name == 'hg' else .05)
        result.pool_scores[name] = PoolScore(name, n, possible, n / possible if possible else 0.,
                                             quality, weight / 100, quality * weight,
                                             tuple(item[2] for item in picked))
    result.pool_scores.pop('b', None)
    # The common row schema always exports a numeric FAH column.
    result.pool_scores['fah'] = PoolScore('fah', 0, 0, 0., 0., 0., 0.)
    result.total_score = round(_total_score_from_pool_scores(result.pool_scores), 4)


def _score_other_support(result):
    """Credit distinct supporting peaks above 0.5% of the raw base peak."""
    record = result.record
    old = result.pool_scores.get('other')
    if old is None or not old.total_count:
        return
    by_pool = _library_fragments_by_pool(record)
    other = {id(f) for f in by_pool.get('other', ())}
    auxiliary = set()
    protected = set()
    if _uses_lcb_other_support(record):
        auxiliary = {id(f) for f in by_pool.get('lcb', ())}
        lcb = result.pool_scores.get('lcb')
        protected = set(lcb.scoring_fragment_indexes) if lcb else set()
    indexes = {id(fragment): index for index, fragment in enumerate(record.fragments)}
    seen = {match.experimental_peak.mz for match in result.matched_fragments
            if indexes.get(id(match.fragment)) in protected}
    credited, secondary = [], []
    for match in result.matched_fragments:
        key = id(match.fragment)
        peak = match.experimental_peak
        routes = _scoring_pools_for_fragment(record, match.fragment, matched=True)
        index = indexes.get(key)
        if (index is None or peak.relative_intensity <= .005 or peak.mz in seen
                or (key not in other and 'other' not in routes and key not in auxiliary)):
            continue
        seen.add(peak.mz)
        credited.append(index)
        if key in auxiliary:
            secondary.append(index)
    count = len(credited)
    result.pool_scores['other'] = replace(
        old, matched_count=count, count_ratio=count / old.total_count,
        intensity_ratio=min(count / min(old.total_count, 3), 1.0),
        scoring_fragment_indexes=tuple(credited), auxiliary_fragment_indexes=tuple(secondary),
    )


class AcidicSearcher(LipidMS2Searcher):
    """Reuse reference admission and ranking, with an explicit scoring profile."""

    def __init__(self, *args, acidic_weights, **kwargs):
        super().__init__(*args, **kwargs)
        self.acidic_weights = acidic_weights

    def _match_fragments_for_record(self, spectrum, record, experimental_mz=None):
        cls = record.compound_class
        if cls not in _ORDERED_MATCHING or not self._uses_structural_score(record):
            return super()._match_fragments_for_record(spectrum, record, experimental_mz)
        groups = [[], [], []]
        for fragment in record.fragments:
            groups[_THREE_POOLS.index(_structural_pool(fragment, cls))].append(fragment)
        remaining = spectrum.peaks
        matches = []
        for i, group in enumerate(groups):
            if i == 2:
                remaining = [p for p in remaining if p.relative_intensity > .005]
            found = match_fragments(remaining, group, self._fragment_window_da)
            used = {m.experimental_peak.mz for m in found}
            matches.extend(found)
            remaining = [p for p in remaining if p.mz not in used]
        return sorted(matches, key=lambda m: m.fragment.mz)

    @staticmethod
    def _uses_structural_score(record):
        cls = record.compound_class
        if cls in {'Cer1P', 'SM'}:
            return record.adduct == '[M+H]+'
        if cls in _CHAIN_PHOSPHOCERAMIDES:
            return record.adduct == '[M+H]+' and '/' in record.lipid_chain_name
        return record.adduct.endswith('+') and cls in _STRUCTURAL_CLASSES

    def _profile(self, result):
        record = result.record
        cls, adduct = record.compound_class, record.adduct
        if not adduct.endswith('+'):
            return None
        if cls == 'Cer':
            if adduct == '[M+Na]+':
                cls = 'Cer [M+Na]+'
            elif is_esterified_ceramide(cls, record.lipid_chain_name):
                cls = 'Esterified Cer'
        if cls == 'LacCer':
            cls = 'Hex2Cer'
        elif cls in {'SM', 'Cer1P'} and adduct == '[M+Na]+':
            cls += ' [M+Na]+'
        elif cls in _CHAIN_PHOSPHOCERAMIDES and '/' in record.lipid_chain_name:
            cls += ' (chain)'
        profile = SCORING_PROFILES.get(cls)
        return (profile[0], self.acidic_weights.get(cls, profile[1])) if profile else None

    def _apply_profile(self, result, profile, overrides=None):
        pools, weights = profile
        if self._uses_structural_score(result.record):
            _structural_score(result, overrides or {}, weights)
        else:
            _score_other_support(result)
            for pool, weight in zip(pools, weights):
                old = result.pool_scores.get(pool)
                if old is not None:
                    result.pool_scores[pool] = replace(old, weight_ratio=weight / 100,
                                                       pool_score=old.intensity_ratio * weight)
            result.total_score = round(_total_score_from_pool_scores(result.pool_scores), 4)

    def _append_optional_fragment(self, spectrum, result, fragment):
        record = result.record
        result.record = replace(record, fragments=[*record.fragments, fragment])
        used = {match.experimental_peak.mz for match in result.matched_fragments}
        available = [peak for peak in spectrum.peaks if peak.mz not in used]
        extra = match_fragments(available, [fragment], self._fragment_window_da)
        if extra and extra[0].experimental_peak.relative_intensity > 0:
            result.matched_fragments.extend(extra)
            peak = extra[0].experimental_peak
            result.matched_intensity_sum += peak.intensity
            result.matched_relative_intensity_sum += peak.relative_intensity

    def _finalize_candidate(self, spectrum, result):
        record = result.record
        profile = self._profile(result)
        if profile is None or not result.passed_required_gates:
            return result
        overrides = None
        if record.compound_class == 'GM3':
            z0 = next((fragment for fragment in record.fragments if fragment.name == 'M+H-C11H17NO8-2Hex-H2O'), None)
            if z0 is not None:
                fragment = FragmentRecord(round(z0.mz - WATER, 4), 'Z0-H2O', 'Common')
                self._append_optional_fragment(spectrum, result, fragment)
        if record.compound_class == 'SM' and record.adduct == '[M+H]+':
            if not any(_is_sm_phosphocholine_loss(record, fragment) for fragment in record.fragments):
                fragment = FragmentRecord(record.precursor_mz - PHOSPHOCHOLINE, _SM_PHOSPHOCHOLINE_LOSS, 'NL')
                self._append_optional_fragment(spectrum, result, fragment)
                other = result.pool_scores['other']
                result.pool_scores['other'] = replace(other, total_count=other.total_count + 1)
        structural = self._uses_structural_score(record)
        if structural:
            overrides = _non_precursor_quality_overrides(spectrum, result.record, result.matched_fragments)
        self._apply_profile(result, profile, overrides)
        if structural and not result.pool_scores['lcb'].matched_count:
            if record.compound_class in _PHOSPHOCERAMIDES:
                return None
            if result.resolution_level == 'chain_level' and record.compound_class in _REASSIGNED:
                result.resolution_level = 'species_level'
                result.downgrade_reason = 'missing_lcb_chain_evidence'
        if (candidates := _SCORED.get()) is not None:
            candidates.append(result)
        return result

    def _rescore_with_shared_chain_peak_penalty(self, spectrum, result, selected_chain_peak_groups):
        profile = self._profile(result)
        if profile is None:
            return super()._rescore_with_shared_chain_peak_penalty(spectrum, result, selected_chain_peak_groups)
        if not self._uses_structural_score(result.record):
            super()._rescore_with_shared_chain_peak_penalty(spectrum, result, selected_chain_peak_groups)
            return self._apply_profile(result, profile)
        overrides = _non_precursor_quality_overrides(spectrum, result.record, result.matched_fragments)
        matches = [match for match in result.matched_fragments if _structural_pool(match.fragment, result.record.compound_class) == 'lcb'
                   and match.experimental_peak.relative_intensity >= .01]
        counts = {}
        for group in selected_chain_peak_groups:
            shared = [match for match in matches if id(match.experimental_peak) in group]
            if shared:
                strongest = max(shared, key=lambda match: overrides.get(id(match), match.experimental_peak.relative_intensity))
                counts[id(strongest)] = counts.get(id(strongest), 0) + 1
        for match in matches:
            key = id(match)
            if key in counts:
                overrides[key] = overrides.get(key, match.experimental_peak.relative_intensity) * self.SHARED_CHAIN_PEAK_INTENSITY_FACTOR ** counts[key]
        if counts:
            self._apply_profile(result, profile, overrides)

    def score_spectrum(self, spectrum, top_n=3):
        candidates = []
        token = _SCORED.set(candidates)
        try:
            rows = super().score_spectrum(spectrum, top_n)
        finally:
            _SCORED.reset(token)
        lookup = {(result.record.lipid_name, result.record.adduct, result.total_score,
                   self._format_matched_fragments(result.matched_fragments)): result for result in candidates}
        for row in rows:
            key = (row['library_species_name'], row['adduct'], row['total_score'], row['matched_fragments'])
            result = lookup.get(key)
            if result is None:
                continue
            pools, weights = self._profile(result)
            evidence = json.loads(row['ms2_evidence_json'])
            evidence['score_profile'] = {'weights_pct': dict(zip(pools, weights)),
                                         'pool_scores': {p: result.pool_scores[p].intensity_ratio
                                                         for p in pools if p in result.pool_scores}}
            compound_class = result.record.compound_class
            if compound_class in _NAMED_HEADGROUPS:
                hg_pool = result.pool_scores['hg']
                gate = SPHINGOLIPID_RULEBOOK[self._sphingo_rule_key(result.record)]
                minimum_fraction = gate.required_type_fraction_groups[0][1]
                evidence['hg_admission'] = {
                    'matched_channels': hg_pool.matched_count,
                    'library_channels': hg_pool.total_count,
                    'required_channels': max(1, math.ceil(hg_pool.total_count * minimum_fraction)),
                    'minimum_fraction': minimum_fraction,
                }
            if self._uses_structural_score(result.record):
                for i, (fragment, exported) in enumerate(zip(result.record.fragments, evidence['fragments'])):
                    pool = _structural_pool(fragment, compound_class)
                    exported.update(scoring_pool=pool, contributes_to_score=i in result.pool_scores[pool].scoring_fragment_indexes)
                    if compound_class in _TYPED_GLYCAN_HEADGROUPS | _REASSIGNED | _PHOSPHOCERAMIDES:
                        exported['role'] = {'hg': 'headgroup', 'lcb': 'chain', 'other': 'other'}[pool]
            elif compound_class == 'SL':
                for i, (fragment, exported) in enumerate(zip(result.record.fragments, evidence['fragments'])):
                    pool, = _scoring_pools_for_fragment(result.record, fragment, matched=True)
                    exported.update(scoring_pool=pool, role={'hg':'headgroup', 'fah':'chain', 'other':'other'}[pool])
                    if pool == 'other' or result.pool_scores[pool].scoring_fragment_indexes:
                        exported['contributes_to_score'] = i in result.pool_scores[pool].scoring_fragment_indexes
            row['ms2_evidence_json'] = json.dumps(evidence, ensure_ascii=False, separators=(',', ':'))
        return rows
