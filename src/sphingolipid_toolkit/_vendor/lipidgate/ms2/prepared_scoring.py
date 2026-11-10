"""Partial evaluation of ordinary phospholipid policy, without changing rules.

Library-only branch choices and chain-count predicates are prepared once in the
bounded fact scope. Each query still computes its own evidence, final gates,
fallbacks, resolution and numerical score. Unusual inputs use the reference.
"""
import math

from . import scoring as reference
from .chain_utils import chain_token_counts, extract_chain_tokens, extract_fragment_chain_token, FRAGMENT_CHAIN_TOKEN_RE, canonical_chain_token
from .models import CandidateScore, FragmentRecord, LibraryRecord
from .negative_lpc import is_negative_lpc
from .positive_lpc_ether import is_positive_ether_lpc
from .native_bridge import load_bridge
from .record_facts import cached_rule_fact, with_record_facts, get_native_query
from . import native_engine
from .scoring_policy import (FAH_ONLY_FALLBACK_MIN_MATCHES, FAH_ONLY_FALLBACK_MIN_POOL_SCORE,
    FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE, HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY)
from .resolution_policy import (
    get_chain_evidence_profile, required_chain_evidence_hits, fragment_is_chain_evidence,
    POSITIVE_DIACYL_CHAIN_CONFIRMATION_CLASSES, SINGLE_CHAIN_COVERAGE_REASON,
    ChainEvidenceProfile,
)


CLASSES = frozenset({'PC','PCO','PCP','PE','PEO','PG','PGO','PI','PIO','PA','PAO','PS','PSO',
                     'LPC','LPCO','LPE','LPEO','LPEP','LPG','LPGO','LPI','LPIO','LPA','LPAO','LPS'})
ADDUCTS = frozenset({'[M+H]+','[M-H]-','[M+CH3COO]-','[M+HCOO]-'})
_NATIVE_RESOLUTIONS = ('class_level', 'species_level', 'chain_level',
                       'tentative_species_level', 'tentative_chain_level')


def _decision_signature(plan):
    profile = plan.get('profile')
    return (FAH_ONLY_FALLBACK_MIN_MATCHES, FAH_ONLY_FALLBACK_MIN_POOL_SCORE,
            FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE,
            reference.POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY,
            HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY,
            None if profile is None else (profile.applies_without_library_evidence,
                                           profile.species_fallback_without_chain))


def _decision_configuration(plan):
    if plan.get('glyceride'):
        return None
    profile = plan['profile']
    if profile is not None and type(profile) is not ChainEvidenceProfile:
        return None
    result = {name: plan[name] for name in (
        'positive', 'ratio_uses_fa', 'loss_resolves', 'confirm', 'allow_hg', 'hg_complete',
        'single_species', 'hg_dominant', 'support_scoring', 'lyso', 'loss_only', 'require_loss',
        'hg_total', 'loss_total', 'hg_required', 'negative_hg_required', 'fah_required',
        'chain_fah_required', 'loss_required', 'ratio_required', 'profile_required',
        'expected_chain_count', 'signature_required')}
    result.update(has_fah=bool(plan['fah']), has_tokens=bool(plan['tokens']),
        has_profile=profile is not None,
        profile_applies=False if profile is None else profile.applies_without_library_evidence,
        profile_species_fallback=False if profile is None else profile.species_fallback_without_chain,
        fah_min_matches=FAH_ONLY_FALLBACK_MIN_MATCHES,
        fah_min_pool=FAH_ONLY_FALLBACK_MIN_POOL_SCORE,
        fah_hg_dominant_min_pool=FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE,
        fah_min_quality=reference.POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY,
        hg_min_relative=HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY)
    return result if all(type(v) in (bool, int, float) for v in result.values()) else None


def _decision_reason(code, missing):
    if code == 9:
        return 'missing_required_' + ','.join(missing)
    return ('', reference.PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON, SINGLE_CHAIN_COVERAGE_REASON,
            reference.CHAIN_LEVEL_INFO_MISSING_REASON, 'missing_chain_annotation',
            'lyso_hg_only_fallback', 'missing_positive_gate', reference.FAH_ONLY_FALLBACK_REASON,
            reference.HG_ONLY_FALLBACK_REASON)[code]


def _mask(record, predicate):
    return sum(1 << i for i, f in enumerate(record.fragments) if predicate(f))


def _count_program(record, predicate):
    multiplicity = chain_token_counts(record)
    groups, unassigned = {}, 0
    for i, fragment in enumerate(record.fragments):
        if not predicate(record, fragment):
            continue
        token = extract_fragment_chain_token(fragment)
        if token is not None and token in multiplicity:
            groups[token] = groups.get(token, 0) | (1 << i)
        else:
            unassigned |= 1 << i
    return (sum(multiplicity.values()) if multiplicity else -1, unassigned,
            tuple((mask, multiplicity[token]) for token, mask in groups.items()))


def _count(mask, program):
    cap, unassigned, groups = program
    value = (mask & unassigned).bit_count() + sum(weight for group, weight in groups if mask & group)
    return min(value, cap) if cap >= 0 else value


def _glyceride_resolution_is_implied(record, derived, groups):
    """Prove the reference chain-loss branch from mandatory library groups.

    A group implies one chain token only if every possible ion in that group
    carries that token. Other records still run the complete resolution code.
    """
    tokens = extract_chain_tokens(record.lipid_chain_name)
    if not tokens or not derived[5] or get_chain_evidence_profile(record) is not None:
        return False
    cap, _, count_groups = _count_program(record, reference._fragment_counts_as_chain_resolving_loss)
    minimum = sum(weight for mask, weight in count_groups
                  if any(group and group & mask == group for group in groups))
    if cap >= 0:
        minimum = min(minimum, cap)
    return minimum >= reference._positive_required_loss_hits(tokens, derived[3])


def _attach_native(record, rule, plan):
    kernel = load_bridge()
    if kernel is not None and hasattr(kernel, 'create_pool_plan'):
        configuration, routes, groups = reference._numeric_pool_plan(record, rule)
        statistics = (plan.get('programs', ()), plan['routes'], plan.get('hg_mask', 0),
                      plan.get('chain_mask', 0), plan.get('expected_chain_masks', ()),
                      plan.get('raw_loss_mask', 0), plan.get('signature_mask', 0),
                      _mask(record, lambda f: f.fragment_type == 'Precursor Ion'))
        decision = _decision_configuration(plan) if hasattr(kernel, 'score_pool_policy') else None
        capsule = (kernel.create_pool_plan(routes, configuration, groups, statistics, decision)
                   if decision is not None else kernel.create_pool_plan(routes, configuration, groups, statistics))
        if capsule is not None:
            plan['native'] = (capsule, configuration, bytes(2048 + len(routes)*64 +
                sum(64+32*len(p[2]) for p in plan.get('programs', ()))))
            if decision is not None:
                plan['native_policy_signature'] = _decision_signature(plan)
            plan['indexed'] = (type(record.precursor_mz) is float and all(
                type(f) is FragmentRecord and type(f.mz) is float and type(f.name) is str
                and type(f.fragment_type) is str and type(f.required_group) in (str,type(None))
                and type(f.intensity) in (float,int) and type(f.weight) in (float,int)
                for f in record.fragments))
    return plan


@cached_rule_fact
def _prepare(record, rule):
    r = reference
    key = r._normalized_class_key(record.compound_class)
    glyceride = r.is_positive_three_substituent_glyceride(record)
    if ((not glyceride and (key not in CLASSES or record.adduct not in ADDUCTS)) or len(record.fragments) > 64
            or len({id(f) for f in record.fragments}) != len(record.fragments)
            or is_negative_lpc(record) or is_positive_ether_lpc(record)
            or r._requires_precursor_fragment(record) or r.oxidized_chain_requirements(record) is not None):
        return None
    derived = r._derive_required_groups(record, rule)
    fah, hg_total, loss_total, fa_loss_total, positive, loss_resolves, lyso, loss_only, require_loss, supports, names = derived
    if not supports and names:
        return None
    if glyceride:
        if hg_total or lyso:
            return None
        groups = r.glyceride_substituent_group_multiplicity(record)
        if not groups:
            return None
        masks = tuple(_mask(record,lambda f:r.glyceride_substituent_evidence_group(record,f)==group)
                      for group in groups)
        return _attach_native(record, rule, dict(glyceride=True, derived=derived,
                    routes={id(f):1 << i for i,f in enumerate(record.fragments)},
                    groups=masks, implied_chain_resolution=_glyceride_resolution_is_implied(record, derived, masks)))
    tokens = extract_chain_tokens(record.lipid_chain_name)
    profile = get_chain_evidence_profile(record)
    expected = set(fah)
    multiplicity = chain_token_counts(record)
    fa_groups = {}
    for i, fragment in enumerate(record.fragments):
        token = extract_fragment_chain_token(fragment)
        if fragment.fragment_type == 'Diagnostic_FA' and token is not None and (not expected or token in expected):
            fa_groups[token] = fa_groups.get(token, 0) | (1 << i)
    fah_program = (-1, 0, tuple((mask, multiplicity.get(token, 1) if positive else 1)
                              for token, mask in fa_groups.items()))
    confirm = (positive and key in POSITIVE_DIACYL_CHAIN_CONFIRMATION_CLASSES and len(tokens) == 2
               and not any(t.upper().startswith(('O-', 'P-')) for t in tokens))
    expected_chains = {canonical_chain_token(m) for m in FRAGMENT_CHAIN_TOKEN_RE.finditer(record.lipid_chain_name)}
    chain_mask = _mask(record, lambda f: fragment_is_chain_evidence(record, f))
    return _attach_native(record, rule, dict(
        positive=positive, fah=fah, hg_total=hg_total, loss_total=loss_total,
        hg_required=r._required_positive_hg_hits(rule, fah, hg_total),
        negative_hg_required=r._required_negative_hg_hits(record, hg_total, rule),
        fah_required=len(set(fah)) if rule.negative_require_all_fah else r._required_fah_gate_hits(record, fah),
        chain_fah_required=r._required_chain_level_fah_hits(record, fah, tokens, rule),
        loss_required=r._positive_required_loss_hits(tokens, fa_loss_total),
        ratio_required=r._required_fractional_hits(fa_loss_total if fa_loss_total>0 else loss_total,
                                                 rule.positive_loss_min_fraction, rule.positive_loss_min_matches),
        ratio_uses_fa=fa_loss_total>0, loss_resolves=loss_resolves,
        tokens=tokens, profile=profile, profile_required=required_chain_evidence_hits(record, profile) if profile else 0,
        confirm=confirm, expected_chain_count=len(expected_chains),
        expected_chain_masks=tuple(_mask(record, lambda f: fragment_is_chain_evidence(record, f)
                                        and extract_fragment_chain_token(f)==token) for token in expected_chains),
        routes={id(f): 1 << i for i, f in enumerate(record.fragments)},
        chain_mask=chain_mask, hg_mask=_mask(record, lambda f: r._fragment_counts_as_hg(record, f)),
        raw_loss_mask=_mask(record, lambda f: f.fragment_type=='Diagnostic_FA_Loss'),
        signature_mask=_mask(record, lambda f: f.fragment_type not in {'Diagnostic_HG',r.CANDIDATE_HG_FRAGMENT_TYPE,'Precursor Ion'}),
        programs=(fah_program, _count_program(record,r._fragment_counts_as_effective_loss),
                  _count_program(record,r._fragment_counts_as_fa_loss_gate),
                  _count_program(record,r._fragment_counts_as_chain_resolving_loss),
                  _count_program(record,fragment_is_chain_evidence)),
        lyso=lyso, loss_only=loss_only, require_loss=require_loss,
        allow_hg=rule.allow_hg_only_if_no_fah, hg_complete=rule.positive_hg_complete_can_resolve_chain,
        signature_required=min(max(int(rule.positive_signature_min_matches_if_no_fah),0),
                               r._record_positive_signature_fragment_count(record)),
        single_species=r._positive_single_chain_species_hg_only(tokens,positive,lyso,fah),
        hg_dominant=r._is_positive_hg_dominant_phospholipid(record),
        support_scoring=r._uses_positive_support_pool_scoring(record)))


def score_prepared(spectrum, record, rule, matches, tolerance_da, tolerance_ppm):
    r = reference
    if r._uses_lcb_other_support(record):
        return NotImplemented
    plan = _prepare(record, rule)
    if plan is None:
        return NotImplemented
    ppm = ((spectrum.precursor_mz-record.precursor_mz)/record.precursor_mz)*1e6
    if (abs(spectrum.precursor_mz-record.precursor_mz)>tolerance_da if tolerance_da is not None else abs(ppm)>tolerance_ppm):
        return NotImplemented
    if not matches:
        return None
    pools, summary = None, None
    if plan.get('native') is not None:
        bridge = load_bridge()
        normalization = (r._native_pool_normalization(spectrum, record)
                         if hasattr(bridge, 'pool_metrics_normalized')
                         and hasattr(bridge, 'score_pool_policy_normalized') else None)
        overrides = None if normalization is not None else r._non_precursor_quality_overrides(spectrum, record, matches)
        if ('native_policy_signature' in plan and
                _decision_signature(plan) == plan['native_policy_signature']):
            native = (bridge.score_pool_policy_normalized(matches, plan['native'][0], normalization)
                      if normalization is not None else bridge.score_pool_policy(matches, plan['native'][0], overrides))
            if native is False:
                return None
            if native is not None:
                metrics, missing_mask, resolution_code, reason_code = native
                order = ((4, 'fah'), (1, 'hg'), (2, 'loss')) if not plan['positive'] else ((1, 'hg'), (2, 'loss'), (4, 'fah'))
                missing = (['chain_confirmation'] if missing_mask & 8 else
                           [name for bit, name in order if missing_mask & bit])
                pools = {item[0]: r.PoolScore(item[0],values[0],item[1],values[1],values[2],0.,values[3])
                         for item,values in zip(plan['native'][1],metrics)}
                return CandidateScore(record,round(r._total_score_from_pool_scores(pools),4),not missing,missing,ppm,
                    _NATIVE_RESOLUTIONS[resolution_code],list(matches),pools,
                    *r._matched_intensity_sums(spectrum,matches),_decision_reason(reason_code,missing))
        native = (bridge.pool_metrics_normalized(matches, plan['native'][0], normalization)
                  if normalization is not None else bridge.pool_metrics(matches, plan['native'][0], overrides))
        if native is None and normalization is not None:
            overrides = r._non_precursor_quality_overrides(spectrum, record, matches)
            native = bridge.pool_metrics(matches, plan['native'][0], overrides)
        if native is not None:
            metrics, summary = native
            pools = {item[0]: r.PoolScore(item[0],values[0],item[1],values[1],values[2],0.,values[3])
                     for item,values in zip(plan['native'][1],metrics)}
    if plan.get('glyceride'):
        mask = summary[0] if summary is not None else 0
        if summary is None:
            for match in matches:
                bit = plan['routes'].get(id(match.fragment))
                if bit is None:
                    return NotImplemented
                mask |= bit
        if not all(mask & group for group in plan['groups']):
            return None
        if pools is None:
            pools = r._calculate_pool_scores(matches,record,rule,spectrum=spectrum)
        resolution, reason = (('chain_level','') if plan['implied_chain_resolution'] else
                              r._determine_resolution(record,matches,rule,plan['derived'][0],True,plan['derived'][5]))
        return CandidateScore(record,round(r._total_score_from_pool_scores(pools),4),True,[],ppm,
            resolution,list(matches),pools,*r._matched_intensity_sums(spectrum,matches),reason)
    if summary is not None:
        mask, counts, hg_count, raw_loss_count, signature_count, hg_max, chain_peak_count, observed = summary
        fah_count, loss_count, fa_loss_count, chain_loss_count, profile_count = counts
    else:
        mask, hg_max, has_hg, chain_peaks = 0, 0., False, set()
        for match in matches:
            bit = plan['routes'].get(id(match.fragment))
            if bit is None:
                return NotImplemented
            mask |= bit
            if bit & plan['hg_mask']:
                relative = float(match.experimental_peak.relative_intensity)
                if not math.isfinite(relative):
                    return NotImplemented
                if not has_hg or relative > hg_max:
                    hg_max = relative
                has_hg = True
            if bit & plan['chain_mask']:
                chain_peaks.add(float(match.experimental_peak.mz))
        fah_count, loss_count, fa_loss_count, chain_loss_count, profile_count = (
            _count(mask, p) for p in plan['programs'])
        hg_count = (mask & plan['hg_mask']).bit_count()
        signature_count = (mask & plan['signature_mask']).bit_count()
        chain_peak_count = len(chain_peaks)
        observed = sum(bool(mask & m) for m in plan['expected_chain_masks'])
    missing = []
    if plan['positive']:
        if plan['hg_total']>0 and hg_count<plan['hg_required']:
            missing.append('hg')
        if plan['hg_total']==0 and plan['loss_total']>0 and loss_count<plan['loss_required']:
            missing.append('loss')
        if (plan['hg_total']>0 and plan['loss_total']>0 and plan['ratio_required']>0
                and (fa_loss_count if plan['ratio_uses_fa'] else loss_count)<plan['ratio_required']):
            missing.append('loss')
        if plan['hg_total']==0 and plan['loss_total']==0 and plan['fah'] and fah_count<plan['fah_required']:
            missing.append('fah')
    else:
        if plan['fah'] and fah_count<plan['fah_required']:
            missing.append('fah')
        if plan['hg_total']>0 and hg_count<plan['negative_hg_required']:
            missing.append('hg')
        if (plan['loss_only'] or plan['require_loss']) and loss_count<1:
            missing.append('loss')
    if missing:
        recoverable = (set(missing).issubset({'hg','sphingo_rule'}) or set(missing).issubset({'fah','loss'})) if plan['positive'] else missing==['hg']
        if not recoverable:
            return None
    chain_decision = None
    if plan['confirm'] and mask & plan['chain_mask']:
        if chain_peak_count<2:
            chain_decision = ('tentative_species_level',r.PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON)
        elif observed<2 or observed<plan['expected_chain_count']:
            chain_decision = ('chain_level',SINGLE_CHAIN_COVERAGE_REASON)
        else:
            chain_decision = ('chain_level','')
    configured = chain_decision
    profile = plan['profile']
    if configured is None and profile is not None:
        if plan['profile_required']>0 and profile_count>=plan['profile_required']:
            configured = ('chain_level','')
        elif ((plan['profile_required']>0 or profile.applies_without_library_evidence)
              and profile.species_fallback_without_chain):
            configured = ('species_level',r.CHAIN_LEVEL_INFO_MISSING_REASON)
    if not plan['tokens']:
        resolution, reason = 'species_level','missing_chain_annotation'
    elif configured is not None and configured[0]=='chain_level':
        resolution, reason = configured
    elif plan['positive']:
        if plan['fah'] and fah_count>=plan['chain_fah_required']:
            resolution, reason = 'chain_level',''
        elif plan['loss_resolves'] and chain_loss_count>=plan['loss_required']:
            resolution, reason = 'chain_level',''
        elif hg_count>=1:
            if configured is not None:
                resolution, reason = configured
            elif not plan['fah'] and not plan['loss_total']:
                if (not plan['allow_hg'] or plan['hg_complete'] or
                        (plan['signature_required']>0 and signature_count>=plan['signature_required'])):
                    resolution, reason = 'chain_level',''
                else:
                    resolution, reason = 'class_level','lyso_hg_only_fallback'
            else:
                resolution, reason = 'species_level',r.CHAIN_LEVEL_INFO_MISSING_REASON
        elif loss_count>=1:
            resolution, reason = 'species_level',r.CHAIN_LEVEL_INFO_MISSING_REASON
        else:
            resolution, reason = 'class_level','missing_positive_gate'
    elif not plan['fah'] and plan['loss_only']:
        resolution, reason = ('chain_level','') if mask & plan['raw_loss_mask'] else ('species_level',r.CHAIN_LEVEL_INFO_MISSING_REASON)
    else:
        resolution, reason = ('chain_level','') if fah_count>=plan['chain_fah_required'] else ('species_level',r.CHAIN_LEVEL_INFO_MISSING_REASON)
    if pools is None:
        pools = r._calculate_pool_scores(matches,record,rule,spectrum=spectrum)
    fah, hg = pools['fah'], pools['hg']
    fah_low = (plan['positive'] and bool(missing) and set(missing).issubset({'hg','sphingo_rule'})
               and fah.total_count>0 and hg.total_count>0 and hg.matched_count==0
               and fah.matched_count>=FAH_ONLY_FALLBACK_MIN_MATCHES and (not plan['support_scoring'] or fah.intensity_ratio>=r.POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY)
               and fah.pool_score >= (FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE if plan['hg_dominant'] else FAH_ONLY_FALLBACK_MIN_POOL_SCORE))
    hg_low = (resolution!='chain_level' and plan['hg_dominant'] and set(missing).issubset({'fah','loss'})
              and hg.total_count>0 and hg.matched_count>0 and hg_max>=HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY)
    weak_hg = resolution!='chain_level' and not missing and plan['hg_dominant'] and hg_count>0 and not hg_low
    if fah_low:
        missing, resolution, reason = [], 'tentative_chain_level',r.FAH_ONLY_FALLBACK_REASON
    elif hg_low:
        missing, resolution, reason = [], 'tentative_species_level',r.HG_ONLY_FALLBACK_REASON
    elif weak_hg:
        missing, resolution, reason = ['hg'],'class_level','missing_required_hg'
    elif missing:
        resolution, reason = 'class_level','missing_required_'+','.join(missing)
    elif plan['lyso'] and not plan['fah']:
        if resolution=='chain_level':
            reason = ''
        elif plan['single_species']:
            resolution, reason = 'species_level','lyso_hg_only_fallback'
        elif plan['hg_complete'] and hg_count>=plan['hg_required']:
            resolution, reason = 'chain_level',''
        else:
            resolution, reason = 'class_level','lyso_hg_only_fallback'
    if not missing and chain_decision is not None and chain_decision[1]:
        resolution, reason = chain_decision
        if reason==r.PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON:
            missing = ['chain_confirmation']
    return CandidateScore(record,round(r._total_score_from_pool_scores(pools),4),not missing,missing,ppm,
        resolution,list(matches),pools,*r._matched_intensity_sums(spectrum,matches),reason)


def supports_index_scoring(record):
    if type(record) is not LibraryRecord:
        return False
    if reference._uses_lcb_other_support(record):
        return False
    key = reference._normalized_class_key(record.compound_class)
    return ((key in CLASSES and record.adduct in ADDUCTS)
            or (key in {'TG','TGO','TGEST'} and str(record.adduct).endswith('+')))


@with_record_facts
def score_indexes(spectrum, record, rule, indexes, tolerance_da, tolerance_ppm):
    """Score numeric occurrences first, materializing only retained evidence."""
    r = reference
    if r._uses_lcb_other_support(record):
        return NotImplemented
    if (r._missing_required_groups is not r._PRODUCTION_MISSING_GROUPS
            or r._determine_resolution is not r._PRODUCTION_RESOLUTION
            or r._non_precursor_quality_overrides is not r._PRODUCTION_QUALITY_OVERRIDES
            or r._normalized_class_key(record.compound_class) in r.NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS):
        return NotImplemented
    query = get_native_query(spectrum, exclusion=r.PRECURSOR_CLUSTER_EXCLUSION_DA)
    bridge = load_bridge()
    if (query is None or bridge is None or not hasattr(bridge,'score_pool_indexes')
            or not hasattr(bridge,'pool_metrics_indexes')
            or not hasattr(bridge,'indexed_scoring_compatible')
            or not bridge.indexed_scoring_compatible(native_engine.FragmentMatch)):
        return NotImplemented
    rule = r.single_chain_rule(record, rule)
    plan = _prepare(record, rule)
    if plan is None or not plan.get('indexed') or plan.get('native') is None:
        return NotImplemented
    ppm = ((spectrum.precursor_mz-record.precursor_mz)/record.precursor_mz)*1e6
    if (abs(spectrum.precursor_mz-record.precursor_mz)>tolerance_da if tolerance_da is not None else abs(ppm)>tolerance_ppm):
        return NotImplemented
    if plan.get('glyceride'):
        native = bridge.pool_metrics_indexes(indexes,plan['native'][0],query)
        if native is None:
            return NotImplemented
        if native is False:
            return None
        (metrics, summary), sums = native
        if not all(summary[0] & group for group in plan['groups']):
            return None
        matches = native_engine.materialize_matches(spectrum,record,indexes)
        resolution,reason = (('chain_level','') if plan['implied_chain_resolution'] else
                            r._determine_resolution(record,matches,rule,plan['derived'][0],True,plan['derived'][5]))
        missing = []
    else:
        if ('native_policy_signature' not in plan or
                _decision_signature(plan)!=plan['native_policy_signature']):
            return NotImplemented
        native = bridge.score_pool_indexes(indexes,plan['native'][0],query)
        if native is None:
            return NotImplemented
        if native is False:
            return None
        (metrics,mask,resolution_code,reason_code), sums = native
        order = ((1,'hg'),(2,'loss'),(4,'fah')) if plan['positive'] else ((4,'fah'),(1,'hg'),(2,'loss'))
        missing = ['chain_confirmation'] if mask & 8 else [name for bit,name in order if mask & bit]
        resolution,reason = _NATIVE_RESOLUTIONS[resolution_code],_decision_reason(reason_code,missing)
        matches = native_engine.materialize_matches(spectrum,record,indexes)
    pools = {item[0]: r.PoolScore(item[0],values[0],item[1],values[1],values[2],0.,values[3])
             for item,values in zip(plan['native'][1],metrics)}
    return CandidateScore(record,round(r._total_score_from_pool_scores(pools),4),not missing,
                          missing,ppm,resolution,matches,pools,*sums,reason)
