"""Skip duplicate evaluations that the reporting policy provably collapses.

Supported protonated phospholipids and d-series SM without chain hits qualify. A profile
records exact pool counts, gates and occurrence roles. The surviving record is
the same ppm-preferred representative chosen by the final reporting collapse;
it still passes through the original full scorer and exporter.
"""
from array import array
from collections import Counter
import bisect
import math
import pickle
import re
import sqlite3
import zlib

from . import safe_pickle
from . import native_bridge
from .chain_utils import extract_chain_tokens
from .evidence_policy import (_fragment_counts_as_hg, _scoring_pools_for_fragment,
                              _record_expected_fah_tokens, _record_loss_fragment_count)
from .gate_policy import required_positive_hg_hits
from .record_facts import _rule_revision, with_record_facts
from .resolution_policy import get_chain_evidence_profile, required_chain_evidence_hits
from .rules import DEFAULT_RULES
from .scoring import _pool_weights_for_record
from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK


VERSION = 2
CLASSES = ('PC', 'PC-O', 'PC-P', 'PE', 'PE-O', 'PG', 'PG-O', 'PI', 'PI-O', 'PA', 'PA-O', 'PS', 'PS-O')
SM_ADDUCTS = ('[M+H]+', '[M+Na]+', '[M+HCOO]-', '[M+CH3COO]-')
POOL_BITS = {'fah': 1, 'hg': 2, 'other': 4, 'lcb': 8}
INELIGIBLE_MATCH = 64


def _sphingo_revision(rule):
    return (rule.lipid_class, rule.adduct, frozenset(rule.required_all),
            tuple(frozenset(s) for s in rule.required_any_groups),
            tuple((k, frozenset(v)) for k, v in rule.required_any_by_series.items()),
            tuple(frozenset(s) for s in rule.required_type_any_groups),
            tuple((frozenset(s), n) for s, n in rule.required_type_count_groups),
            tuple(tuple((frozenset(s), n) for s, n in group) for group in rule.required_type_count_any_groups),
            tuple((frozenset(s), n) for s, n in rule.required_type_fraction_groups),
            tuple(rule.required_name_pattern_fraction_groups),
            rule.only_non_hydroxy_fa, rule.prefer_d_series_if_ambiguous,
            rule.allow_hg_only_fallback, rule.allow_fah_only_fallback)


def policy_revision(key):
    if key in CLASSES:
        return _rule_revision(DEFAULT_RULES.get(key))
    if key.startswith('SM_') and key in SPHINGOLIPID_RULEBOOK:
        return _rule_revision(DEFAULT_RULES.get('SM')), _sphingo_revision(SPHINGOLIPID_RULEBOOK[key])
    return None


_REFERENCE_POLICY_REVISION = policy_revision
_REFERENCE_RULE_REVISION = _rule_revision
_REFERENCE_SPHINGO_REVISION = _sphingo_revision


def policy_unchanged(key, revision):
    if (policy_revision is _REFERENCE_POLICY_REVISION and _rule_revision is _REFERENCE_RULE_REVISION
            and (key in CLASSES or key.startswith('SM_') and key in SPHINGOLIPID_RULEBOOK)):
        bridge = native_bridge.load_bridge()
        if bridge is not None and hasattr(bridge, 'rule_unchanged'):
            if key in CLASSES:
                unchanged = bridge.rule_unchanged(DEFAULT_RULES.get(key), revision)
            elif (_sphingo_revision is _REFERENCE_SPHINGO_REVISION
                  and hasattr(bridge,'sphingo_rule_unchanged')
                  and type(revision) is tuple and len(revision)==2):
                unchanged = bridge.rule_unchanged(DEFAULT_RULES.get('SM'),revision[0])
                if unchanged is True:
                    unchanged=bridge.sphingo_rule_unchanged(SPHINGOLIPID_RULEBOOK[key],revision[1])
            else:
                unchanged = None
            if unchanged is not None:
                return unchanged
    return policy_revision(key) == revision


def rule_revisions():
    return {key: policy_revision(key) for key in (*CLASSES, *('SM_'+a for a in SM_ADDUCTS))}


def _sm_profile(record):
    # d/non-hydroxy records do not participate in d/h versus t isomer ranking.
    # Hydroxy/trihydroxy candidates remain fully scored and separately ranked.
    if (record.adduct not in SM_ADDUCTS or
            re.fullmatch(r'SM\(d\d+:\d+/\d+:\d+\)', record.lipid_chain_name) is None):
        return None
    roles, allowed_labels = [], []
    totals = [0, 0, 0, 0]
    hg_count = 0
    for i, fragment in enumerate(record.fragments):
        pools = _scoring_pools_for_fragment(record, fragment, matched=False)
        role = sum(POOL_BITS[name] for name in pools)
        for index, name in enumerate(POOL_BITS):
            totals[index] += int(name in pools)
        hg = fragment.fragment_type == 'Diagnostic_HG'
        hg_count += int(hg)
        role |= 16 * int(hg)
        role |= 32 * int(fragment.fragment_type == 'Precursor Ion')
        # All chain roles, including experimental aliases, keep the full path.
        if fragment.fragment_type not in {'Diagnostic_HG', 'Common', 'Precursor Ion', 'C类碎片', 'Neutral_Loss'}:
            role |= INELIGIBLE_MATCH
        else:
            allowed_labels.append((i, fragment.name, fragment.fragment_type))
        roles.append(role)
    if not hg_count:
        return None
    rule = DEFAULT_RULES.get('SM')
    key = ('SM', record.adduct, record.lipid_name, tuple(totals),
           tuple(_pool_weights_for_record(record, rule).items()),
           tuple(sorted(Counter(f.fragment_type for f in record.fragments).items())),
           tuple(allowed_labels), 'd')
    return key, bytes(roles), 'SM_'+record.adduct


@with_record_facts
def _profile(spectrum, record):
    if record.compound_class == 'SM':
        if (type(record.precursor_mz) not in (float, int)
                or not math.isfinite(record.precursor_mz) or record.precursor_mz <= 0):
            return None
        return _sm_profile(record)
    if (record.compound_class not in CLASSES or record.adduct != '[M+H]+'
            or type(record.precursor_mz) not in (float, int)
            or not math.isfinite(record.precursor_mz) or record.precursor_mz <= 0):
        return None
    rule = DEFAULT_RULES.get(record.compound_class)
    if rule.positive_loss_min_fraction or rule.positive_loss_min_matches:
        return None
    profile = get_chain_evidence_profile(record)
    configured_species = (profile is not None and profile.species_fallback_without_chain
                          and (required_chain_evidence_hits(record, profile) > 0
                               or profile.applies_without_library_evidence))
    if (len(extract_chain_tokens(record.lipid_chain_name)) < 2 or not (
            configured_species or _record_expected_fah_tokens(record) or _record_loss_fragment_count(record) > 0)):
        return None
    roles = []
    totals = [0, 0, 0, 0]
    hg_count = 0
    for fragment in record.fragments:
        pools = _scoring_pools_for_fragment(record, fragment, matched=False)
        role = sum(POOL_BITS[name] for name in pools)
        for index, name in enumerate(POOL_BITS):
            if name in pools:
                totals[index] += 1
        hg = _fragment_counts_as_hg(record, fragment)
        hg_count += int(hg)
        role |= 16 * int(hg)
        role |= 32 * int(fragment.fragment_type == 'Precursor Ion')
        if fragment.fragment_type not in {'Diagnostic_HG', 'Candidate_HG', 'Common', 'Precursor Ion'}:
            role |= INELIGIBLE_MATCH
        roles.append(role)
    if not hg_count:
        return None
    key = (record.compound_class, record.adduct, record.lipid_name, tuple(totals),
           tuple(_pool_weights_for_record(record, rule).items()),
           hg_count, required_positive_hg_hits(rule, _record_expected_fah_tokens(record), hg_count))
    return key, bytes(roles), record.compound_class


class EquivalentCandidateWriter:
    def __init__(self, connection):
        self.connection = connection
        self.groups, self.profiles, self.profile_ids = {}, [(0, b'', '')], {}
        connection.execute('CREATE TABLE equivalent_candidates (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL)')

    def add_block(self, block_id, records):
        ids = array('I')
        for record in records:
            result = _profile(None, record)
            if result is None:
                ids.append(0)
                continue
            key, roles, policy = result
            group = self.groups.setdefault(key, len(self.groups) + 1)
            profile = (group, roles, policy)
            if profile not in self.profile_ids:
                self.profile_ids[profile] = len(self.profiles)
                self.profiles.append(profile)
            ids.append(self.profile_ids[profile])
        payload = zlib.compress(pickle.dumps(ids, protocol=5), level=6)
        self.connection.execute('INSERT INTO equivalent_candidates VALUES (?, ?)', (block_id, payload))

    def finish(self, metadata):
        header = (VERSION, rule_revisions(), self.profiles)
        payload = zlib.compress(pickle.dumps(header, protocol=5), level=6)
        self.connection.execute('INSERT INTO equivalent_candidates VALUES (-1, ?)', (payload,))
        metadata['equivalent_candidate_version'] = VERSION


class EquivalentCandidates:
    def __init__(self, library):
        self.library = library
        row = library._connection.execute('SELECT payload FROM equivalent_candidates WHERE block_id=-1').fetchone()
        if row is None:
            raise ValueError('Missing equivalence header')
        version, self.revisions, self.profiles = safe_pickle.loads(zlib.decompress(row[0]))
        if version != VERSION or not isinstance(self.revisions, dict):
            raise ValueError('Unsupported equivalence policy')
        if not self.profiles or self.profiles[0] != (0, b'', ''):
            raise ValueError('Invalid equivalence profiles')
        if any(type(group) is not int or group < 0 or type(roles) is not bytes
               or any(role > 127 for role in roles) or (group and policy not in self.revisions)
               for group, roles, policy in self.profiles):
            raise ValueError('Invalid equivalence role')
        self.ids = array('I')
        positions, size = library._positions, library.metadata['block_size']
        if len(self.profiles) == 1:
            return
        first = positions[0] // size if positions else 0
        last = positions[-1] // size if positions else -1
        for block_id, payload in library._connection.execute(
                'SELECT block_id,payload FROM equivalent_candidates WHERE block_id>=? AND block_id<=? ORDER BY block_id',
                (first, last)):
            left = bisect.bisect_left(positions, block_id * size)
            right = bisect.bisect_left(positions, (block_id + 1) * size, left)
            if left == right:
                continue
            ids = safe_pickle.loads(zlib.decompress(payload))
            if type(ids) is not array or ids.typecode != 'I' or any(i >= len(self.profiles) for i in ids):
                raise ValueError('Invalid equivalence block')
            for visible in range(left, right):
                self.ids.append(ids[positions[visible] - block_id * size])
        if len(self.ids) != len(positions):
            raise ValueError('Incomplete equivalence columns')

    @property
    def bytes(self):
        return len(self.ids) * self.ids.itemsize + sum(40 + len(roles) for _, roles, _ in self.profiles)

    def collapse(self, spectrum, matches):
        if not self.ids or type(spectrum.precursor_mz) not in (float, int):
            return matches
        keep, best, checked = set(), {}, {}
        for index, hits in matches.items():
            profile = self.ids[index]
            if not profile or hits is None:
                keep.add(index)
                continue
            group, roles, policy = self.profiles[profile]
            if policy not in checked:
                checked[policy] = policy_unchanged(policy, self.revisions.get(policy))
            if not checked[policy]:
                keep.add(index)
                continue
            if len(roles) != len(hits):
                return matches
            signature = []
            for hit, role in zip(hits, roles):
                if hit >= 0:
                    if role & INELIGIBLE_MATCH:
                        break
                    signature.append((hit, role))
            else:
                key = (group, tuple(signature))
                precursor = self.library.precursors[index]
                ppm = abs(((spectrum.precursor_mz - precursor) / precursor) * 1e6)
                prior = best.get(key)
                if prior is None:
                    best[key] = (index, index, ppm)
                    keep.add(index)
                elif ppm < prior[2]:
                    if prior[1] != prior[0]:
                        keep.remove(prior[1])
                    keep.add(index)
                    best[key] = (prior[0], index, ppm)
                continue
            keep.add(index)
        # Keep both the first and ppm-best member, in original library order.
        # The first preserves a reporting group's insertion slot. The best
        # remains at its original position, preserving ties with other profiles
        # of the same reported identity (different pool counts can score alike).
        return {index: hits for index, hits in matches.items() if index in keep}


def load_equivalent_candidates(library):
    if library.metadata.get('equivalent_candidate_version') != VERSION:
        return None
    try:
        return EquivalentCandidates(library)
    except (ValueError, TypeError, IndexError, OverflowError, pickle.UnpicklingError,
            EOFError, sqlite3.DatabaseError, zlib.error):
        return None
