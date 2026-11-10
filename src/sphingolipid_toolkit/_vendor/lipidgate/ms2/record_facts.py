"""Cache library-only policy facts, with bounded canonical-search reuse.

Public scoring calls get a fresh scope. The canonical indexed search may retain
immutable library facts after checking every input field for edits. Experimental
matches, intensities, scores, gates and output rankings are never cached.
"""
from __future__ import annotations

from contextvars import ContextVar
from collections import OrderedDict
from dataclasses import fields, is_dataclass
from functools import wraps
import sys
from . import native_bridge


def record_revision(record):
    """Immutable input fields used by the scoring policy, including occurrences."""
    return (record.record_id, record.compound_class, record.lipid_name,
            record.lipid_chain_name, type(record.precursor_mz), record.precursor_mz,
            record.adduct, record.formula, record.polarity,
            tuple((id(f), type(f.mz), f.mz, f.name, f.fragment_type, type(f.intensity),
                   f.intensity, type(f.weight), f.weight, f.required_group)
                  for f in record.fragments), tuple(record.metadata.items()))


def record_unchanged(record, revision):
    kernel = native_bridge.load_bridge()
    if kernel is not None and hasattr(kernel, "record_unchanged"):
        answer = kernel.record_unchanged(record, revision)
        if answer is not None:
            return answer
    return record_revision(record) == revision


_current_facts = ContextVar("lipidgate_record_facts", default=None)
_current_library_facts = ContextVar("lipidgate_library_facts", default=None)
_current_spectrum_facts = ContextVar("lipidgate_spectrum_facts", default=None)
_unspecified = object()
MAX_LIBRARY_FACT_RECORDS = 4096
MAX_LIBRARY_FACT_BYTES = 64 * 1024 * 1024


def _fact_bytes(value, seen):
    identity = id(value)
    if identity in seen or callable(value) or isinstance(value, type):
        return 0
    seen.add(identity)
    size = sys.getsizeof(value)
    if isinstance(value, dict):
        return size + sum(_fact_bytes(k, seen) + _fact_bytes(v, seen) for k, v in value.items())
    if isinstance(value, (tuple, list, set, frozenset)):
        return size + sum(_fact_bytes(v, seen) for v in value)
    if is_dataclass(value):
        return size + sum(_fact_bytes(getattr(value, f.name), seen) for f in fields(value))
    return size


class LibraryFactCache:
    """An LRU of source-record facts; no query/result objects are retained."""
    def __init__(self):
        self.entries = OrderedDict()
        self.bytes = 0

    def values(self, record):
        key = id(record)
        entry = self.entries.get(key)
        if entry is not None and entry[0] is record and record_unchanged(record, entry[1]):
            self.entries.move_to_end(key)
            return entry[2]
        revision = record_revision(record)
        if entry is not None:
            self.bytes -= entry[3]
        values = {}
        self.entries[key] = [record, revision, values, 0, -1]
        self.entries.move_to_end(key)
        return values

    def trim(self, record, values):
        entry = self.entries.get(id(record))
        if entry is None or entry[2] is not values:
            return
        if len(values) != entry[4]:
            self.bytes -= entry[3]
            entry[3] = _fact_bytes((record, entry[1], values), set()) + 256
            entry[4] = len(values)
            self.bytes += entry[3]
        while len(self.entries) > MAX_LIBRARY_FACT_RECORDS or self.bytes > MAX_LIBRARY_FACT_BYTES:
            _, removed = self.entries.popitem(last=False)
            self.bytes -= removed[3]

    def clear(self):
        self.entries.clear()
        self.bytes = 0


def with_search_facts(function):
    """Let an explicitly prepared searcher reuse only library facts."""
    @wraps(function)
    def search(self, *args, **kwargs):
        cache = getattr(self, "_library_fact_cache", None)
        query_token = _current_spectrum_facts.set(None)
        if cache is None:
            try:
                return function(self, *args, **kwargs)
            finally:
                _current_spectrum_facts.reset(query_token)
        token = _current_library_facts.set(cache)
        try:
            return function(self, *args, **kwargs)
        finally:
            _current_library_facts.reset(token)
            _current_spectrum_facts.reset(query_token)
    return search


def enable_spectrum_facts(spectrum):
    """Enable one invocation's read-only peak facts for canonical search only."""
    _current_spectrum_facts.set((spectrum, {}))


_native_query_key = object()


def set_native_query(spectrum, query, bridge, exclusion):
    """Retain one numeric snapshot until this canonical search returns."""
    current = _current_spectrum_facts.get()
    if current is not None and current[0] is spectrum:
        current[1][_native_query_key] = (query, bridge, exclusion)
        # Query-frame admission validates the exact builtin mass/raw fields.
        current[1][_safe_spectrum_key] = True


def get_native_query(spectrum=None, *, exclusion=None):
    current = _current_spectrum_facts.get()
    if current is None or (spectrum is not None and current[0] is not spectrum):
        return None
    found = current[1].get(_native_query_key)
    if (found is None or found[1] is not native_bridge.load_bridge()
            or (exclusion is not None and exclusion != found[2])):
        return None
    return found[0]


def cached_spectrum_fact(function):
    @wraps(function)
    def fact(spectrum):
        current = _current_spectrum_facts.get()
        if current is None or current[0] is not spectrum:
            return function(spectrum)
        values = current[1]
        if _safe_spectrum_key not in values:
            kernel = native_bridge.load_bridge()
            safe = (kernel.spectrum_is_builtin(spectrum)
                    if kernel is not None and hasattr(kernel, 'spectrum_is_builtin') else None)
            from .models import ExperimentalPeak, ExperimentalSpectrum
            values[_safe_spectrum_key] = safe if safe is not None else (
                type(spectrum) is ExperimentalSpectrum
                and type(spectrum.precursor_mz) in (float, int)
                and type(spectrum.peaks) in (list, tuple)
                and all(type(p) is ExperimentalPeak and type(p.mz) in (float, int)
                        and type(p.intensity) in (float, int) for p in spectrum.peaks))
        if not values[_safe_spectrum_key]:
            return function(spectrum)
        if function not in values:
            values[function] = function(spectrum)
        return values[function]
    return fact


_safe_spectrum_key = object()


def _rule_revision(rule):
    # ClassRule and ScoreProfile are frozen, but their dict fields are mutable.
    profile = rule.score_profile
    return (rule.lipid_class, tuple(rule.required_groups.items()), rule.chain_level_min_fah,
            rule.allow_hg_only_if_no_fah, rule.positive_hg_min_matches_if_no_fah,
            rule.positive_hg_min_fraction, rule.positive_hg_min_matches,
            rule.positive_loss_min_fraction, rule.positive_loss_min_matches,
            rule.positive_hg_complete_can_resolve_chain,
            rule.positive_signature_min_matches_if_no_fah, rule.allow_loss_only_if_no_fah,
            rule.require_loss_with_fah_only, rule.negative_hg_min_matches,
            rule.negative_require_all_fah, tuple(profile.pool_weights.items()),
            tuple(profile.metric_weights.items()), profile.single_group_fallback_min_matches,
            profile.single_group_fallback_min_relative_intensity_sum)


def cached_rule_fact(function):
    """A library fact parameterized by the exact, revision-checked class rule."""
    last_key = object()
    @wraps(function)
    def fact(record, rule):
        current = _current_facts.get()
        if current is None or current[0] is not record:
            return function(record, rule)
        values = current[1]
        last = values.get(last_key)
        revision = None
        if last is not None and last[0] == id(rule):
            kernel = native_bridge.load_bridge()
            unchanged = (kernel.rule_unchanged(rule, last[1])
                         if kernel is not None and hasattr(kernel, 'rule_unchanged') else None)
            if unchanged is True:
                return last[2]
            if unchanged is None:
                revision = _rule_revision(rule)
                if revision == last[1]:
                    return last[2]
        if revision is None:
            revision = _rule_revision(rule)
        key = (function, revision)
        if key not in values:
            values[key] = function(record, rule)
        # The numeric identity is only a lookup hint. Even if Python reuses an
        # ID, every field/dict is verified against the retained revision first.
        values[last_key] = (id(rule), revision, values[key])
        return values[key]
    return fact


def cached_record_fact(function):
    """Cache a read-only bool/int fact about the active candidate record."""
    @wraps(function)
    def fact(record):
        current = _current_facts.get()
        if current is None or current[0] is not record:
            return function(record)
        values = current[1]
        if function not in values:
            values[function] = function(record)
        return values[function]
    return fact


def cached_match_fact(function):
    """Reuse a decision only within this score, on an immutable match tuple.

    This separate scratch dictionary is discarded before the next candidate.
    It never enters the persistent library facts. Mutable public match lists
    retain their reference behavior.
    """
    @wraps(function)
    def fact(record, matches, *args, **kwargs):
        current = _current_facts.get()
        if current is None or current[0] is not record or type(matches) is not tuple:
            return function(record, matches, *args, **kwargs)
        extra = tuple(tuple(v) if isinstance(v, list) else v for v in args)
        key = (function, id(matches), extra, tuple(kwargs.items()))
        values = current[3]
        try:
            found = values.get(key, _unspecified)
        except TypeError:
            return function(record, matches, *args, **kwargs)
        if found is not _unspecified:
            return found[1]
        result = function(record, matches, *args, **kwargs)
        values[key] = (matches, result)
        return result
    return fact


def with_record_facts(function):
    """Give each candidate invocation a fresh cache; restore nested callers."""
    @wraps(function)
    def score(spectrum, record, *args, **kwargs):
        adduct = getattr(record, "adduct", None)
        # Negative spectra usually have short, cheap routing branches: tuple
        # cache keys cost more than recomputing them. Record facts still cache.
        fragment_cache = adduct is None or str(adduct).strip().endswith("+")
        cache = _current_library_facts.get()
        values = {} if cache is None else cache.values(record)
        token = _current_facts.set((record, values, fragment_cache, {}))
        try:
            return function(spectrum, record, *args, **kwargs)
        finally:
            _current_facts.reset(token)
            if cache is not None:
                cache.trim(record, values)
    return score


def with_result_facts(function):
    """Reuse the same library scope when reranking a scored candidate."""
    @with_record_facts
    def invoke(spectrum, record, self, result, *args, **kwargs):
        return function(self, spectrum, result, *args, **kwargs)

    @wraps(function)
    def rerank(self, spectrum, result, *args, **kwargs):
        return invoke(spectrum, result.record, self, result, *args, **kwargs)
    return rerank


def with_output_facts(function):
    """Reuse validated library facts in positive secondary-result selection.

    Negative and standalone calls keep their cheaper reference route. Query
    matches and decisions remain fresh; only library fields enter the LRU.
    """
    @with_record_facts
    def invoke(spectrum, record, self, result):
        return function(self, result)

    @wraps(function)
    def select(self, result):
        if (_current_library_facts.get() is None
                or not str(result.record.adduct or '').strip().endswith('+')):
            return function(self, result)
        return invoke(None, result.record, self, result)
    return select


def cached_fragment_fact(function=None, *, matched_independent=False):
    """Reuse immutable fragment routing facts within the active candidate only.

    The matched/unmatched flag is part of the key unless a helper explicitly
    declares that its classification does not depend on it. Fragment identities
    remain distinct, including repeated occurrences with identical masses/labels.
    """
    if function is None:
        return lambda wrapped: cached_fragment_fact(
            wrapped, matched_independent=matched_independent,
        )

    @wraps(function)
    def fact(record, fragment, *, matched=_unspecified):
        current = _current_facts.get()
        if current is not None and current[0] is record and (current[2] or matched_independent):
            key = (function, id(fragment), _unspecified if matched_independent else matched)
            result = current[1].get(key, _unspecified)
            if result is not _unspecified:
                return result
        else:
            key = None
        result = (function(record, fragment) if matched is _unspecified else
                  function(record, fragment, matched=matched))
        if key is not None:
            current[1][key] = result
        return result
    return fact
