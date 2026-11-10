"""Compact, conservative search plans for immutable on-disk library records.

Plans encode necessary evidence as ORs of ANDs of occurrence-count tests.
They never award a score, choose a candidate, or replace a final evidence gate.
Unusual records and custom policies retain the reference scoring path.
"""
from __future__ import annotations

from array import array
import math
import pickle
import re
import zlib

from . import safe_pickle
from .chain_utils import extract_fragment_chain_token
from .esterified_ceramide import is_esterified_ceramide
from .negative_lpe import is_negative_lpe
from .negative_lpc import is_negative_lpc
from .positive_lpc_ether import is_positive_ether_lpc
from .evidence_policy import (
    LOSS_FRAGMENT_TYPES, POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
    _fragment_counts_as_effective_loss, _fragment_counts_as_hg,
    _is_positive_adduct, _is_positive_glyceride_chain_gate_fragment,
    _normalized_class_key, _record_expected_fah_tokens,
    _requires_precursor_fragment, _scoring_pools_for_fragment,
    sphingolipid_fragment_types,
)
from .gate_policy import (
    FAH_ONLY_FALLBACK_BLOCKED_CLASSES, expected_glyceride_substituent_groups,
    glyceride_substituent_evidence_group, is_positive_three_substituent_glyceride,
)
from .positive_pc_sodium import is_positive_pc_sodium, HG_NAMES as PC_HG_NAMES
from .positive_pe_sodium import is_positive_pe_sodium, HG_NAMES as PE_HG_NAMES
from .oxidized_chain_evidence import oxidized_chain_requirements, oxygen_count
from .resolution_policy import fragment_is_chain_evidence
from .scoring_policy import FAH_ONLY_FALLBACK_MIN_MATCHES, SPHINGO_HG_ONLY_FALLBACK_CLASS_KEYS
from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK, _series_key


PLAN_VERSION = 1
# One unconstrained alternative accepts everything; no alternatives rejects all.
ACCEPT = ((),)
REJECT = ()


def _mask(record, predicate):
    return sum(1 << i for i, fragment in enumerate(record.fragments) if predicate(fragment))


def _and(program, mask, count=1):
    if count <= 0:
        return program
    if mask.bit_count() < count:
        return REJECT
    clause = (mask, count)
    return tuple(arm if clause in arm else (*arm, clause) for arm in program)


def _any(program, alternatives):
    result = tuple(arm for mask, count in alternatives for arm in _and(program, mask, count))
    return tuple(dict.fromkeys(result))


def _sphingo_program(record, rule):
    """The standard gate OR necessary conditions for either existing fallback."""
    fragments = record.fragments
    types = [sphingolipid_fragment_types(record, fragment) for fragment in fragments]
    library_types = set().union(*types)
    def typed_mask(group):
        return sum(1 << i for i, fragment_types in enumerate(types) if fragment_types & group)
    def name_mask(names):
        return _mask(record, lambda f: str(f.name).strip() in names)
    standard = ACCEPT
    for name in sorted(rule.required_all):
        standard = _and(standard, name_mask({name}))
    for names in rule.required_any_groups:
        standard = _and(standard, name_mask(names))
    found = re.search(r'(?:\(|\s)([mdt])\d', str(record.lipid_chain_name or record.lipid_name), flags=re.IGNORECASE)
    series = _series_key(found.group(1).lower() if found else 'other')
    for key, names in rule.required_any_by_series.items():
        if key == series:
            standard = _and(standard, name_mask(names))
    for group in rule.required_type_any_groups:
        standard = _and(standard, typed_mask(group))
    for group, minimum in rule.required_type_count_groups:
        standard = _and(standard, typed_mask(group), minimum)
    for alternatives in rule.required_type_count_any_groups:
        standard = _any(standard, [(typed_mask(group), minimum) for group, minimum in alternatives])
    for group, fraction in rule.required_type_fraction_groups:
        mask = typed_mask(group)
        standard = _and(standard, mask, max(1, math.ceil(mask.bit_count() * float(fraction))))
    for pattern, fraction in rule.required_name_pattern_fraction_groups:
        mask = _mask(record, lambda f: re.search(pattern, str(f.name or ''), flags=re.IGNORECASE))
        standard = _and(standard, mask, max(1, math.ceil(mask.bit_count() * float(fraction))))
    explicit = bool(rule.required_type_any_groups or rule.required_type_count_groups
                    or rule.required_type_count_any_groups or rule.required_type_fraction_groups
                    or rule.required_name_pattern_fraction_groups)
    if not explicit:
        for fragment_type in ('Diagnostic_HG', 'C类碎片', 'LCB碎片'):
            if fragment_type in library_types:
                standard = _and(standard, typed_mask({fragment_type}))
    adduct = str(record.adduct or '').strip()
    key = re.sub(r'[^A-Za-z0-9]+', '', str(record.compound_class or '').upper())
    tiered = ((adduct == '[M-H]-' and key in {'CER1P', 'CERP'})
              or (adduct in {'[M+CH3COO]-', '[M+HCOO]-'} and key == 'SM'))
    if tiered:
        count = sum(f.fragment_type == 'Diagnostic_HG' for f in fragments)
        standard = _and(standard, typed_mask({'Diagnostic_HG'}), max(1, (count + 1) // 2)) if count else REJECT
    # Omitting identity and intensity checks is conservative: the full scorer
    # still verifies them. Only impossible alternatives disappear here.
    result = standard
    if (rule.allow_hg_only_fallback and key in SPHINGO_HG_ONLY_FALLBACK_CLASS_KEYS and not tiered
            and 'Diagnostic_HG' in library_types):
        hg = _mask(record, lambda f: f.fragment_type == 'Diagnostic_HG')
        result += _and(_and(ACCEPT, hg), (1 << len(fragments)) - 1, 2)
    if (rule.allow_fah_only_fallback and not tiered and _is_positive_adduct(adduct)
            and key not in FAH_ONLY_FALLBACK_BLOCKED_CLASSES
            and library_types & {'Diagnostic_HG', 'C类碎片'}):
        fah = _mask(record, lambda f: 'fah' in _scoring_pools_for_fragment(record, f, matched=True))
        result += _and(ACCEPT, fah, FAH_ONLY_FALLBACK_MIN_MATCHES)
    return tuple(dict.fromkeys(result))


def evidence_program(record):
    """Compile only checks that every emitted result must satisfy."""
    if len(record.fragments) > 64:
        return ACCEPT
    if is_negative_lpc(record) or is_positive_ether_lpc(record):
        return _and(ACCEPT, _mask(record, lambda f: _fragment_counts_as_hg(record, f)))
    if is_negative_lpe(record):
        precursor = _mask(record, lambda f: f.fragment_type == 'Precursor Ion'
                          and abs(f.mz-record.precursor_mz) <= 1e-6)
        hg = _mask(record, lambda f: _fragment_counts_as_hg(record, f))
        return _and(_and(ACCEPT, precursor), hg)
    key = f'{record.compound_class}_{record.adduct}'
    if is_esterified_ceramide(record.compound_class, record.lipid_chain_name):
        key = f'Cer-esterified_{record.adduct}'
    if key in SPHINGOLIPID_RULEBOOK:
        return _sphingo_program(record, SPHINGOLIPID_RULEBOOK[key])
    class_key = _normalized_class_key(record.compound_class)
    if class_key in {'FA', 'NAE'}:
        return ACCEPT
    program = ACCEPT
    positive = _is_positive_adduct(record.adduct)
    expected_fah = _record_expected_fah_tokens(record)
    if not positive and oxidized_chain_requirements(record) is None:
        for token in expected_fah:
            program = _and(program, _mask(record, lambda f: f.fragment_type == 'Diagnostic_FA'
                                         and extract_fragment_chain_token(f.name) == token))
    oxidized = oxidized_chain_requirements(record)
    if oxidized is not None:
        # Follow chain_evidence_groups: a chain-labelled loss can replace an
        # absent RCOO ion. Oxygenated chains also need an undehydrated anchor.
        def chain_fragment(fragment):
            return fragment.fragment_type in {'Diagnostic_FA', 'Diagnostic_FA_Loss', 'FA_Frag'} or (
                fragment.fragment_type == 'Neutral_Loss'
                and re.search(r'^\[M.*(?:R\d?(?:COOH|=O)|ROOH)', fragment.name))
        for token, minimum in oxidized:
            mask = _mask(record, lambda f: chain_fragment(f) and extract_fragment_chain_token(f) == token)
            program = _and(program, mask, minimum)
            if oxygen_count(token) > 0:
                anchors = _mask(record, lambda f: chain_fragment(f) and extract_fragment_chain_token(f) == token
                                and re.search(r'-(?:\d*)H2O', f.name.replace(' ', ''), re.IGNORECASE) is None)
                program = _and(program, anchors)
    hg = _mask(record, lambda f: _fragment_counts_as_hg(record, f))
    if is_positive_pe_sodium(record):
        for name in sorted(PE_HG_NAMES):
            program = _and(program, _mask(record, lambda f: f.fragment_type == 'Diagnostic_HG' and f.name == name))
    if is_positive_pc_sodium(record):
        program = _and(program, _mask(record, lambda f: f.fragment_type == 'Diagnostic_HG' and f.name in PC_HG_NAMES), 2)
        program = _and(program, _mask(record, lambda f: fragment_is_chain_evidence(record, f)), 2)
    if positive and class_key == 'PEP':
        program = _and(program, hg)
        chains = _mask(record, lambda f: fragment_is_chain_evidence(record, f))
        count = len({float(f.mz) for f in record.fragments if fragment_is_chain_evidence(record, f)})
        program = _and(program, chains, max(1, (count + 1) // 2))
    if is_positive_three_substituent_glyceride(record):
        groups = expected_glyceride_substituent_groups(record)
        if not groups:
            return REJECT
        for group in sorted(groups):
            program = _and(program, _mask(record, lambda f: glyceride_substituent_evidence_group(record, f) == group))
    # A support-only match cannot pass a declared structural gate. Preserve
    # unusual records with no such gate, and precursor-mandatory exceptions.
    if class_key != 'OXFA' and not _requires_precursor_fragment(record):
        loss = _mask(record, lambda f: _fragment_counts_as_effective_loss(record, f))
        if hg or expected_fah or (positive and loss):
            nonsupport = _mask(record, lambda f: (
                _fragment_counts_as_hg(record, f)
                or f.fragment_type not in {'Common', 'Precursor Ion', 'FA_Frag'}
                or (f.fragment_type == 'FA_Frag'
                    and (_is_positive_glyceride_chain_gate_fragment(record, f)
                         or _fragment_counts_as_effective_loss(record, f)))
            ))
            program = _and(program, nonsupport)
    return program


class SearchPlanWriter:
    """Write numeric evidence columns alongside the authoritative record blocks."""
    def __init__(self, connection):
        self.connection = connection
        self.programs, self.program_ids = [], {}
        self.supported = True
        from .equivalent_candidates import EquivalentCandidateWriter
        self.equivalent = EquivalentCandidateWriter(connection)
        connection.execute('CREATE TABLE search_plans (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL)')

    def add_block(self, block_id, records):
        self.equivalent.add_block(block_id, records)
        masses, flags, counts, record_flags, programs = array('d'), array('B'), array('I'), array('B'), array('I')
        for record in records:
            for f in record.fragments:
                if type(f.mz) not in (float, int) or not math.isfinite(f.mz) or float(f.mz) != f.mz:
                    self.supported = False
                    return
                masses.append(f.mz)
                flags.append(int(f.fragment_type == 'FA_Frag' or f.fragment_type in LOSS_FRAGMENT_TYPES)
                             | (2 * int(f.fragment_type == 'Precursor Ion')))
            counts.append(len(record.fragments))
            class_key = re.sub(r'[^A-Za-z0-9]+', '', str(record.compound_class or '').upper())
            positive = str(record.adduct or '').strip().endswith('+')
            # NAE re-matches a synthetic precursor in the full scorer; never
            # replace that match by a library occurrence, even for custom data.
            mode = int(positive and class_key in POSITIVE_GLYCERIDE_RCO_GATE_CLASSES)
            mode |= 2 * int(str(record.compound_class or '').strip().upper() == 'FA'
                            and f'{record.compound_class}_{record.adduct}' not in SPHINGOLIPID_RULEBOOK)
            mode |= 4 * int(class_key == 'NAE')
            record_flags.append(mode)
            program = evidence_program(record)
            plan_id = self.program_ids.get(program)
            if plan_id is None:
                self.program_ids[program] = plan_id = len(self.programs)
                self.programs.append(program)
            programs.append(plan_id)
        payload = zlib.compress(pickle.dumps((masses, flags, counts, record_flags, programs), protocol=5), level=6)
        self.connection.execute('INSERT INTO search_plans VALUES (?, ?)', (block_id, payload))

    def finish(self, metadata):
        self.equivalent.finish(metadata)
        if not self.supported:
            self.connection.execute('DROP TABLE search_plans')
            metadata.pop('search_plan_version', None)
            return
        plan_offsets, arm_offsets, masks, counts = array('I', [0]), array('I', [0]), array('Q'), array('B')
        for program in self.programs:
            for arm in program:
                for mask, count in arm:
                    masks.append(mask); counts.append(count)
                arm_offsets.append(len(masks))
            plan_offsets.append(len(arm_offsets) - 1)
        payload = zlib.compress(pickle.dumps((plan_offsets, arm_offsets, masks, counts), protocol=5), level=6)
        self.connection.execute('INSERT INTO search_plans VALUES (-1, ?)', (payload,))
        metadata['search_plan_version'] = PLAN_VERSION


def read_programs(connection):
    row = connection.execute('SELECT payload FROM search_plans WHERE block_id=-1').fetchone()
    if row is None:
        raise ValueError('Missing search evidence programs')
    return safe_pickle.loads(zlib.decompress(row[0]))
