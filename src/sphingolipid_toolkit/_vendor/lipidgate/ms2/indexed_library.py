"""Read mass-sorted library records from bounded, compressed disk blocks.

The authoritative blocks contain exact normalized LibraryRecord objects, in
the same stable precursor order used by the list-based searcher. Optional
numeric search columns derive necessary conditions without replacing scores.
"""

from __future__ import annotations

from array import array
from collections import OrderedDict
from collections.abc import Iterable, Sequence
import bisect
import json
from itertools import accumulate
import math
import operator
from pathlib import Path
import pickle
import re
import sqlite3
import sys
import tempfile
import zlib

from .models import FragmentRecord, LibraryRecord
from .record_facts import record_revision
from . import safe_pickle


INDEX_FORMAT_VERSION = 3
BLOCK_SIZE = 1024
MAX_CACHED_BLOCKS = 4
MAX_SELECTED_RECORDS = 4096
MAX_SELECTED_RECORD_BYTES = 32 * 1024 * 1024
SQL_CACHE_KIB = 4096
LEGACY_RECORD_ENCODING = "tuples-v1"
TEMPLATE_RECORD_ENCODING = "fragment-templates-v2"
RECORD_ENCODING = "fragment-columns-v3"


def _compact_unsigned(values):
    """Use the narrowest exact integer column; widen only at the C ABI."""
    highest = max(values, default=0)
    return array("B" if highest <= 255 else "H" if highest <= 65535 else "I", values)


def _column_offsets(columns):
    columns["offsets"] = array("I", accumulate(columns["counts"], initial=0))
    return columns


def _pack_record(record):
    return (record.record_id, record.compound_class, record.lipid_name,
            record.lipid_chain_name, record.precursor_mz, record.adduct,
            record.formula, record.polarity,
            [(f.mz, f.name, f.fragment_type, f.intensity, f.weight, f.required_group)
             for f in record.fragments], record.metadata)


def _unpack_record(value):
    return LibraryRecord(*value[:8], fragments=[FragmentRecord(*f) for f in value[8]], metadata=value[9])


def _record_cache_bytes(record):
    # Count shared strings separately for each record: this intentionally
    # overestimates ownership and keeps the selected cache within its budget.
    size = sys.getsizeof(record) + sys.getsizeof(record.fragments) + sys.getsizeof(record.metadata) + 256
    for value in (record.record_id, record.compound_class, record.lipid_name,
                  record.lipid_chain_name, record.precursor_mz, record.adduct,
                  record.formula, record.polarity):
        size += sys.getsizeof(value)
    for fragment in record.fragments:
        size += sys.getsizeof(fragment)
        for value in (fragment.mz, fragment.name, fragment.fragment_type,
                      fragment.intensity, fragment.weight, fragment.required_group):
            size += sys.getsizeof(value)
    for key, value in record.metadata.items():
        size += sys.getsizeof(key) + sys.getsizeof(value)
    return size


def _template_fragment(template, position, columns=None):
    if columns is None or columns["masses"] is None:
        return template
    mass = columns["masses"][position]
    if columns["mass_kinds"][position]:
        mass = int(mass)
    return (mass, *template)


def _unpack_template_record(value, templates, occurrences=None, columns=None):
    # Repeated fragment occurrences must remain distinct objects: scoring uses
    # fragment identity as well as mass, including repeated-chain evidence.
    return LibraryRecord(*value[:8], fragments=[FragmentRecord(*_template_fragment(templates[i], i, columns)) for i in
                         (value[8] if occurrences is None else occurrences)], metadata=value[9])


def _encode_block(records):
    """Store each exact fragment definition once, retaining every occurrence."""
    templates, template_ids, values = [], {}, []
    occurrences, offsets = array("I"), array("I", [0])
    pairs, pair_ids, record_pairs = [], {}, array("I")
    for record in records:
        value = _pack_record(record)
        record_occurrences = []
        for fragment in value[8]:
            # Byte keys also distinguish numeric types and signed zero. No
            # rounding, m/z-only merging or fragment removal is permitted.
            key = pickle.dumps(fragment, protocol=pickle.HIGHEST_PROTOCOL)
            position = template_ids.get(key)
            if position is None:
                position = len(templates)
                template_ids[key] = position
                templates.append(fragment)
            record_occurrences.append(position)
        occurrences.extend(record_occurrences)
        offsets.append(len(occurrences))
        # The local row number replaces per-record lists; the flat column
        # retains every occurrence and its exact order without duplicating it.
        values.append((*value[:8], None, value[9]))
        pair = (value[1], value[5])
        pair_key = pickle.dumps(pair, protocol=pickle.HIGHEST_PROTOCOL)
        if pair_key not in pair_ids:
            pair_ids[pair_key] = len(pairs)
            pairs.append(pair)
        record_pairs.append(pair_ids[pair_key])
    types, type_ids, fragment_types = [], {}, array("I")
    for fragment in templates:
        key = pickle.dumps(fragment[2], protocol=pickle.HIGHEST_PROTOCOL)
        if key not in type_ids:
            type_ids[key] = len(types)
            types.append(fragment[2])
        fragment_types.append(type_ids[key])
    # Original tuple masses/types remain authoritative and unchanged. Only
    # exactly representable built-in numbers receive an optional C ABI view.
    numeric = all(type(f[0]) in (float, int) and math.isfinite(f[0]) and float(f[0]) == f[0]
                  for f in templates)
    columns = dict(masses=array("d", (f[0] for f in templates)) if numeric else None,
                   mass_kinds=array("B", (int(type(f[0]) is int) for f in templates)) if numeric else None,
                   type_names=types, type_ids=_compact_unsigned(fragment_types), occurrences=_compact_unsigned(occurrences),
                   counts=_compact_unsigned([b-a for a,b in zip(offsets, offsets[1:])]),
                   class_adduct_names=pairs, class_adduct_ids=_compact_unsigned(record_pairs))
    # Store mass once in the column, with its original float/int type. The
    # remaining annotations/weights/groups stay in exact template tuples.
    stored_templates = [f[1:] for f in templates] if numeric else templates
    return pickle.dumps((stored_templates, values, columns), protocol=pickle.HIGHEST_PROTOCOL)


def _decode_block(payload, metadata):
    values = safe_pickle.loads(zlib.decompress(payload))
    if metadata.get("record_encoding") == RECORD_ENCODING:
        templates, records, columns = values
        _column_offsets(columns)
        occurrences = memoryview(columns["occurrences"])
        offsets = columns["offsets"]
        return [_unpack_template_record(value, templates, occurrences[offsets[i]:offsets[i+1]], columns)
                for i, value in enumerate(records)]
    if metadata.get("record_encoding") == TEMPLATE_RECORD_ENCODING:
        templates, records = values
        return [_unpack_template_record(value, templates) for value in records]
    if metadata.get("record_encoding") == LEGACY_RECORD_ENCODING:
        return [_unpack_record(value) for value in values]
    return values


class _RecordBlock:
    __slots__ = ("values", "records", "packed", "templates", "column_buffers", "numeric_view")

    def __init__(self, payload, metadata):
        self.values = safe_pickle.loads(zlib.decompress(payload))
        encoding = metadata.get("record_encoding")
        self.templates = None
        self.column_buffers = None
        if encoding == RECORD_ENCODING:
            self.templates, self.values, self.column_buffers = self.values
            _column_offsets(self.column_buffers)
        elif encoding == TEMPLATE_RECORD_ENCODING:
            self.templates, self.values = self.values
        self.packed = encoding in {RECORD_ENCODING, TEMPLATE_RECORD_ENCODING, LEGACY_RECORD_ENCODING}
        self.records = {}
        # Optional numeric view shares the same four-block lifetime. It never
        # retains a LibraryRecord or changes its fields/fragment identities.
        self.numeric_view = None

    def record(self, offset):
        if not self.packed:
            return self.values[offset]
        result = self.records.get(offset)
        if result is None:
            result = (_unpack_template_record(self.values[offset], self.templates, self.occurrence_ids(offset), self.column_buffers)
                      if self.templates is not None else _unpack_record(self.values[offset]))
            self.records[offset] = result
        return result

    def occurrence_ids(self, offset):
        if self.column_buffers is None:
            return self.values[offset][8]
        columns = self.column_buffers
        return memoryview(columns["occurrences"])[columns["offsets"][offset]:columns["offsets"][offset+1]]

    def fragment_template(self, position):
        return _template_fragment(self.templates[position], position, self.column_buffers)


def _readonly_connection(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
    connection.execute(f"PRAGMA cache_size=-{SQL_CACHE_KIB}")
    connection.execute("PRAGMA mmap_size=0")
    return connection


def _read_metadata(connection: sqlite3.Connection) -> dict:
    row = connection.execute("SELECT value FROM metadata WHERE key='library'").fetchone()
    if row is None:
        raise ValueError("Library index is missing metadata")
    metadata = json.loads(row[0])
    supported = {1: {None, LEGACY_RECORD_ENCODING}, 2: {TEMPLATE_RECORD_ENCODING},
                 INDEX_FORMAT_VERSION: {RECORD_ENCODING}}
    if metadata.get("index_format") not in supported:
        raise ValueError("Library index format is incompatible")
    if metadata.get("record_encoding") not in supported[metadata["index_format"]]:
        raise ValueError("Library record encoding is incompatible")
    return metadata


def validate_precursor_range(lower=None, upper=None):
    values = []
    for value in (lower, upper):
        if value is None:
            values.append(None)
        else:
            value = float(value)
            if not math.isfinite(value) or value <= 0:
                raise ValueError("前体 m/z 范围必须是大于 0 的有限数值")
            values.append(value)
    lower, upper = values
    if lower is not None and upper is not None and lower > upper:
        raise ValueError("前体 m/z 下限不能大于上限")
    return lower, upper


def write_library_index(path: Path, records: Sequence[LibraryRecord], identity: dict) -> dict:
    """Build an atomic, verified index at packaging time, never at GUI startup."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda record: record.precursor_mz)
    metadata = dict(identity, index_format=INDEX_FORMAT_VERSION, block_size=BLOCK_SIZE,
                    record_encoding=RECORD_ENCODING,
                    record_count=len(ordered),
                    classes=sorted({r.compound_class for r in ordered}),
                    adducts=sorted({r.adduct for r in ordered}))
    class_ids = {value: index for index, value in enumerate(metadata["classes"])}
    adduct_ids = {value: index for index, value in enumerate(metadata["adducts"])}
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.stem, suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    connection = None
    try:
        connection = sqlite3.connect(temporary)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute(f"PRAGMA cache_size=-{SQL_CACHE_KIB}")
        connection.executescript(
            "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);"
            "CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, precursor REAL NOT NULL, "
            "class_id INTEGER NOT NULL, adduct_id INTEGER NOT NULL);"
            "CREATE TABLE blocks (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL);"
        )
        from .search_plan import SearchPlanWriter
        search_plans = SearchPlanWriter(connection)
        for start in range(0, len(ordered), BLOCK_SIZE):
            block = ordered[start:start + BLOCK_SIZE]
            payload = zlib.compress(_encode_block(block), level=6)
            # Every normalized field and stable record order must survive storage.
            if _decode_block(payload, metadata) != block:
                raise ValueError(f"Library block {start // BLOCK_SIZE} changed during serialization")
            connection.execute("INSERT INTO blocks VALUES (?, ?)",
                               (start // BLOCK_SIZE, payload))
            search_plans.add_block(start // BLOCK_SIZE, block)
            connection.executemany("INSERT INTO entries VALUES (?, ?, ?, ?)",
                                   ((start + offset, record.precursor_mz, class_ids[record.compound_class], adduct_ids[record.adduct])
                                    for offset, record in enumerate(block)))
        search_plans.finish(metadata)
        connection.execute("INSERT INTO metadata VALUES ('library', ?)", (json.dumps(metadata),))
        connection.commit()
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Library index failed SQLite integrity verification")
        connection.close()
        connection = None
        temporary.replace(path)
        return metadata
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)


def write_streaming_library_index(path: Path, records: Iterable[LibraryRecord], identity: dict) -> dict:
    """Freshly parse/normalize records, sort on disk, and verify one block at a time.

    Equal precursor masses retain input order exactly as in write_library_index.
    This builder changes neither the runtime format nor the four-block cache.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix=path.stem + '_stage_') as directory:
        stage = sqlite3.connect(Path(directory) / 'records.sqlite')
        output = None
        temporary = Path(directory) / 'verified.sqlite'
        try:
            stage.execute('PRAGMA journal_mode=OFF')
            stage.execute('PRAGMA synchronous=OFF')
            stage.execute('PRAGMA temp_store=FILE')
            stage.execute(f'PRAGMA cache_size=-{SQL_CACHE_KIB}')
            stage.execute('CREATE TABLE records (position INTEGER PRIMARY KEY, precursor REAL NOT NULL, payload BLOB NOT NULL)')
            count = 0
            classes, adducts = set(), set()
            batch = []
            for record in records:
                classes.add(record.compound_class)
                adducts.add(record.adduct)
                batch.append((count, record.precursor_mz, pickle.dumps(_pack_record(record), protocol=5)))
                count += 1
                if len(batch) == BLOCK_SIZE:
                    stage.executemany('INSERT INTO records VALUES (?, ?, ?)', batch)
                    stage.commit()
                    batch.clear()
            if batch:
                stage.executemany('INSERT INTO records VALUES (?, ?, ?)', batch)
                stage.commit()
                batch.clear()
            stage.execute('CREATE INDEX mass_order ON records (precursor, position)')
            stage.commit()
            metadata = dict(identity, index_format=INDEX_FORMAT_VERSION, block_size=BLOCK_SIZE,
                            record_encoding=RECORD_ENCODING, record_count=count,
                            classes=sorted(classes), adducts=sorted(adducts))
            class_ids = {value: index for index, value in enumerate(metadata['classes'])}
            adduct_ids = {value: index for index, value in enumerate(metadata['adducts'])}
            output = sqlite3.connect(temporary)
            output.execute('PRAGMA journal_mode=OFF')
            output.execute('PRAGMA synchronous=OFF')
            output.execute(f'PRAGMA cache_size=-{SQL_CACHE_KIB}')
            output.executescript(
                'CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);'
                'CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, precursor REAL NOT NULL, '
                'class_id INTEGER NOT NULL, adduct_id INTEGER NOT NULL);'
                'CREATE TABLE blocks (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL);')
            from .search_plan import SearchPlanWriter
            search_plans = SearchPlanWriter(output)
            cursor = stage.execute('SELECT payload FROM records ORDER BY precursor, position')
            verified = 0
            while True:
                rows = cursor.fetchmany(BLOCK_SIZE)
                if not rows:
                    break
                block = [_unpack_record(safe_pickle.loads(row[0])) for row in rows]
                payload = zlib.compress(_encode_block(block), level=6)
                if _decode_block(payload, metadata) != block:
                    raise ValueError(f'Library block {verified // BLOCK_SIZE} changed during serialization')
                output.execute('INSERT INTO blocks VALUES (?, ?)', (verified // BLOCK_SIZE, payload))
                search_plans.add_block(verified // BLOCK_SIZE, block)
                output.executemany('INSERT INTO entries VALUES (?, ?, ?, ?)',
                    ((verified + offset, r.precursor_mz, class_ids[r.compound_class], adduct_ids[r.adduct])
                     for offset, r in enumerate(block)))
                verified += len(block)
            if verified != count:
                raise ValueError('Streaming library index lost records')
            search_plans.finish(metadata)
            output.execute("INSERT INTO metadata VALUES ('library', ?)", (json.dumps(metadata),))
            output.commit()
            if output.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
                raise ValueError('Streaming library index failed integrity verification')
            output.close()
            output = None
            stage.close()
            stage = None
            temporary.replace(path)
            return metadata
        finally:
            if output is not None:
                output.close()
            if stage is not None:
                stage.close()


def repack_library_index(source: Path, destination: Path, identity: dict, *, expected_identity: dict) -> dict:
    """Stream a storage-only migration of already normalized records.

    Callers must establish that normalization is unchanged. This is not a
    replacement for reparsing a library after library-content or rule changes.
    Pin the input identity explicitly, retain its provenance, verify every
    record and precursor entry, and keep only one block in working memory.
    """
    from .provenance import sha256

    source, destination = Path(source), Path(destination)
    if not {"version", "source_sha256", "rules_sha256"}.issubset(expected_identity):
        raise ValueError("Source library identity must include version, source and code fingerprints")
    if source.resolve() == destination.resolve():
        raise ValueError("Repacking requires a separate source index")
    original = _readonly_connection(source)
    connection = None
    temporary = None
    try:
        prior = _read_metadata(original)
        if not all(prior.get(key) == value for key, value in expected_identity.items()):
            raise ValueError("Source library index identity changed")
        if any(prior.get(key) != identity.get(key) for key in ("version", "source_sha256")):
            raise ValueError("Repacking cannot change the normalized library version or source")
        metadata = dict(prior, **identity, index_format=INDEX_FORMAT_VERSION,
                        record_encoding=RECORD_ENCODING,
                        normalized_from_index_sha256=sha256(source),
                        normalized_from_rules_sha256=prior["rules_sha256"])
        if original.execute("SELECT COUNT(*) FROM entries").fetchone()[0] != metadata["record_count"]:
            raise ValueError("Source index record count is inconsistent")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=destination.stem, suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
        connection = sqlite3.connect(temporary)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute(f"PRAGMA cache_size=-{SQL_CACHE_KIB}")
        connection.executescript(
            "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);"
            "CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, precursor REAL NOT NULL, "
            "class_id INTEGER NOT NULL, adduct_id INTEGER NOT NULL);"
            "CREATE TABLE blocks (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL);"
        )
        from .search_plan import SearchPlanWriter
        search_plans = SearchPlanWriter(connection)
        class_ids = {value: index for index, value in enumerate(metadata["classes"])}
        adduct_ids = {value: index for index, value in enumerate(metadata["adducts"])}
        verified, last_precursor = 0, -math.inf
        for expected_block, (block_id, payload) in enumerate(original.execute("SELECT block_id, payload FROM blocks ORDER BY block_id")):
            if block_id != expected_block or verified != block_id * metadata["block_size"]:
                raise ValueError("Source index block sequence is inconsistent")
            records = _decode_block(payload, prior)
            if not records or len(records) > metadata["block_size"]:
                raise ValueError("Source index block length is inconsistent")
            encoded = zlib.compress(_encode_block(records), level=6)
            if _decode_block(encoded, metadata) != records:
                raise ValueError(f"Library block {block_id} changed during repacking")
            entries = [(verified + offset, r.precursor_mz, class_ids[r.compound_class], adduct_ids[r.adduct])
                       for offset, r in enumerate(records)]
            stored = original.execute(
                "SELECT ordinal, precursor, class_id, adduct_id FROM entries WHERE ordinal>=? AND ordinal<? ORDER BY ordinal",
                (verified, verified + len(records)),
            ).fetchall()
            if entries != stored or any(r.precursor_mz < last_precursor for r in records):
                raise ValueError("Source records disagree with the precursor index")
            if any(a.precursor_mz > b.precursor_mz for a, b in zip(records, records[1:])):
                raise ValueError("Source record order is inconsistent")
            last_precursor = records[-1].precursor_mz
            connection.execute("INSERT INTO blocks VALUES (?, ?)", (block_id, encoded))
            search_plans.add_block(block_id, records)
            connection.executemany("INSERT INTO entries VALUES (?, ?, ?, ?)", entries)
            verified += len(records)
        if verified != metadata["record_count"]:
            raise ValueError("Source index has missing record blocks")
        search_plans.finish(metadata)
        connection.execute("INSERT INTO metadata VALUES ('library', ?)", (json.dumps(metadata),))
        connection.commit()
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Repacked library index failed integrity verification")
        connection.close()
        connection = None
        temporary.replace(destination)
        return metadata
    finally:
        original.close()
        if connection is not None:
            connection.close()
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class IndexedLibrary(Sequence):
    """Four decoded blocks plus a bounded cache of individually used records."""

    def __init__(self, path: Path, *, allowed_adducts=(), allowed_class_keys=(), mz_min=None, mz_max=None):
        self.path = Path(path)
        self._connection = _readonly_connection(self.path)
        self._blocks = OrderedDict()
        self._selected_records = OrderedDict()
        self._selected_positions = []
        self._selected_bytes = 0
        self.precursors = array("d")
        self._positions = array("I")
        try:
            self.metadata = _read_metadata(self._connection)
            classes, adducts = self.metadata["classes"], self.metadata["adducts"]
            self.available_classes = tuple(value for value in classes if value)
            self.available_adducts = tuple(value for value in adducts if value)
            conditions = ["1=1"]
            parameters = []
            excluded = [i for i, value in enumerate(classes) if value in {"PS", "LPS"}]
            if excluded and "[M+NH4]+" in adducts:
                conditions.append("NOT (class_id IN (" + ",".join("?" for _ in excluded) + ") AND adduct_id=?)")
                parameters.extend([*excluded, adducts.index("[M+NH4]+")])
            if allowed_adducts:
                values = [i for i, value in enumerate(adducts) if value in allowed_adducts]
                conditions.append("adduct_id IN (" + ",".join("?" for _ in values) + ")")
                parameters.extend(values)
            if allowed_class_keys:
                values = [i for i, value in enumerate(classes)
                          if re.sub(r"[^A-Za-z0-9]+", "", value.upper()) in allowed_class_keys]
                conditions.append("class_id IN (" + ",".join("?" for _ in values) + ")")
                parameters.extend(values)
            for value, comparison in ((mz_min, ">="), (mz_max, "<=")):
                if value is not None:
                    conditions.append("precursor " + comparison + " ?")
                    parameters.append(value)
            query = "SELECT ordinal, precursor FROM entries WHERE " + " AND ".join(conditions) + " ORDER BY ordinal"
            for ordinal, precursor in self._connection.execute(query, parameters):
                self._positions.append(ordinal)
                self.precursors.append(precursor)
        except Exception:
            self.close()
            raise

    def __len__(self):
        return len(self._positions)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[position] for position in range(*index.indices(len(self)))]
        position = self._positions[operator.index(index)]
        selected = self._selected_records.get(position)
        if selected is not None:
            self._selected_records.move_to_end(position)
            return selected[0]
        block_size = self.metadata["block_size"]
        block_id, offset = divmod(position, block_size)
        block = self._get_block(block_id)
        record = block.record(offset)
        # Keep only selected annotations, never an experimental spectrum or
        # query result. This avoids decoding an entire block on repeat visits.
        revision = record_revision(record)
        size = (_record_cache_bytes(record) + sys.getsizeof(revision)
                + sum(sys.getsizeof(f) for f in revision[9]) + sys.getsizeof(position) + 8)
        if size <= MAX_SELECTED_RECORD_BYTES:
            self._selected_records[position] = (record, size, revision)
            bisect.insort(self._selected_positions, position)
            self._selected_bytes += size
            while (len(self._selected_records) > MAX_SELECTED_RECORDS
                   or self._selected_bytes > MAX_SELECTED_RECORD_BYTES):
                removed_position, (_, removed, _) = self._selected_records.popitem(last=False)
                self._selected_positions.pop(bisect.bisect_left(self._selected_positions, removed_position))
                self._selected_bytes -= removed
        return record

    def _get_block(self, block_id):
        block = self._blocks.get(block_id)
        if block is None:
            if self._connection is None:
                raise RuntimeError("Library index is closed")
            row = self._connection.execute("SELECT payload FROM blocks WHERE block_id=?", (block_id,)).fetchone()
            if row is None:
                raise ValueError(f"Library index is missing block {block_id}")
            block = _RecordBlock(row[0], self.metadata)
            self._blocks[block_id] = block
            if len(self._blocks) > MAX_CACHED_BLOCKS:
                self._blocks.popitem(last=False)
        else:
            self._blocks.move_to_end(block_id)
        return block

    def candidate_record_blocks(self, left, right):
        """Yield packed blocks and visible/local indexes without decoding records."""
        block_size = self.metadata["block_size"]
        while left < right:
            block_id = self._positions[left] // block_size
            indexes, offsets = [], []
            while left < right and self._positions[left] // block_size == block_id:
                indexes.append(left)
                offsets.append(self._positions[left] % block_size)
                left += 1
            yield self._get_block(block_id), indexes, offsets

    def close(self):
        if self._connection is not None:
            self._connection.close()
            self._connection = None
        self._blocks.clear()
        self._selected_records.clear()
        self._selected_positions.clear()
        self._selected_bytes = 0

    def __del__(self):
        if hasattr(self, "_connection"):
            self.close()


def prebuilt_index_path(library_path: Path) -> Path | None:
    """Validate source and code identity without decoding any library records."""
    from .library import LIBRARY_CACHE_VERSION, _is_bundled_library
    from .provenance import code_fingerprint, sha256

    library_path = Path(library_path)
    if not library_path.is_file():
        return None
    if library_path.suffix.lower() == ".sqlite":
        bank = library_path
        expected = {"version": LIBRARY_CACHE_VERSION,
                    "rules_sha256": code_fingerprint(Path(__file__).parent)}
    elif library_path.name not in {"current_positive.msp.gz", "current_negative.msp.gz", "sphingolipid_positive.msp.gz"}:
        return None
    elif _is_bundled_library(library_path):
        base = Path(sys._MEIPASS) / "libraries" / "ms2" / "prebuilt"
    else:
        root = Path(__file__).resolve().parents[3]
        if library_path.resolve().parent != (root / "libraries" / "ms2").resolve():
            return None
        base = root / "build" / "prebuilt_libraries"
    if library_path.suffix.lower() != ".sqlite":
        bank = base / (library_path.name.removesuffix(".msp.gz") + ".sqlite")
        expected = {"version": LIBRARY_CACHE_VERSION, "source_sha256": sha256(library_path),
                    "rules_sha256": code_fingerprint(Path(__file__).parent)}
    if not bank.is_file():
        return None
    try:
        connection = _readonly_connection(bank)
        try:
            metadata = _read_metadata(connection)
            if all(metadata.get(key) == value for key, value in expected.items()):
                return bank
        finally:
            connection.close()
    except (OSError, ValueError, sqlite3.Error):
        pass
    if _is_bundled_library(library_path) or library_path.suffix.lower() == ".sqlite":
        raise ValueError("内置谱库索引校验失败，请使用完整的 LipidGate 发布包")
    return None


def open_indexed_library(library_path, *, allowed_adducts=(), allowed_class_keys=(), mz_min=None, mz_max=None):
    bank = prebuilt_index_path(Path(library_path))
    if bank is None:
        return None
    return IndexedLibrary(bank, allowed_adducts=allowed_adducts, allowed_class_keys=allowed_class_keys,
                          mz_min=mz_min, mz_max=mz_max)
