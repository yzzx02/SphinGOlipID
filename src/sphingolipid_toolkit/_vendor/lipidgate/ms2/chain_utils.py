from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache
from typing import Callable, Sequence

from .models import FragmentMatch, FragmentRecord, LibraryRecord


CHAIN_TOKEN_RE = re.compile(
    r"(?:[mdt]|O-|P-)?\d+:\d+(?:\((?:\d+)?OH\)|;\d*O(?:H)?|,O\d*)?",
    flags=re.IGNORECASE,
)
FRAGMENT_CHAIN_TOKEN_RE = re.compile(
    r"(?P<prefix>[OP]-)?(?P<base>\d+:\d+)"
    r"(?:(?:\((?P<paren_ox>\d*)O\))|(?:,O(?P<comma_ox>\d*))|"
    r"(?:;O(?P<semicolon_ox>\d*))|"
    r"(?:;\(?(?P<oh_count>\d+)OH\)?))?"
)


def canonical_chain_token(match: re.Match[str]) -> str:
    """Normalize oxygen suffixes used by fragment annotations."""

    prefix = match.group("prefix") or ""
    token = f"{prefix}{match.group('base')}"
    oxygen_count = match.group("paren_ox")
    if oxygen_count is None:
        oxygen_count = match.group("comma_ox")
    if oxygen_count is None:
        oxygen_count = match.group("semicolon_ox")
    if oxygen_count is not None:
        return f"{token};O{oxygen_count or '1'}"
    hydroxy_count = match.group("oh_count")
    if hydroxy_count is not None:
        return f"{token};{hydroxy_count}OH"
    return token


def extract_chain_tokens(name: object) -> list[str]:
    """Extract every explicit non-zero chain from a lipid name."""
    # Cache immutable text results; callers still receive independent lists.
    return list(_chain_tokens_text(str(name or "")))


@lru_cache(maxsize=4096)
def _chain_tokens_text(text: str) -> tuple[str, ...]:
    return tuple(
        match.group(0)
        for match in CHAIN_TOKEN_RE.finditer(text)
        if match.group(0) != "0:0"
    )


def extract_fragment_chain_token(fragment_or_name: FragmentRecord | object) -> str | None:
    """Return the first chain token carried by a fragment annotation."""

    name = (
        fragment_or_name.name
        if isinstance(fragment_or_name, FragmentRecord)
        else fragment_or_name
    )
    return _fragment_chain_token_text(str(name or ""))


@lru_cache(maxsize=4096)
def _fragment_chain_token_text(text: str) -> str | None:
    matched = FRAGMENT_CHAIN_TOKEN_RE.search(text)
    return canonical_chain_token(matched) if matched is not None else None


def chain_token_multiplicity(record_or_name: LibraryRecord | object) -> Counter[str]:
    """Count structural chain positions while retaining repeated chains."""

    name = (
        record_or_name.lipid_chain_name
        if isinstance(record_or_name, LibraryRecord)
        else record_or_name
    )
    return Counter(dict(_chain_count_items_text(str(name or ""))))


@lru_cache(maxsize=4096)
def _chain_count_items_text(text: str) -> tuple[tuple[str, int], ...]:
    counts: dict[str, int] = {}
    for token in _chain_tokens_text(text):
        counts[token] = counts.get(token, 0) + 1
    return tuple(counts.items())


def chain_token_counts(record: LibraryRecord) -> dict[str, int]:
    """Return independent counts for read-only scoring, without Counter setup."""
    return dict(_chain_count_items_text(str(record.lipid_chain_name or "")))


def count_chain_evidence(
    record: LibraryRecord,
    fragments_or_matches: Sequence[FragmentRecord | FragmentMatch],
    predicate: Callable[[LibraryRecord, FragmentRecord], bool],
) -> int:
    """Count chain positions using a caller's evidence definition.

    Alternative ions of one chain count once; repeated structural positions
    retain their multiplicity. Unlabelled qualifying ions retain the existing
    count, capped by the number of structural positions when known.
    """
    multiplicity = None
    matched_tokens: set[str] = set()
    unassigned = 0
    for item in fragments_or_matches:
        fragment = item.fragment if isinstance(item, FragmentMatch) else item
        if not predicate(record, fragment):
            continue
        if multiplicity is None:
            multiplicity = chain_token_counts(record)
        token = extract_fragment_chain_token(fragment)
        if token is not None and token in multiplicity:
            matched_tokens.add(token)
        else:
            unassigned += 1
    if not multiplicity:
        return unassigned
    count = sum(multiplicity[token] for token in matched_tokens) + unassigned
    return min(count, sum(multiplicity.values()))
