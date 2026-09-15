"""Parse the final manuscript's residue/parenthesis notation, not invented #.

Source: 中文初稿改1(1).docx, Figure 4 and adjacent paragraph:
Gal-Gal(-Fuc)-GlcNAc-Gal-Glc. Parentheses attach a branch to the immediately
preceding main-chain residue; order is non-reducing end towards ceramide.
Nested branching is not defined in that source and is explicitly rejected.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re

RESIDUES = frozenset({"Fuc","Gal","Glc","Hex","GalNAc","GlcNAc","KDN","NeuGc","NeuAc"})


@dataclass(frozen=True)
class GlycanBranch:
    position: int  # zero-based index of attachment residue in main_chain
    residues: tuple[str,...]


@dataclass(frozen=True)
class GlycanEncoding:
    main_chain: tuple[str,...]
    branches: tuple[GlycanBranch,...] = ()
    source_format: str = "residue_sequence"


def _residue(value: str) -> str:
    if value not in RESIDUES:
        raise ValueError(f"Unknown glycan residue: {value!r}")
    return value


def parse_glycan_encoding(text: str, legacy_locations: str = "") -> GlycanEncoding:
    """Normalize final-paper notation and historical integer-position input.

    Historical indices refer to literal split(' ') tokens, including a leading
    empty token. A nonempty branch token attaches to the following main residue;
    index 0 on the empty token is the historical sentinel, not a residue.
    """
    if not isinstance(text,str) or not text.strip():
        raise ValueError("A glycan sequence must not be empty")
    if "#" in text:
        raise ValueError("# grammar is unverified; use final-manuscript parentheses")
    if legacy_locations or (text.lstrip().startswith("-") and "(" not in text):
        if "(" in text or ")" in text:
            raise ValueError("Do not mix integer positions and parenthesis notation")
        tokens = text.split(" ")
        for token in tokens:
            if token:
                if not token.startswith("-"):
                    raise ValueError("Legacy residues must start with '-'")
                _residue(token[1:])
        if legacy_locations and not re.fullmatch(r"\d+(?:,\d+)*",str(legacy_locations)):
            raise ValueError("Legacy positions must be comma-separated non-negative integers")
        locations = [int(x) for x in str(legacy_locations).split(",")] if legacy_locations else []
        if len(locations) != len(set(locations)) or any(i >= len(tokens) for i in locations):
            raise ValueError("Duplicate or out-of-range legacy branch position")
        branch_indices = {i for i in locations if tokens[i]}
        main_indices = [i for i,t in enumerate(tokens) if t and i not in branch_indices]
        branches = []
        for i in sorted(branch_indices):
            following = next((j for j in main_indices if j > i),None)
            if following is None:
                raise ValueError("Legacy branch has no following attachment residue")
            branches.append(GlycanBranch(main_indices.index(following),(tokens[i][1:],)))
        return GlycanEncoding(tuple(tokens[i][1:] for i in main_indices),tuple(branches),"legacy_positions")

    sequence = re.sub(r"\s+","",text)
    main,branches = [],[]
    i = 0
    while i < len(sequence):
        match = re.match(r"[A-Za-z]+",sequence[i:])
        if not match:
            raise ValueError(f"Expected glycan residue at position {i}")
        main.append(_residue(match.group()))
        i += len(match.group())
        while i < len(sequence) and sequence[i] == "(":
            end = sequence.find(")",i+1)
            if end == -1 or "(" in sequence[i+1:end]:
                raise ValueError("Unclosed or nested glycan branch")
            branch = sequence[i+1:end]
            if not branch.startswith("-") or not branch[1:]:
                raise ValueError("Branch must have the manuscript form (-Residue)")
            if "-" in branch[1:] or any(b.position == len(main)-1 for b in branches):
                raise ValueError("Only one single-residue branch per attachment is documented")
            residues = (_residue(branch[1:]),)
            branches.append(GlycanBranch(len(main)-1,residues))
            i = end+1
        if i < len(sequence):
            if sequence[i] != "-" or i == len(sequence)-1:
                raise ValueError(f"Expected '-' and a residue at position {i}")
            i += 1
    if not main:
        raise ValueError("Missing glycan main chain")
    return GlycanEncoding(tuple(main),tuple(branches))


def generate_glycan_losses(encoding: GlycanEncoding, precursor: float, name: str, masses: dict) -> list[tuple[str,float]]:
    """Enumerate main-chain losses and branch-first/retaining loss paths.

    Masses are supplied by GSL_fragments, preserving its residue/water values.
    A branch leaves with its attachment residue unless already removed. Partial
    branch-chain losses run from the non-reducing end, without double subtraction.
    The caller retains legacy NeuAc diagnostic masses and deduplicates m/z.
    """
    output = []
    branches = encoding.branches
    def emit(losses, path):
        mass = sum(masses['-'+r] for r in losses)
        label = name + "-" + "-".join(losses) + " ["+path+"]"
        output.extend([(label,precursor-mass),(label+"-H2O",precursor-mass-18.0106)])
    # Complete branch removal combinations are sufficient to define residual
    # main chains; branch-local partial losses are also explicit.
    for branch_index,branch in enumerate(branches):
        for count in range(1,len(branch.residues)):
            emit(branch.residues[:count],f"branch {branch_index+1} partial")
    for size in range(len(branches)+1):
        for removed_tuple in combinations(range(len(branches)),size):
            removed = set(removed_tuple)
            branch_losses = [r for j in removed_tuple for r in branches[j].residues]
            if removed:
                emit(branch_losses,"branch-first "+','.join(str(j+1) for j in removed_tuple))
            for count in range(1,len(encoding.main_chain)+1):
                lost_branches = removed | {j for j,b in enumerate(branches) if b.position < count}
                losses = list(encoding.main_chain[:count]) + [r for j in sorted(lost_branches) for r in branches[j].residues]
                retained = sorted(set(range(len(branches)))-lost_branches)
                path = f"main {count}; removed branches {sorted(removed)}; retained branches {retained}"
                emit(losses,path)
    return output
