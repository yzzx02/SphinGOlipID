"""Shared deterministic centroid-to-theory assignment for all MS/MS runners.

Each candidate is a separate structural hypothesis. Within a candidate, solve
maximum-cardinality matching, then minimum total absolute ppm error, then
maximum summed observed intensity. Exact rational costs avoid epsilon weights
that could silently trade ppm accuracy for intensity. No scientific threshold
other than the caller's inclusive ppm window is introduced.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Sequence

import numpy as np
import pandas as pd


def _assignment(observed, theoretical, pairs):
    """Successive shortest augmenting paths with lexicographic exact costs."""
    source = 0
    obs_offset = 1
    theo_offset = obs_offset + len(observed)
    sink = theo_offset + len(theoretical)
    graph = [[] for _ in range(sink+1)]
    zero = (Fraction(0), Fraction(0))

    def edge(a,b,cost):
        forward = [b,len(graph[b]),1,cost]
        backward = [a,len(graph[a]),0,tuple(-c for c in cost)]
        graph[a].append(forward)
        graph[b].append(backward)
        return forward

    for i in range(len(observed)):
        edge(source,obs_offset+i,zero)
    for j in range(len(theoretical)):
        edge(theo_offset+j,sink,zero)
    tracked = []
    for i,j,row_index in pairs:
        mz,intensity = observed[i]
        tmz = theoretical[j]
        ppm = abs(Fraction(mz)-Fraction(tmz))/abs(Fraction(tmz))*1_000_000
        forward = edge(obs_offset+i,theo_offset+j,(ppm,-Fraction(intensity)))
        tracked.append((forward,row_index))

    while True:
        distance = [None]*len(graph)
        parent = [None]*len(graph)
        distance[source] = zero
        # Reverse edges can have negative costs. Bellman-Ford also allows an
        # earlier assignment to move, unlike greedy nearest-peak selection.
        for _ in range(len(graph)-1):
            changed = False
            for a,neighbors in enumerate(graph):
                if distance[a] is None:
                    continue
                for k,(b,rev,capacity,cost) in enumerate(neighbors):
                    if not capacity:
                        continue
                    proposal = tuple(x+y for x,y in zip(distance[a],cost))
                    if distance[b] is None or proposal < distance[b]:
                        distance[b] = proposal
                        parent[b] = (a,k)
                        changed = True
            if not changed:
                break
        if parent[sink] is None:
            break
        b = sink
        while b != source:
            a,k = parent[b]
            forward = graph[a][k]
            forward[2] = 0
            graph[b][forward[1]][2] = 1
            b = a
    return [row_index for forward,row_index in tracked if forward[2] == 0]


def select_one_to_one_edges(
    edges: pd.DataFrame,
    *,
    observed_cols: Sequence[str],
    observed_mz_col: str,
    intensity_col: str,
    theoretical_mz_col: str,
    candidate_cols: Sequence[str] = ("anno",),
    label_col: str | None = None,
) -> pd.DataFrame:
    """Select admissible edges; callers construct inclusive ppm bounds first.

    ``observed_cols`` identify a centroid, including scan/position where known.
    Identical theoretical masses form one node, even with multiple labels.
    Labels are joined in lexical order for explanation only, never extra score.
    """
    if edges.empty:
        return edges.copy()
    required = set(observed_cols) | set(candidate_cols) | {observed_mz_col,intensity_col,theoretical_mz_col}
    if not required.issubset(edges):
        raise KeyError(f"Missing assignment columns: {sorted(required-set(edges))}")
    data = edges.reset_index(drop=True).copy()
    numeric = data[[observed_mz_col,intensity_col,theoretical_mz_col]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or (numeric[:,2] <= 0).any():
        raise ValueError("Assignment requires finite masses/intensities and positive theoretical m/z")
    groups = data.groupby(list(candidate_cols),sort=True,dropna=False) if candidate_cols else [(None,data)]
    outputs = []
    for _,group in groups:
        group = group.copy()
        if label_col and label_col in group:
            labels = group.groupby(theoretical_mz_col)[label_col].agg(
                lambda values: " | ".join(sorted({str(v) for v in values if pd.notna(v)})))
            group[label_col] = group[theoretical_mz_col].map(labels)
        obs_map = {}
        for idx,row in group.iterrows():
            key = tuple(row[c] for c in observed_cols)
            values = (float(row[observed_mz_col]),float(row[intensity_col]))
            if key in obs_map and obs_map[key] != values:
                raise ValueError("An observed centroid identity has inconsistent mass/intensity")
            obs_map[key] = values
        obs_keys = sorted(obs_map,key=lambda key:(obs_map[key][0],-obs_map[key][1],repr(key)))
        theo_keys = sorted(set(group[theoretical_mz_col].astype(float)))
        oi = {k:i for i,k in enumerate(obs_keys)}
        ti = {k:i for i,k in enumerate(theo_keys)}
        pair_rows = {}
        # Canonical row content resolves equivalent metadata rows, never order.
        for idx in sorted(group.index,key=lambda i:repr(tuple(group.loc[i]))):
            row = group.loc[idx]
            pair = (oi[tuple(row[c] for c in observed_cols)],ti[float(row[theoretical_mz_col])])
            pair_rows.setdefault(pair,idx)
        pairs = [(i,j,pair_rows[i,j]) for i,j in sorted(pair_rows)]
        selected = _assignment([obs_map[k] for k in obs_keys],theo_keys,pairs)
        outputs.append(group.loc[selected].sort_values([theoretical_mz_col,observed_mz_col],kind="stable"))
    return pd.concat(outputs,ignore_index=True) if outputs else data.iloc[:0]
