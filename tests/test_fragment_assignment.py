"""Exhaustive small-graph oracle; no experimental spectra."""
from fractions import Fraction
from itertools import combinations

import pandas as pd
from sphingolipid_toolkit.fragment_assignment import _assignment, select_one_to_one_edges


def test_all_three_by_three_graphs_against_exhaustive_oracle():
    observed = [(100., 40.), (100.0002, 70.), (100.0004, 90.)]
    theoretical = [100.0001, 100.0003, 100.0005]
    full = [(i,j,i*3+j) for i in range(3) for j in range(3)]
    def cost(edges):
        return (-len(edges), sum((abs(Fraction(observed[i][0])-Fraction(theoretical[j])) /
            Fraction(theoretical[j])*1_000_000 for i,j,_ in edges), Fraction()),
            -sum(observed[i][1] for i,j,_ in edges))
    for mask in range(512):
        edges = [e for k,e in enumerate(full) if mask & (1 << k)]
        feasible = [subset for n in range(4) for subset in combinations(edges,n)
            if len({e[0] for e in subset}) == n == len({e[1] for e in subset})]
        selected = set(_assignment(observed,theoretical,edges))
        assert cost([e for e in edges if e[2] in selected]) == min(map(cost,feasible))


def test_intensity_tie_and_duplicate_labels_are_deterministic():
    edges = pd.DataFrame({"obs":[0,0,1,1],"omz":[100.]*4,"intensity":[40.,40.,80.,80.],
        "tmz":[100.]*4,"label":["b","a","b","a"],"anno":["A"]*4})
    for seed in range(5):
        selected = select_one_to_one_edges(edges.sample(frac=1,random_state=seed),
            observed_cols=["obs"],observed_mz_col="omz",intensity_col="intensity",
            theoretical_mz_col="tmz",label_col="label")
        assert selected.obs.tolist() == [1]
        assert selected.label.tolist() == ["a | b"]
