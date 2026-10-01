"""Outcome-blind label hierarchies for the adult variables, for tree-structured searches.

The divisive search (``SplitSearch``) needs ``Ordinal`` or a **binary** ``Tree`` whose
groups are contiguous in label order, so an unordered variable has to be given a
hierarchy and recoded so its labels follow the tree's leaf order. Recoding is a
relabelling of the population, applied before anything is simulated, so every
exactness property of :mod:`adult_semisynth` survives unchanged.

Occupation
----------
Adult's 14 occupation labels are the detailed groups of the **1990 Census
occupational classification** (the one the 1994 CPS used), which arranges them in
six major groups. Above those, the standard white-collar / manual-and-service split.
Two three-member major groups need one binary choice each; those are the only
choices not taken from the classification, and are marked below. Armed-Forces
(9 rows) sits with Protective-serv.

    root
    +- White collar
    |  +- Managerial & professional:   Exec-managerial | Prof-specialty
    |  +- Technical, sales & admin:    Sales | (Tech-support | Adm-clerical)   [*]
    +- Manual & service
       +- Service:                     (Protective-serv | Armed-Forces) | (Other-service | Priv-house-serv)
       +- Manual:                      Farming-fishing | (Craft-repair | Operators)
          Operators, fabricators, laborers: Machine-op-inspct | (Transport-moving | Handlers-cleaners)   [*]

[*] binarisation choices. Support (technical + clerical) vs sales; machine operators
vs the non-machine manual groups.

A level pooled by :func:`adult_semisynth.pool_rare` (``"Other (a, b, ...)"``) takes the
place of its most frequent member.
"""

from __future__ import annotations

import numpy as np

import adult
from catci.structure import Tree

OCCUPATION_PARENTS = {
    # leaves
    "Exec-managerial": "MgrProf", "Prof-specialty": "MgrProf",
    "Sales": "TSA", "Tech-support": "Support", "Adm-clerical": "Support",
    "Protective-serv": "Protect", "Armed-Forces": "Protect",
    "Other-service": "SvcOther", "Priv-house-serv": "SvcOther",
    "Farming-fishing": "Manual", "Craft-repair": "BlueCollar",
    "Machine-op-inspct": "Operators", "Transport-moving": "Laborers",
    "Handlers-cleaners": "Laborers",
    # internal nodes
    "Support": "TSA", "MgrProf": "WhiteCollar", "TSA": "WhiteCollar",
    "Protect": "Service", "SvcOther": "Service",
    "Laborers": "Operators", "Operators": "BlueCollar", "BlueCollar": "Manual",
    "Service": "ManualService", "Manual": "ManualService",
    "WhiteCollar": None, "ManualService": None,
}

TAXONOMIES = {"Occupation": OCCUPATION_PARENTS}


def _members(level: str) -> list[str]:
    if level.startswith("Other (") and level.endswith(")"):
        return [m.strip() for m in level[len("Other ("):-1].split(",")]
    return [level]


def _dfs_leaves(parents: dict, present: set[str]) -> list[str]:
    children: dict = {}
    for code, parent in parents.items():
        children.setdefault(parent, []).append(code)
    out, stack = [], [None]
    while stack:
        node = stack.pop()
        if node in present:
            out.append(node)
        stack.extend(reversed(children.get(node, [])))
    return out


def apply_taxonomy(data: adult.AdultData, name: str) -> tuple[adult.AdultData, Tree]:
    """Recode ``name`` so its labels follow the taxonomy's leaf order; return the tree.

    The returned tree is binary and its groups are contiguous label ranges, as the
    divisive search requires. A pooled level stands in for its most frequent member
    (frequencies from the unpooled data).
    """
    parents = TAXONOMIES[name]
    leaves = {k for k in parents if k not in set(parents.values())}
    var = data.variables[name]

    raw = adult.load()
    raw_freq = raw.codes[name].value_counts()
    raw_levels = raw.variables[name].levels

    def stand_in(level: str) -> str:
        members = _members(level)
        best = max(members, key=lambda m: int(raw_freq.get(raw_levels.index(m) + 1, 0)))
        if best not in leaves:
            raise ValueError(f"{name} level {level!r} has no place in the taxonomy")
        return best

    leaf_of = {lev: stand_in(lev) for lev in var.levels}
    order = _dfs_leaves(parents, set(leaf_of.values()))
    level_of = {leaf: lev for lev, leaf in leaf_of.items()}
    new_levels = tuple(level_of[leaf] for leaf in order)
    old_to_new = {var.levels.index(lev) + 1: k for k, lev in enumerate(new_levels, start=1)}

    codes = data.codes.copy()
    codes[name] = codes[name].map(old_to_new).astype(np.int64)
    variables = dict(data.variables)
    variables[name] = adult.Variable(var.name, var.code, var.kind, new_levels)

    # Labels 1..d replace the leaf names; internal codes are strings, so no clash.
    tree_parents = {k: parents[leaf_of[lev]] for k, lev in enumerate(new_levels, start=1)}
    tree_parents.update({c: p for c, p in parents.items() if c not in leaves})
    tree = Tree.from_parents(tree_parents, list(range(1, len(new_levels) + 1)))
    return adult.AdultData(codes=codes, variables=variables), tree
