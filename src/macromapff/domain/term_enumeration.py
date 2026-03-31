#!/usr/bin/env python3
from itertools import combinations


def _norm_rev_tuple(values):
    """Return lexicographically smaller ordering between tuple and reverse."""
    vals = tuple(values)
    rev = tuple(reversed(vals))
    return vals if vals <= rev else rev


def enumerate_terms(mol):
    """Enumerate bonds, angles, dihedrals, and impropers from one molecule."""
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx() + 1
        j = bond.GetEndAtomIdx() + 1
        bonds.append((i, j))

    angles = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        for i, k in combinations(sorted(nbs), 2):
            angles.append((i, center + 1, k))

    dihedral_set = set()
    for bond in mol.GetBonds():
        j = bond.GetBeginAtomIdx() + 1
        k = bond.GetEndAtomIdx() + 1
        j_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(j - 1).GetNeighbors()
            if x.GetIdx() + 1 != k
        ]
        k_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(k - 1).GetNeighbors()
            if x.GetIdx() + 1 != j
        ]
        for i in j_neighbors:
            for l in k_neighbors:
                t = (i, j, k, l)
                dihedral_set.add(_norm_rev_tuple(t))
    dihedrals = sorted(dihedral_set)

    impropers = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        if len(nbs) < 3:
            continue
        for a, b, c in combinations(sorted(nbs), 3):
            impropers.append((center + 1, a, b, c))

    return bonds, angles, dihedrals, impropers
