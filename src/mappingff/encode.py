"""Environment key encoding for mappingff.

This module provides functions for computing multi-level chemical environment
descriptors (hop0, hop1, hop2) and encoding them as SHA-256 keys.

The three levels of environment description:
    - hop0: Center atom properties + neighbor signatures (coarse-grained)
    - hop1: hop0 + detailed hop1 shell (first fallback)
    - hop2: hop1 + detailed hop2 shell (exact match)

Each environment is encoded as a SHA-256 hex string for use as a dictionary key
in the parameter database.
"""

from __future__ import annotations

import hashlib
import json
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit.Chem import rdchem


# ── Bond Type Helper ────────────────────────────────────────────────────────────


def _bondTypeCode(bond: rdchem.Bond) -> str:
    """Get single-character bond type code.

    Codes: S=SINGLE, D=DOUBLE, T=TRIPLE, A=AROMATIC, U=UNKNOWN
    """
    from rdkit import Chem
    bt = bond.GetBondType()
    if bt == Chem.rdchem.BondType.SINGLE:
        return "S"
    elif bt == Chem.rdchem.BondType.DOUBLE:
        return "D"
    elif bt == Chem.rdchem.BondType.TRIPLE:
        return "T"
    elif bt == Chem.rdchem.BondType.AROMATIC:
        return "A"
    else:
        return "U"


# ── Graph-Based Subgraph Extraction ────────────────────────────────────────────


def _atomProps(atom: rdchem.Atom) -> dict:
    """Get canonical properties dict for an atom."""
    return {
        "z": atom.GetAtomicNum(),
        "fc": atom.GetFormalCharge(),
        "ar": int(atom.GetIsAromatic()),
        "deg": atom.GetTotalDegree(),
        "h": atom.GetTotalNumHs(includeNeighbors=True),
        "ring": int(atom.IsInRing()),
    }


def getHop1Subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Extract hop1 subgraph centered on an atom.

    The subgraph includes:
    - center: center atom properties
    - hop1: {idx: {z, fc, ar, deg, h, ring, bt_to_center}}
    - bonds: [(center_idx, hop1_idx, bt), ...] including hop1-hop1 bonds

    Args:
        mol: RDKit molecule.
        atom: Center atom (RDKit 0-based index).

    Returns:
        dict with center, hop1 dict, and bonds list.
    """
    center_idx = atom.GetIdx()
    ring_info = mol.GetRingInfo()

    # Center properties
    center_props = _atomProps(atom)
    center_props["ring_count"] = ring_info.NumAtomRings(center_idx)

    # Hop1 atoms: neighbors of center
    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _atomProps(n)
        props["ring_count"] = ring_info.NumAtomRings(idx)
        # Bond type from hop1 to center
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bondTypeCode(bond)
        hop1[idx] = props

    # Bonds: center-hop1 and hop1-hop1
    bonds = []
    for hop1_idx in hop1_idx_set:
        bonds.append((center_idx, hop1_idx, hop1[hop1_idx]["bt_to_center"]))
        for n2 in mol.GetAtomWithIdx(hop1_idx).GetNeighbors():
            idx2 = n2.GetIdx()
            if idx2 in hop1_idx_set and idx2 > hop1_idx:
                bond = mol.GetBondBetweenAtoms(hop1_idx, idx2)
                bonds.append((hop1_idx, idx2, _bondTypeCode(bond)))

    return {
        "center": center_props,
        "hop1": hop1,
        "bonds": bonds,
    }


def getHop2Subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Extract hop2 subgraph centered on an atom.

    The subgraph includes:
    - center: center atom properties
    - hop1: {idx: {z, fc, ar, deg, h, ring, bt_to_center}}
    - hop2: {idx: {z, fc, ar, deg, h, ring, bt_to_parent, parent_idx}}
    - bonds: all bonds within hop2 scope (center-hop1, hop1-hop2, hop2-hop2, hop1-hop1)

    Args:
        mol: RDKit molecule.
        atom: Center atom (RDKit 0-based index).

    Returns:
        dict with center, hop1 dict, hop2 dict, and bonds list.
    """
    center_idx = atom.GetIdx()
    ring_info = mol.GetRingInfo()

    # Center properties
    center_props = _atomProps(atom)
    center_props["ring_count"] = ring_info.NumAtomRings(center_idx)

    # Hop1 atoms: neighbors of center
    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _atomProps(n)
        props["ring_count"] = ring_info.NumAtomRings(idx)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bondTypeCode(bond)
        hop1[idx] = props

    # Hop2 atoms: neighbors of hop1, excluding center and hop1
    hop2_idx_set = set()
    hop2 = {}
    for hop1_idx in hop1_idx_set:
        hop1_atom = mol.GetAtomWithIdx(hop1_idx)
        for n in hop1_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop1_idx_set or idx == center_idx:
                continue
            if idx not in hop2_idx_set:
                hop2_idx_set.add(idx)
                props = _atomProps(n)
                props["ring_count"] = ring_info.NumAtomRings(idx)
                # Bond to parent hop1
                bond = mol.GetBondBetweenAtoms(hop1_idx, idx)
                props["bt_to_parent"] = _bondTypeCode(bond)
                props["parent_idx"] = hop1_idx
                hop2[idx] = props

    bonds = []

    # center-hop1 bonds
    for hop1_idx in hop1_idx_set:
        bonds.append((center_idx, hop1_idx, hop1[hop1_idx]["bt_to_center"]))

    # hop1-hop1 bonds
    for hop1_idx in hop1_idx_set:
        for n2 in mol.GetAtomWithIdx(hop1_idx).GetNeighbors():
            idx2 = n2.GetIdx()
            if idx2 in hop1_idx_set and idx2 > hop1_idx:
                bond = mol.GetBondBetweenAtoms(hop1_idx, idx2)
                bonds.append((hop1_idx, idx2, _bondTypeCode(bond)))

    # hop1-hop2 bonds
    for hop2_idx in hop2_idx_set:
        hop2_atom = mol.GetAtomWithIdx(hop2_idx)
        for n in hop2_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop1_idx_set:
                bond = mol.GetBondBetweenAtoms(hop2_idx, idx)
                bonds.append((hop2_idx, idx, _bondTypeCode(bond)))

    # hop2-hop2 bonds
    for hop2_idx in hop2_idx_set:
        hop2_atom = mol.GetAtomWithIdx(hop2_idx)
        for n in hop2_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop2_idx_set and idx > hop2_idx:
                bond = mol.GetBondBetweenAtoms(hop2_idx, idx)
                bonds.append((hop2_idx, idx, _bondTypeCode(bond)))

    return {
        "center": center_props,
        "hop1": hop1,
        "hop2": hop2,
        "bonds": bonds,
    }


def getHop3Subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Extract hop3 subgraph centered on an atom.

    The subgraph includes:
    - center: center atom properties
    - hop1: {idx: {z, fc, ar, deg, h, ring, bt_to_center}}
    - hop2: {idx: {z, fc, ar, deg, h, ring, bt_to_parent, parent_idx}}
    - hop3: {idx: {z, fc, ar, deg, h, ring, bt_to_parent, parent_idx}}
    - bonds: all bonds within hop3 scope (center-hop1, hop1-hop2, hop2-hop3, hop1-hop1, hop2-hop2, hop3-hop3)

    Args:
        mol: RDKit molecule.
        atom: Center atom (RDKit 0-based index).

    Returns:
        dict with center, hop1, hop2, hop3 dicts, and bonds list.
    """
    center_idx = atom.GetIdx()
    ring_info = mol.GetRingInfo()

    # Center properties
    center_props = _atomProps(atom)
    center_props["ring_count"] = ring_info.NumAtomRings(center_idx)

    # Hop1 atoms: neighbors of center
    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _atomProps(n)
        props["ring_count"] = ring_info.NumAtomRings(idx)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bondTypeCode(bond)
        hop1[idx] = props

    # Hop2 atoms: neighbors of hop1, excluding center and hop1
    hop2_idx_set = set()
    hop2 = {}
    for hop1_idx in hop1_idx_set:
        hop1_atom = mol.GetAtomWithIdx(hop1_idx)
        for n in hop1_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop1_idx_set or idx == center_idx:
                continue
            if idx not in hop2_idx_set:
                hop2_idx_set.add(idx)
                props = _atomProps(n)
                props["ring_count"] = ring_info.NumAtomRings(idx)
                bond = mol.GetBondBetweenAtoms(hop1_idx, idx)
                props["bt_to_parent"] = _bondTypeCode(bond)
                props["parent_idx"] = hop1_idx
                hop2[idx] = props

    # Hop3 atoms: neighbors of hop2, excluding center, hop1, hop2
    hop3_idx_set = set()
    hop3 = {}
    for hop2_idx in hop2_idx_set:
        hop2_atom = mol.GetAtomWithIdx(hop2_idx)
        for n in hop2_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop1_idx_set or idx in hop2_idx_set or idx == center_idx:
                continue
            if idx not in hop3_idx_set:
                hop3_idx_set.add(idx)
                props = _atomProps(n)
                props["ring_count"] = ring_info.NumAtomRings(idx)
                bond = mol.GetBondBetweenAtoms(hop2_idx, idx)
                props["bt_to_parent"] = _bondTypeCode(bond)
                props["parent_idx"] = hop2_idx
                hop3[idx] = props

    # All bonds within hop3 scope
    bonds = []

    # center-hop1 bonds
    for hop1_idx in hop1_idx_set:
        bonds.append((center_idx, hop1_idx, hop1[hop1_idx]["bt_to_center"]))

    # hop1-hop1 bonds
    for hop1_idx in hop1_idx_set:
        for n2 in mol.GetAtomWithIdx(hop1_idx).GetNeighbors():
            idx2 = n2.GetIdx()
            if idx2 in hop1_idx_set and idx2 > hop1_idx:
                bond = mol.GetBondBetweenAtoms(hop1_idx, idx2)
                bonds.append((hop1_idx, idx2, _bondTypeCode(bond)))

    # hop1-hop2 bonds
    for hop2_idx in hop2_idx_set:
        hop2_atom = mol.GetAtomWithIdx(hop2_idx)
        for n in hop2_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop1_idx_set:
                bond = mol.GetBondBetweenAtoms(hop2_idx, idx)
                bonds.append((hop2_idx, idx, _bondTypeCode(bond)))

    # hop2-hop2 bonds
    for hop2_idx in hop2_idx_set:
        hop2_atom = mol.GetAtomWithIdx(hop2_idx)
        for n in hop2_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop2_idx_set and idx > hop2_idx:
                bond = mol.GetBondBetweenAtoms(hop2_idx, idx)
                bonds.append((hop2_idx, idx, _bondTypeCode(bond)))

    # hop2-hop3 bonds
    for hop3_idx in hop3_idx_set:
        hop3_atom = mol.GetAtomWithIdx(hop3_idx)
        for n in hop3_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop2_idx_set:
                bond = mol.GetBondBetweenAtoms(hop3_idx, idx)
                bonds.append((hop3_idx, idx, _bondTypeCode(bond)))

    # hop3-hop3 bonds
    for hop3_idx in hop3_idx_set:
        hop3_atom = mol.GetAtomWithIdx(hop3_idx)
        for n in hop3_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in hop3_idx_set and idx > hop3_idx:
                bond = mol.GetBondBetweenAtoms(hop3_idx, idx)
                bonds.append((hop3_idx, idx, _bondTypeCode(bond)))

    return {
        "center": center_props,
        "hop1": hop1,
        "hop2": hop2,
        "hop3": hop3,
        "bonds": bonds,
    }


# ── Graph Hashing ────────────────────────────────────────────────────────────────


def _canonicalizeSubgraph(subgraph: dict, center_idx: int, include_hop3: bool = False) -> dict:
    """Canonicalize a subgraph for hashing.

    Assigns canonical 0-based indices to atoms based on sorted properties.

    Args:
        subgraph: Graph dict with center, hop1/hop2 (and optionally hop3).
        center_idx: Original RDKit index of the center atom.
        include_hop3: If True, also canonicalize hop3 atoms.

    Returns:
        Canonicalized subgraph with remapped indices.
    """
    canonical_atoms = []

    # Center is always index 0
    canonical_atoms.append((0, center_idx, subgraph["center"]))

    # Collect hop1 atoms with their original indices
    hop1_items = list(subgraph["hop1"].items())

    # Sort hop1 atoms by their properties for canonical ordering
    def hop1_sort_key(item):
        orig_idx, props = item
        return (props["z"], props["fc"], props["ar"], props["deg"], props["h"], props["ring"], props["bt_to_center"])

    hop1_items.sort(key=hop1_sort_key)

    hop1_remap = {}  # original_idx -> canonical_idx
    for canon_idx, (orig_idx, props) in enumerate(hop1_items, start=1):
        hop1_remap[orig_idx] = canon_idx
        canonical_atoms.append((canon_idx, orig_idx, props))

    # Collect hop2 atoms if present
    hop2_items = list(subgraph.get("hop2", {}).items())

    def hop2_sort_key(item):
        orig_idx, props = item
        return (props["z"], props["fc"], props["ar"], props["deg"], props["h"], props["ring"], props["bt_to_parent"], props["parent_idx"])

    hop2_items.sort(key=hop2_sort_key)

    hop2_remap = {}  # original_idx -> canonical_idx
    for canon_idx, (orig_idx, props) in enumerate(hop2_items, start=len(hop1_items) + 1):
        hop2_remap[orig_idx] = canon_idx
        canonical_atoms.append((canon_idx, orig_idx, props))

    # Collect hop3 atoms if present and include_hop3 is True
    if include_hop3:
        hop3_items = list(subgraph.get("hop3", {}).items())

        def hop3_sort_key(item):
            orig_idx, props = item
            return (props["z"], props["fc"], props["ar"], props["deg"], props["h"], props["ring"], props["bt_to_parent"], props["parent_idx"])

        hop3_items.sort(key=hop3_sort_key)

        hop3_remap = {}  # original_idx -> canonical_idx
        for canon_idx, (orig_idx, props) in enumerate(hop3_items, start=len(hop1_items) + len(hop2_items) + 1):
            hop3_remap[orig_idx] = canon_idx
            canonical_atoms.append((canon_idx, orig_idx, props))
    else:
        hop3_items = []
        hop3_remap = {}

    # Build canonical representation
    canon_hop1 = {}
    for orig_idx, props in hop1_items:
        canon_idx = hop1_remap[orig_idx]
        canon_hop1[canon_idx] = {k: v for k, v in props.items() if k != "bt_to_center"}
        canon_hop1[canon_idx]["bt_to_center"] = props["bt_to_center"]

    canon_hop2 = {}
    for orig_idx, props in hop2_items:
        canon_idx = hop2_remap[orig_idx]
        canon_hop2[canon_idx] = {
            k: v for k, v in props.items()
            if k not in ("bt_to_parent", "parent_idx")
        }
        canon_hop2[canon_idx]["bt_to_parent"] = props["bt_to_parent"]
        canon_hop2[canon_idx]["parent_idx"] = hop1_remap[props["parent_idx"]]

    canon_hop3 = {}
    if include_hop3:
        for orig_idx, props in hop3_items:
            canon_idx = hop3_remap[orig_idx]
            canon_hop3[canon_idx] = {
                k: v for k, v in props.items()
                if k not in ("bt_to_parent", "parent_idx")
            }
            canon_hop3[canon_idx]["bt_to_parent"] = props["bt_to_parent"]
            canon_hop3[canon_idx]["parent_idx"] = hop2_remap[props["parent_idx"]]

    # Remap bonds to canonical indices
    canon_bonds = []
    for idx1, idx2, bt in subgraph["bonds"]:
        if idx1 == center_idx:
            c1 = 0
        elif idx1 in hop1_remap:
            c1 = hop1_remap[idx1]
        elif idx1 in hop2_remap:
            c1 = hop2_remap[idx1]
        else:
            c1 = hop3_remap[idx1]

        if idx2 == center_idx:
            c2 = 0
        elif idx2 in hop1_remap:
            c2 = hop1_remap[idx2]
        elif idx2 in hop2_remap:
            c2 = hop2_remap[idx2]
        else:
            c2 = hop3_remap[idx2]

        if c1 > c2:
            c1, c2 = c2, c1
        canon_bonds.append((c1, c2, bt))

    canon_bonds.sort()

    result = {
        "center": canonical_atoms[0][2],
        "hop1": canon_hop1,
        "hop2": canon_hop2,
        "bonds": canon_bonds,
    }
    if include_hop3:
        result["hop3"] = canon_hop3

    return result


def encodeHop1Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    """Encode hop1 subgraph as SHA-256 hash.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        64-character SHA-256 hex string.
    """
    subgraph = getHop1Subgraph(mol, atom)
    center_idx = atom.GetIdx()
    canonical = _canonicalizeSubgraph(subgraph, center_idx)
    normalized = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode()).hexdigest()


def encodeHop2Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    """Encode hop2 subgraph as SHA-256 hash.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        64-character SHA-256 hex string.
    """
    subgraph = getHop2Subgraph(mol, atom)
    center_idx = atom.GetIdx()
    canonical = _canonicalizeSubgraph(subgraph, center_idx, include_hop3=False)
    normalized = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode()).hexdigest()


def encodeHop3Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    """Encode hop3 subgraph as SHA-256 hash.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        64-character SHA-256 hex string.
    """
    subgraph = getHop3Subgraph(mol, atom)
    center_idx = atom.GetIdx()
    canonical = _canonicalizeSubgraph(subgraph, center_idx, include_hop3=True)
    normalized = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode()).hexdigest()


def computeGraphHopKeys(mol: rdchem.Mol, atom: rdchem.Atom) -> tuple[str, str, str, str]:
    """Compute hop3/hop2/hop1/hop0 keys using graph-based encoding.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        Tuple of (hop3Key, hop2Key, hop1Key, hop0Key).
    """
    # hop3Key uses full hop3 graph
    hop3_key = encodeHop3Graph(mol, atom)

    # hop2Key uses hop2 subgraph
    hop2_key = encodeHop2Graph(mol, atom)

    # hop1Key uses hop1 subgraph
    hop1_key = encodeHop1Graph(mol, atom)

    # hop0Key uses only center atom properties
    center_props = _atomProps(atom)
    ring_info = mol.GetRingInfo()
    center_props["ring_count"] = ring_info.NumAtomRings(atom.GetIdx())
    hop0_key = hashlib.sha256(
        json.dumps(center_props, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return hop3_key, hop2_key, hop1_key, hop0_key
