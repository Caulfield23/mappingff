"""Environment key encoding for mappingff.

This module provides functions for computing multi-level chemical environment
descriptors (hop0, hop1, hop2, hop3) and encoding them as SHA-256 keys.

The four levels of environment description:
    - hop0: Center atom + neighbor signatures + bond kinds (coarse-grained)
    - hop1: hop0 + first neighbor details (first fallback)
    - hop2: hop1 + second neighbor details (second fallback)
    - hop3: hop2 + third neighbor details (finest classification)

Each environment is encoded as a SHA-256 hash for use as a dictionary key
in the parameter database.
"""

from __future__ import annotations

import hashlib
import json

from rdkit import Chem
from rdkit.Chem import rdchem

_BOND_TYPE_CODE = {
    Chem.rdchem.BondType.SINGLE: "S",
    Chem.rdchem.BondType.DOUBLE: "D",
    Chem.rdchem.BondType.TRIPLE: "T",
    Chem.rdchem.BondType.AROMATIC: "A",
}

def _bond_type_code(bond: rdchem.Bond) -> str:
    return _BOND_TYPE_CODE.get(bond.GetBondType(), "U")


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


def _propsWithRing(atom: rdchem.Atom, mol: rdchem.Mol) -> dict:
    props = _atomProps(atom)
    props["ring_count"] = mol.GetRingInfo().NumAtomRings(atom.GetIdx())
    return props


def _collectBonds(atom_idx_set: set, mol: rdchem.Mol, exclude_self: bool) -> list:
    """Collect bonds among atoms in atom_idx_set via mol.GetAtomWithIdx."""
    bonds = []
    for idx in atom_idx_set:
        atom = mol.GetAtomWithIdx(idx)
        for n in atom.GetNeighbors():
            idx2 = n.GetIdx()
            if idx2 not in atom_idx_set:
                continue
            if exclude_self and idx2 <= idx:
                continue
            bond = mol.GetBondBetweenAtoms(idx, idx2)
            bonds.append((idx, idx2, _bond_type_code(bond)))
    return bonds


def get_hop0_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Extract hop0 subgraph centered on an atom.

    The subgraph includes:
    - center: center atom properties (z, formal_charge, aromatic, hybridization,
              degree, total_hs, in_ring, ring_count)
    - neighbor_sig: list of "element:bond_type:formal_charge" for each neighbor
    - bond_kinds: sorted list of bond types (S, D, A, T) to neighbors

    Args:
        mol: RDKit molecule.
        atom: Center atom (RDKit 0-based index).

    Returns:
        dict with center properties, neighbor_sig, and bond_kinds.
    """
    center_idx = atom.GetIdx()
    ring_info = mol.GetRingInfo()

    center_props = {
        "z": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "aromatic": int(atom.GetIsAromatic()),
        "hybridization": str(atom.GetHybridization()),
        "degree": atom.GetTotalDegree(),
        "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
        "in_ring": int(atom.IsInRing()),
        "ring_count": ring_info.NumAtomRings(center_idx),
    }

    neighbor_sig = []
    bond_kinds = []
    for neighbor in atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        bond = mol.GetBondBetweenAtoms(center_idx, n_idx)
        bt_code = _bond_type_code(bond)
        sig = f"{neighbor.GetAtomicNum()}:{bt_code}:{neighbor.GetFormalCharge()}"
        neighbor_sig.append(sig)
        bond_kinds.append(bt_code)

    neighbor_sig.sort()
    bond_kinds.sort()

    center_props["neighbor_sig"] = neighbor_sig
    center_props["bond_kinds"] = bond_kinds

    return {"center": center_props, "hop1": {}, "hop2": {}, "bonds": []}


def _buildHopLevel(
    mol: rdchem.Mol,
    parent_idx_set: set,
    exclude: set,
    parent_key: str,
) -> tuple[set, dict]:
    """Build one hop level from parent atom indices.

    Args:
        mol: RDKit molecule.
        parent_idx_set: Set of parent atom indices.
        exclude: Indices to skip (center + previous levels).
        parent_key: "bt_to_center" for hop1, "bt_to_parent" for hop2/hop3.

    Returns:
        (idx_set, atom_dict) for this hop level.
    """
    idx_set = set()
    atom_dict = {}
    for p_idx in parent_idx_set:
        parent_atom = mol.GetAtomWithIdx(p_idx)
        for n in parent_atom.GetNeighbors():
            idx = n.GetIdx()
            if idx in exclude or idx in idx_set:
                continue
            idx_set.add(idx)
            props = _propsWithRing(n, mol)
            bond = mol.GetBondBetweenAtoms(p_idx, idx)
            props["bt_to_parent"] = _bond_type_code(bond)
            props["parent_idx"] = p_idx
            atom_dict[idx] = props
    return idx_set, atom_dict


def _getHop1Subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    center_props = _atomProps(atom)
    center_props["ring_count"] = mol.GetRingInfo().NumAtomRings(center_idx)

    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _propsWithRing(n, mol)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bond_type_code(bond)
        hop1[idx] = props

    bonds = [(center_idx, h1, hop1[h1]["bt_to_center"]) for h1 in hop1_idx_set]
    bonds += _collectBonds(hop1_idx_set, mol, exclude_self=True)

    return {"center": center_props, "hop1": hop1, "hop2": {}, "bonds": bonds}


def _getHop2Subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    center_props = _atomProps(atom)
    center_props["ring_count"] = mol.GetRingInfo().NumAtomRings(center_idx)

    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _propsWithRing(n, mol)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bond_type_code(bond)
        hop1[idx] = props

    hop2_idx_set, hop2 = _buildHopLevel(mol, hop1_idx_set, hop1_idx_set | {center_idx}, "bt_to_parent")

    bonds = [(center_idx, h1, hop1[h1]["bt_to_center"]) for h1 in hop1_idx_set]
    bonds += _collectBonds(hop1_idx_set, mol, exclude_self=True)
    bonds += _collectBonds(hop2_idx_set, mol, exclude_self=True)

    return {"center": center_props, "hop1": hop1, "hop2": hop2, "bonds": bonds}


def get_hop3_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    center_props = _atomProps(atom)
    center_props["ring_count"] = mol.GetRingInfo().NumAtomRings(center_idx)

    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for n in atom.GetNeighbors():
        idx = n.GetIdx()
        props = _propsWithRing(n, mol)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bond_type_code(bond)
        hop1[idx] = props

    hop2_idx_set, hop2 = _buildHopLevel(mol, hop1_idx_set, hop1_idx_set | {center_idx}, "bt_to_parent")
    hop3_idx_set, hop3 = _buildHopLevel(mol, hop2_idx_set, hop1_idx_set | hop2_idx_set | {center_idx}, "bt_to_parent")

    bonds = [(center_idx, h1, hop1[h1]["bt_to_center"]) for h1 in hop1_idx_set]
    for h2 in hop2_idx_set:
        for n in mol.GetAtomWithIdx(h2).GetNeighbors():
            if n.GetIdx() in hop1_idx_set:
                bonds.append((h2, n.GetIdx(), _bond_type_code(mol.GetBondBetweenAtoms(h2, n.GetIdx()))))
    for h3 in hop3_idx_set:
        for n in mol.GetAtomWithIdx(h3).GetNeighbors():
            if n.GetIdx() in hop2_idx_set:
                bonds.append((h3, n.GetIdx(), _bond_type_code(mol.GetBondBetweenAtoms(h3, n.GetIdx()))))
    bonds += _collectBonds(hop1_idx_set, mol, exclude_self=True)
    bonds += _collectBonds(hop2_idx_set, mol, exclude_self=True)
    bonds += _collectBonds(hop3_idx_set, mol, exclude_self=True)

    return {"center": center_props, "hop1": hop1, "hop2": hop2, "hop3": hop3, "bonds": bonds}


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


def _encodeHop1Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = _getHop1Subgraph(mol, atom)
    canonical = _canonicalizeSubgraph(subgraph, atom.GetIdx())
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _encodeHop2Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = _getHop2Subgraph(mol, atom)
    canonical = _canonicalizeSubgraph(subgraph, atom.GetIdx(), include_hop3=False)
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _encodeHop3Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = get_hop3_subgraph(mol, atom)
    canonical = _canonicalizeSubgraph(subgraph, atom.GetIdx(), include_hop3=True)
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _encodeHop0Graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = get_hop0_subgraph(mol, atom)
    canonical = _canonicalizeSubgraph(subgraph, atom.GetIdx())
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _computeGraphHopKeys(mol: rdchem.Mol, atom: rdchem.Atom) -> tuple[str, str, str, str]:
    """Compute hop3/hop2/hop1/hop0 keys using graph-based encoding.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        Tuple of (hop3Key, hop2Key, hop1Key, hop0Key).
    """
    hop3_key = _encodeHop3Graph(mol, atom)
    hop2_key = _encodeHop2Graph(mol, atom)
    hop1_key = _encodeHop1Graph(mol, atom)
    hop0_key = _encodeHop0Graph(mol, atom)

    return hop3_key, hop2_key, hop1_key, hop0_key
