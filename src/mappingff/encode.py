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


def _atom_props(atom: rdchem.Atom) -> dict:
    """Get canonical properties dict for an atom."""
    return {
        "z": atom.GetAtomicNum(),
        "fc": atom.GetFormalCharge(),
        "ar": int(atom.GetIsAromatic()),
        "deg": atom.GetTotalDegree(),
        "h": atom.GetTotalNumHs(includeNeighbors=True),
        "ring": int(atom.IsInRing()),
    }


def _props_with_ring(atom: rdchem.Atom, mol: rdchem.Mol) -> dict:
    props = _atom_props(atom)
    props["ring_count"] = mol.GetRingInfo().NumAtomRings(atom.GetIdx())
    return props


def _collect_bonds(atom_idx_set: set, mol: rdchem.Mol, exclude_self: bool) -> list:
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
        bond = mol.GetBondBetweenAtoms(center_idx, neighbor.GetIdx())
        bt_code = _bond_type_code(bond)
        sig = f"{neighbor.GetAtomicNum()}:{bt_code}:{neighbor.GetFormalCharge()}"
        neighbor_sig.append(sig)
        bond_kinds.append(bt_code)

    neighbor_sig.sort()
    bond_kinds.sort()

    center_props["neighbor_sig"] = neighbor_sig
    center_props["bond_kinds"] = bond_kinds

    return _canonicalize_subgraph(
        {"center": center_props, "hop1": {}, "hop2": {}, "hop3": {}, "bonds": []},
        center_idx,
    )


def _build_hop_level(
    mol: rdchem.Mol,
    parent_idx_set: set,
    exclude: set,
) -> tuple[set, dict]:
    """Build one hop level from parent atom indices.

    Args:
        mol: RDKit molecule.
        parent_idx_set: Set of parent atom indices.
        exclude: Indices to skip (center + previous levels).

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
            props = _props_with_ring(n, mol)
            bond = mol.GetBondBetweenAtoms(p_idx, idx)
            props["bt_to_parent"] = _bond_type_code(bond)
            props["parent_idx"] = p_idx
            atom_dict[idx] = props
    return idx_set, atom_dict


def _center_props_with_ring(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    center_props = _atom_props(atom)
    center_props["ring_count"] = mol.GetRingInfo().NumAtomRings(center_idx)
    return center_props


def _build_hop1_shell(
    mol: rdchem.Mol, atom: rdchem.Atom
) -> tuple[set[int], dict[int, dict]]:
    center_idx = atom.GetIdx()
    hop1_idx_set = {n.GetIdx() for n in atom.GetNeighbors()}
    hop1 = {}
    for neighbor in atom.GetNeighbors():
        idx = neighbor.GetIdx()
        props = _props_with_ring(neighbor, mol)
        bond = mol.GetBondBetweenAtoms(center_idx, idx)
        props["bt_to_center"] = _bond_type_code(bond)
        hop1[idx] = props
    return hop1_idx_set, hop1


def _build_subgraph(mol: rdchem.Mol, atom: rdchem.Atom, max_hop: int) -> dict:
    center_idx = atom.GetIdx()
    center_props = _center_props_with_ring(mol, atom)

    hop1_idx_set, hop1 = _build_hop1_shell(mol, atom)
    hop2_idx_set: set[int] = set()
    hop2: dict[int, dict] = {}
    hop3_idx_set: set[int] = set()
    hop3: dict[int, dict] = {}

    if max_hop >= 2:
        hop2_idx_set, hop2 = _build_hop_level(
            mol, hop1_idx_set, hop1_idx_set | {center_idx}
        )
    if max_hop >= 3:
        hop3_idx_set, hop3 = _build_hop_level(
            mol, hop2_idx_set, hop1_idx_set | hop2_idx_set | {center_idx}
        )

    bonds = [(center_idx, h1, hop1[h1]["bt_to_center"]) for h1 in hop1_idx_set]
    if max_hop >= 3:
        bonds += _collect_cross_level_bonds(mol, hop2_idx_set, hop1_idx_set)
        bonds += _collect_cross_level_bonds(mol, hop3_idx_set, hop2_idx_set)
    bonds += _collect_bonds(hop1_idx_set, mol, exclude_self=True)
    if max_hop >= 2:
        bonds += _collect_bonds(hop2_idx_set, mol, exclude_self=True)
    if max_hop >= 3:
        bonds += _collect_bonds(hop3_idx_set, mol, exclude_self=True)

    return {
        "center": center_props,
        "hop1": hop1,
        "hop2": hop2 if max_hop >= 2 else {},
        "hop3": hop3 if max_hop >= 3 else {},
        "bonds": bonds,
    }


def _get_hop1_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    subgraph = _build_subgraph(mol, atom, max_hop=1)
    return _canonicalize_subgraph(subgraph, center_idx)


def _get_hop2_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    subgraph = _build_subgraph(mol, atom, max_hop=2)
    return _canonicalize_subgraph(subgraph, center_idx)


def get_hop3_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    center_idx = atom.GetIdx()
    subgraph = _build_subgraph(mol, atom, max_hop=3)
    return _canonicalize_subgraph(subgraph, center_idx, include_hop3=True)


# ── Graph Hashing ────────────────────────────────────────────────────────────────


def _collect_cross_level_bonds(
    mol: rdchem.Mol, src_idx_set: set[int], dst_idx_set: set[int]
) -> list[tuple[int, int, str]]:
    bonds = []
    for src_idx in src_idx_set:
        for neighbor in mol.GetAtomWithIdx(src_idx).GetNeighbors():
            dst_idx = neighbor.GetIdx()
            if dst_idx in dst_idx_set:
                bond = mol.GetBondBetweenAtoms(src_idx, dst_idx)
                bonds.append((src_idx, dst_idx, _bond_type_code(bond)))
    return bonds


def _hop_sort_key(props: dict, bond_key: str) -> tuple[int, ...]:
    key: tuple[int, ...] = (
        props["z"],
        props["fc"],
        props["ar"],
        props["deg"],
        props["h"],
        props["ring"],
        props[bond_key],
    )
    if "parent_idx" in props:
        key = key + (props["parent_idx"],)
    return key


def _canonicalize_hop(
    items: list[tuple[int, dict]],
    start_idx: int,
    bond_key: str,
    parent_remap: dict[int, int] | None = None,
) -> tuple[dict[int, int], dict[int, dict]]:
    items.sort(key=lambda item: _hop_sort_key(item[1], bond_key))
    remap = {}
    canon_hop = {}
    for canon_idx, (orig_idx, props) in enumerate(items, start=start_idx):
        remap[orig_idx] = canon_idx
        canon_props = {
            k: v for k, v in props.items() if k not in (bond_key, "parent_idx")
        }
        canon_props[bond_key] = props[bond_key]
        if parent_remap is not None:
            canon_props["parent_idx"] = parent_remap[props["parent_idx"]]
        canon_hop[canon_idx] = canon_props
    return remap, canon_hop


def _canonicalize_subgraph(
    subgraph: dict, center_idx: int, include_hop3: bool = False
) -> dict:
    """Canonicalize a subgraph for hashing.

    Assigns canonical 0-based indices to atoms based on sorted properties.

    Args:
        subgraph: Graph dict with center, hop1/hop2 (and optionally hop3).
        center_idx: Original RDKit index of the center atom.
        include_hop3: If True, also canonicalize hop3 atoms.

    Returns:
        Canonicalized subgraph with remapped indices.
    """
    # Collect hop1 atoms with their original indices
    hop1_items = list(subgraph["hop1"].items())

    # Collect hop2 atoms if present
    hop2_items = list(subgraph.get("hop2", {}).items())

    # Collect hop3 atoms if present and include_hop3 is True
    if include_hop3:
        hop3_items = list(subgraph.get("hop3", {}).items())
    else:
        hop3_items = []

    hop1_remap, canon_hop1 = _canonicalize_hop(
        hop1_items, start_idx=1, bond_key="bt_to_center"
    )
    hop2_remap, canon_hop2 = _canonicalize_hop(
        hop2_items,
        start_idx=len(hop1_items) + 1,
        bond_key="bt_to_parent",
        parent_remap=hop1_remap,
    )
    if include_hop3:
        hop3_remap, canon_hop3 = _canonicalize_hop(
            hop3_items,
            start_idx=len(hop1_items) + len(hop2_items) + 1,
            bond_key="bt_to_parent",
            parent_remap=hop2_remap,
        )
    else:
        hop3_remap = {}
        canon_hop3 = {}

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
        "center": subgraph["center"],
        "hop1": canon_hop1,
        "hop2": canon_hop2,
        "bonds": canon_bonds,
    }
    if include_hop3:
        result["hop3"] = canon_hop3

    return result


def _encode_hop1_graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = _get_hop1_subgraph(mol, atom)
    return hashlib.sha256(
        json.dumps(subgraph, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _encode_hop2_graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    subgraph = _get_hop2_subgraph(mol, atom)
    return hashlib.sha256(
        json.dumps(subgraph, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _encode_hop3_graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    canonical = get_hop3_subgraph(mol, atom)  # already canonicalized
    return hashlib.sha256(
        json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _encode_hop0_graph(mol: rdchem.Mol, atom: rdchem.Atom) -> str:
    canonical = get_hop0_subgraph(mol, atom)  # already canonicalized
    return hashlib.sha256(
        json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def compute_graph_hop_keys(
    mol: rdchem.Mol, atom: rdchem.Atom
) -> tuple[str, str, str, str]:
    """Compute hop3/hop2/hop1/hop0 keys using graph-based encoding.

    Args:
        mol: RDKit molecule.
        atom: Center atom.

    Returns:
        Tuple of (hop3_key, hop2_key, hop1_key, hop0_key).
    """
    hop3_key = _encode_hop3_graph(mol, atom)
    hop2_key = _encode_hop2_graph(mol, atom)
    hop1_key = _encode_hop1_graph(mol, atom)
    hop0_key = _encode_hop0_graph(mol, atom)

    return hop3_key, hop2_key, hop1_key, hop0_key
