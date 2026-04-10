#!/usr/bin/env python3
"""RDKit-based atom environment feature extraction and env-key generation."""

import json

from rdkit import Chem

from macromapff.domain.env_key_codec import ordered_env_key_obj


INTERNAL_HOP_DEPTH = 2


def bond_type_code(bond: Chem.Bond) -> str:
    """Map RDKit bond types to compact one-letter codes."""
    bt = bond.GetBondType()
    if bt == Chem.rdchem.BondType.SINGLE:
        return "S"
    if bt == Chem.rdchem.BondType.DOUBLE:
        return "D"
    if bt == Chem.rdchem.BondType.TRIPLE:
        return "T"
    if bt == Chem.rdchem.BondType.AROMATIC:
        return "A"
    return str(bt)


def precompute_atom_context(mol: Chem.Mol):
    """Precompute reusable atom context such as distance matrix."""
    dist_matrix = Chem.GetDistanceMatrix(mol)

    context = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        context[idx] = {
            "dist_matrix": dist_matrix,
        }

    return context


def hop_shell_signatures(mol: Chem.Mol, atom_idx: int, dist_matrix):
    """Collect sorted neighbor shell signatures up to a hop depth."""
    dist = dist_matrix
    shells = {hop: [] for hop in range(1, INTERNAL_HOP_DEPTH + 1)}
    for idx in range(mol.GetNumAtoms()):
        hop = int(dist[atom_idx][idx])
        if hop < 1 or hop > INTERNAL_HOP_DEPTH:
            continue
        atom = mol.GetAtomWithIdx(idx)
        token = {
            "z": atom.GetAtomicNum(),
            "fc": atom.GetFormalCharge(),
            "ar": int(atom.GetIsAromatic()),
            "deg": atom.GetDegree(),
            "h": atom.GetTotalNumHs(includeNeighbors=True),
            "ring": int(atom.IsInRing()),
        }
        shells[hop].append(token)

    out = {}
    for hop in range(1, INTERNAL_HOP_DEPTH + 1):
        out[f"hop{hop}_shell"] = sorted(
            shells[hop],
            key=lambda x: (
                x["z"],
                x["fc"],
                x["ar"],
                x["deg"],
                x["h"],
                x["ring"],
            ),
        )
    return out


def make_env_key(
    mol: Chem.Mol,
    atom: Chem.Atom,
    atom_ctx: dict,
):
    """Create canonical env-key string and raw feature dict for an atom."""
    atom_idx = atom.GetIdx()
    nbr_sigs = []
    bond_kinds = []

    for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        bcode = bond_type_code(bond)
        bond_kinds.append(bcode)
        nbr_sigs.append(f"{nbr.GetAtomicNum()}:{bcode}:{int(nbr.GetIsAromatic())}")

    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumAtomRings(atom_idx)
    dist_matrix = atom_ctx["dist_matrix"]

    features = {
        "z": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "aromatic": int(atom.GetIsAromatic()),
        "hybridization": str(atom.GetHybridization()),
        "in_ring": int(atom.IsInRing()),
        "ring_count": ring_count,
        "degree": atom.GetDegree(),
        "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
        "neighbor_sig": sorted(nbr_sigs),
        "bond_kinds": sorted(bond_kinds),
    }
    features.update(hop_shell_signatures(mol, atom_idx, dist_matrix=dist_matrix))

    key = json.dumps(
        ordered_env_key_obj(features),
        ensure_ascii=False,
        sort_keys=False,
        separators=(",", ":"),
    )
    return key, features
