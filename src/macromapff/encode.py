"""Environment key encoding for MacroMapFF.

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
    from rdkit import Chem
    from rdkit.Chem import rdchem


def encodeAtomEnvHop0(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Compute hop0 level environment for an atom.

    Hop0 is the coarsest environment descriptor. It captures:
        - Center atom properties (z, charge, degree, hybridization, ring membership)
        - Neighbor signatures: "atomicNum:bondType:formalCharge" for each neighbor
        - Bond kind list: sorted list of S/D/T/A/U codes for each bond

    This level is used for bonded parameter matching (bond/angle/dihedral/improper)
    because it captures the immediate chemical environment without 2nd-order details.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object (1-based index in calling code).

    Returns:
        Dictionary with keys:
            - z: atomic number
            - formal_charge: formal charge
            - degree: total degree
            - hybridization: hybridization as string
            - in_ring: 0 or 1
            - ring_count: number of rings atom participates in
            - total_hs: number of implicit hydrogens
            - neighbor_sig: sorted list of "z:bondType:charge" strings
            - bond_kinds: sorted list of S/D/T/A/U codes
    """
    ring_info = mol.GetRingInfo()
    atom_idx = atom.GetIdx()
    ring_count = ring_info.NumAtomRings(atom_idx)

    return {
        "z": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "aromatic": int(atom.GetIsAromatic()),
        "degree": atom.GetDegree(),
        "hybridization": str(atom.GetHybridization()),
        "in_ring": int(atom.IsInRing()),
        "ring_count": ring_count,
        "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
        "neighbor_sig": _neighborSignature(mol, atom),
        "bond_kinds": _bondKindList(mol, atom),
    }


def encodeAtomEnvHop1(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Compute hop1 level environment for an atom.

    Hop1 extends hop0 by adding the hop1_shell, which contains detailed
    properties of first neighbors (their element, aromaticity, degree, etc.).

    This level is used for the first fallback step in atom typing when
    no exact hop2 match is found.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object.

    Returns:
        Dictionary with all hop0 keys plus:
            - hop1_shell: list of neighbor descriptors, each with z, ar, deg, fc, h, ring, bt
    """
    env = encodeAtomEnvHop0(mol, atom)
    env["hop1_shell"] = _shellNeighbors(mol, atom, depth=1)
    return env


def encodeAtomEnvHop2(mol: rdchem.Mol, atom: rdchem.Atom) -> dict:
    """Compute hop2 level environment for an atom.

    Hop2 extends hop1 by adding the hop2_shell, which contains detailed
    properties of second neighbors (neighbors of neighbors, excluding the
    center and hop1 atoms).

    This is the finest-grained environment descriptor and provides the
    most specific atom type identification when an exact match exists.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object.

    Returns:
        Dictionary with all hop1 keys plus:
            - hop2_shell: list of 2nd-order neighbor descriptors
    """
    env = encodeAtomEnvHop1(mol, atom)
    env["hop2_shell"] = _shellNeighbors(mol, atom, depth=2)
    return env


def computeHopKeys(envHop2: dict) -> tuple[str, str, str]:
    """Compute hop2, hop1, and hop0 keys from a hop2 environment dict.

    Derives hop1 key by stripping hop2_shell from hop2 env.
    Derives hop0 key by stripping hop1_shell from hop1 env.
    Each key is a SHA-256 encoding of the respective environment dict.

    Args:
        envHop2: Full hop2 environment dictionary.

    Returns:
        Tuple of (hop2Key, hop1Key, hop0Key), all 64-character hex strings.
    """
    hop1Env = _stripHop2(envHop2)
    hop0Env = _stripHop1(hop1Env)

    hop2Key = encodeEnvKey(envHop2)
    hop1Key = encodeEnvKey(hop1Env)
    hop0Key = encodeEnvKey(hop0Env)

    return hop2Key, hop1Key, hop0Key


# ── Helper Functions ────────────────────────────────────────────────────────────


def _neighborSignature(mol: rdchem.Mol, atom: rdchem.Atom) -> list[str]:
    """Generate neighbor signature list for an atom.

    Each neighbor contributes a signature string in the format:
    "atomicNum:bondType:aromaticFlag"

    The bond type uses single-character codes: S, D, T, A, U
    to match the legacy format used in the reference implementation.

    The aromaticFlag is 1 if the neighbor is aromatic, 0 otherwise.

    The signatures are sorted to ensure canonical ordering regardless of
    the order in which neighbors are enumerated by RDKit.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object (center atom).

    Returns:
        Sorted list of neighbor signature strings.
    """
    from rdkit import Chem
    sigs = []
    for neighbor in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE:
            bond_code = "S"
        elif bt == Chem.rdchem.BondType.DOUBLE:
            bond_code = "D"
        elif bt == Chem.rdchem.BondType.TRIPLE:
            bond_code = "T"
        elif bt == Chem.rdchem.BondType.AROMATIC:
            bond_code = "A"
        else:
            bond_code = "U"
        sig = f"{neighbor.GetAtomicNum()}:{bond_code}:{int(neighbor.GetIsAromatic())}"
        sigs.append(sig)
    sigs.sort()
    return sigs


def _bondKindList(mol: rdchem.Mol, atom: rdchem.Atom) -> list[str]:
    """Get sorted list of bond kind codes for an atom.

    Each connected bond contributes a single character code:
        - S: SINGLE
        - D: DOUBLE
        - T: TRIPLE
        - A: AROMATIC
        - U: UNKNOWN

    The list is sorted to ensure canonical ordering.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object.

    Returns:
        Sorted list of bond kind codes (S, D, T, A, U).
    """
    from rdkit import Chem
    kinds = []
    for neighbor in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE:
            kinds.append("S")
        elif bt == Chem.rdchem.BondType.DOUBLE:
            kinds.append("D")
        elif bt == Chem.rdchem.BondType.TRIPLE:
            kinds.append("T")
        elif bt == Chem.rdchem.BondType.AROMATIC:
            kinds.append("A")
        else:
            kinds.append("U")
    kinds.sort()
    return kinds


def _shellNeighbors(mol: rdchem.Mol, atom: rdchem.Atom, depth: int) -> list[dict]:
    """Get neighboring atoms at specified shell depth.

    For depth=1 (hop1): Returns all first neighbors of the center atom.
    For depth=2 (hop2): Returns neighbors of first neighbors, excluding
                        the center and hop1 atoms.

    Each neighbor is described by a dictionary with:
        z, ar (aromatic), deg (degree), fc (formal charge),
        h (total Hs), ring, bt (bond type)

    The result is sorted by (z, bt, ar, deg) for canonical ordering.

    Args:
        mol: RDKit molecule object.
        atom: RDKit atom object (center atom).
        depth: Shell depth (1 or 2).

    Returns:
        Sorted list of neighbor descriptor dictionaries.

    Raises:
        ValueError: If depth is not 1 or 2.
    """
    if depth == 1:
        shell = []
        for neighbor in atom.GetNeighbors():
            entry = {
                "z": neighbor.GetAtomicNum(),
                "fc": neighbor.GetFormalCharge(),
                "ar": int(neighbor.GetIsAromatic()),
                "deg": neighbor.GetTotalDegree(),
                "h": neighbor.GetTotalNumHs(includeNeighbors=True),
                "ring": int(neighbor.IsInRing()),
            }
            shell.append(entry)
        shell.sort(key=lambda x: (x["z"], x["fc"], x["ar"], x["deg"], x["h"], x["ring"]))
        return shell

    elif depth == 2:
        # Collect hop1 atom indices for exclusion
        hop1_neighbors = {n.GetIdx() for n in atom.GetNeighbors()}
        shell = []
        for hop1_neighbor in atom.GetNeighbors():
            for hop2_neighbor in hop1_neighbor.GetNeighbors():
                # Skip back to center or other hop1 atoms
                if hop2_neighbor.GetIdx() in hop1_neighbors:
                    continue
                if hop2_neighbor.GetIdx() == atom.GetIdx():
                    continue
                entry = {
                    "z": hop2_neighbor.GetAtomicNum(),
                    "fc": hop2_neighbor.GetFormalCharge(),
                    "ar": int(hop2_neighbor.GetIsAromatic()),
                    "deg": hop2_neighbor.GetTotalDegree(),
                    "h": hop2_neighbor.GetTotalNumHs(includeNeighbors=True),
                    "ring": int(hop2_neighbor.IsInRing()),
                }
                shell.append(entry)
        shell.sort(key=lambda x: (x["z"], x["fc"], x["ar"], x["deg"], x["h"], x["ring"]))
        return shell

    else:
        raise ValueError(f"Unsupported depth: {depth}")


def encodeEnvKey(envDict: dict) -> str:
    """Encode an environment dictionary as a SHA-256 hex string.

    The dictionary is normalized by JSON serialization with sorted keys
    and no whitespace, ensuring the same environment always produces
    the same hash regardless of dict ordering.

    Args:
        envDict: Environment dictionary to encode.

    Returns:
        64-character SHA-256 hexadecimal string.
    """
    normalized = json.dumps(envDict, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode()).hexdigest()


def _stripHop2(envHop2: dict) -> dict:
    """Strip hop2_shell from a hop2 environment dict to get hop1 env.

    Args:
        envHop2: Full hop2 environment dictionary.

    Returns:
        Dictionary with hop2_shell key removed.
    """
    return {k: v for k, v in envHop2.items() if k != "hop2_shell"}


def _stripHop1(envHop1: dict) -> dict:
    """Strip hop1_shell from a hop1 environment dict to get hop0 env.

    Args:
        envHop1: Hop1 environment dictionary (with hop1_shell).

    Returns:
        Dictionary with hop1_shell key removed.
    """
    return {k: v for k, v in envHop1.items() if k != "hop1_shell"}


def _hop0Env(envHop2: dict) -> dict:
    """Strip both hop1_shell and hop2_shell to get hop0 env.

    This is equivalent to calling _stripHop2 then _stripHop1.

    Args:
        envHop2: Full hop2 environment dictionary.

    Returns:
        Dictionary with both hop1_shell and hop2_shell removed.
    """
    result = _stripHop2(envHop2)
    result = _stripHop1(result)
    return result
