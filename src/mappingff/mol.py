"""Molecular structure parsing using RDKit.

This module provides the MolReader class for parsing .mol and .pdb files
and extracting molecular graph information including atoms, bonds, and
3D coordinates. It also provides functions for computing chemical environment
descriptors at three levels of granularity (hop0, hop1, hop2).

The three-level environment descriptor:
    - hop0: Center atom properties + neighbor signatures (no 2nd-order neighbors)
    - hop1: hop0 + detailed hop1 shell (first neighbors of neighbors)
    - hop2: hop1 + detailed hop2 shell (second neighbors of neighbors)

Each environment is encoded as a SHA-256 hex string via encode.py.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import SanitizeMol, rdchem, rdDetermineBonds

from mappingff import encode
from mappingff.lmp import parse_lammps
from mappingff.lmp2rdkitmol import lmp_to_rdkit_mol


class MolReader:
    """Parse .mol/.pdb files and extract molecular graph.

    Uses RDKit to parse molecular structure files and provides methods to
    access atoms, bonds, coordinates, and compute chemical environment
    descriptors at multiple granularities.

    Attributes:
        _path: Path to the molecule file.
        _mol: RDKit Mol object.
    """

    def __init__(self, path: Path):
        """Initialize MolReader with a molecule file.

        Args:
            path: Path to .mol or .pdb file.
        """
        self._path = Path(path)
        self._mol: rdchem.Mol
        self._parse()

    def _parse(self) -> None:
        """Parse the molecule file using RDKit.

        Supports both .mol (MDL Molfile) and .pdb (PDB format) files.
        For PDB files, DetermineBondOrders is used to infer bond orders
        from geometry since PDB CONECT records don't encode bond types.

        Sanitization is performed after parsing to ensure implicit valence
        and other computed properties are available.
        """
        # Determine file type by extension and parse accordingly
        suffix = self._path.suffix.lower()
        if suffix == ".mol":
            text = self._path.read_text(encoding="utf-8")
            # Remove OBJ3D blocks which newer RDKit versions reject
            lines = text.splitlines()
            skip = False
            filtered_lines = []
            for line in lines:
                if "M  V30 BEGIN OBJ3D" in line:
                    skip = True
                    continue
                if "M  V30 END OBJ3D" in line:
                    skip = False
                    continue
                if not skip:
                    filtered_lines.append(line)
            text = "\n".join(filtered_lines)
            mol = Chem.MolFromMolBlock(text, sanitize=True, removeHs=False)
            if mol is None:
                raise ValueError(f"Failed to parse MOL file: {self._path}")
            self._mol = mol

        elif suffix == ".pdb":
            pdb_text = self._path.read_text(encoding="utf-8", errors="replace")
            has_conect = "CONECT" in pdb_text
            mol = Chem.MolFromPDBBlock(
                pdb_text,
                sanitize=False,
                removeHs=False,
                proximityBonding=not has_conect,
            )
            if mol is None:
                raise ValueError(f"Failed to parse PDB file: {self._path}")
            rdDetermineBonds.DetermineBondOrders(mol, charge=0)
            SanitizeMol(mol)
            self._mol = mol

        elif suffix == ".lmp" or suffix == ".lammps" or suffix == ".data":
            lmp_data = parse_lammps(self._path)
            self._mol = lmp_to_rdkit_mol(lmp_data)
        else:
            raise ValueError(f"Unsupported file type: {self._path}")

    @property
    def mol(self) -> rdchem.Mol:
        """Get the RDKit Mol object."""
        return self._mol

    def get_atoms(self) -> list[dict]:
        """Get all atoms as a list of dictionaries.

        Each dictionary contains properties for one atom including index,
        element symbol, atomic number, charge, degree, hybridization,
        ring membership, and hydrogen count.

        Returns:
            List of atom records, each with keys:
                - idx: 1-based atom index
                - symbol: element symbol (e.g., 'C', 'O')
                - atomic_num: atomic number
                - formal_charge: formal charge
                - degree: total degree (number of bonds)
                - hybridization: hybridization state as string
                - in_ring: whether atom is in a ring (0 or 1)
                - ring_count: number of rings the atom participates in
                - total_hs: total number of implicit hydrogens
                - aromatic: whether atom is aromatic (0 or 1)
                - chiral_tag: chiral tag as string
        """
        ring_info = self._mol.GetRingInfo()
        atoms = []
        for i, rd_atom in enumerate(self._mol.GetAtoms()):
            # Count how many rings (of size 3-8) this atom participates in
            atom_ring_count = sum(
                1
                for ring_size in range(3, 9)
                if ring_info.IsAtomInRingOfSize(i, ring_size)
            )
            atoms.append(
                {
                    "idx": i + 1,  # 1-based indexing for external use
                    "symbol": rd_atom.GetSymbol(),
                    "atomic_num": rd_atom.GetAtomicNum(),
                    "formal_charge": rd_atom.GetFormalCharge(),
                    "degree": rd_atom.GetTotalDegree(),
                    "hybridization": str(rd_atom.GetHybridization()),
                    "in_ring": int(rd_atom.IsInRing()),
                    "ring_count": atom_ring_count,
                    "total_hs": rd_atom.GetTotalNumHs(includeNeighbors=True),
                    "aromatic": int(rd_atom.GetIsAromatic()),
                    "chiral_tag": str(rd_atom.GetChiralTag()),
                }
            )
        return atoms

    def get_bonds(self) -> list[dict]:
        """Get all bonds as a list of dictionaries.

        Returns:
            List of bond records, each with keys:
                - idx: 1-based bond index
                - a1: 1-based index of first atom
                - a2: 1-based index of second atom
                - bond_type: type as string (e.g., 'SINGLE', 'DOUBLE', 'AROMATIC')
                - aromatic: whether bond is aromatic (0 or 1)
        """
        bonds = []
        for i, rd_bond in enumerate(self._mol.GetBonds()):
            bonds.append(
                {
                    "idx": i + 1,
                    "a1": rd_bond.GetBeginAtomIdx() + 1,  # 1-based
                    "a2": rd_bond.GetEndAtomIdx() + 1,  # 1-based
                    "bond_type": str(rd_bond.GetBondType()),
                    "aromatic": int(rd_bond.GetIsAromatic()),
                }
            )
        return bonds

    def get_coords(self) -> dict[int, tuple[float, float, float]]:
        """Get 3D coordinates for all atoms.

        Returns:
            Dictionary mapping 1-based atom index to (x, y, z) tuple.
        """
        coords = {}
        conf = self._mol.GetConformer()
        for i in range(self._mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords[i + 1] = (pos.x, pos.y, pos.z)
        return coords

    def get_angles(self) -> list[dict]:
        """Get all valence angles as a list of dictionaries.

        An angle is defined by three atoms a1-a2-a3 where a2 is the center.
        Derived by finding all atoms with at least 2 neighbors.

        Returns:
            List of angle records, each with keys:
                - idx: 1-based angle index
                - a1: 1-based index of first atom
                - a2: 1-based index of center atom
                - a3: 1-based index of third atom
        """
        angles = []
        angle_idx = 1
        for center_idx in range(self._mol.GetNumAtoms()):
            center_atom = self._mol.GetAtomWithIdx(center_idx)
            neighbors = center_atom.GetNeighbors()
            if len(neighbors) < 2:
                continue
            # Get all pairs of neighbors to form angles
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    a1 = neighbors[i].GetIdx() + 1
                    a2 = center_idx + 1
                    a3 = neighbors[j].GetIdx() + 1
                    angles.append(
                        {
                            "idx": angle_idx,
                            "a1": a1,
                            "a2": a2,
                            "a3": a3,
                        }
                    )
                    angle_idx += 1
        return angles

    def get_dihedrals(self) -> list[dict]:
        """Get all dihedral angles as a list of dictionaries.

        A dihedral is defined by four atoms a1-a2-a3-a4 where the bonds
        are a1-a2, a2-a3, a3-a4. Derived from all paths of length 3.

        Returns:
            List of dihedral records, each with keys:
                - idx: 1-based dihedral index
                - a1: 1-based index of first atom
                - a2: 1-based index of second atom
                - a3: 1-based index of third atom
                - a4: 1-based index of fourth atom
        """
        dihedrals = []
        dih_idx = 1
        visited = set()
        for bond in self._mol.GetBonds():
            a2_idx = bond.GetBeginAtomIdx()
            a3_idx = bond.GetEndAtomIdx()

            a2_atom = self._mol.GetAtomWithIdx(a2_idx)
            a3_atom = self._mol.GetAtomWithIdx(a3_idx)

            # For a2, get neighbors excluding a3
            for a1_neighbor in a2_atom.GetNeighbors():
                a1_idx = a1_neighbor.GetIdx()
                if a1_idx == a3_idx:
                    continue
                # For a3, get neighbors excluding a2
                for a4_neighbor in a3_atom.GetNeighbors():
                    a4_idx = a4_neighbor.GetIdx()
                    if a4_idx == a2_idx:
                        continue
                    # Create canonical dihedral key (a1 < a4 to avoid duplicates)
                    if a1_idx < a4_idx:
                        key = (a1_idx, a2_idx, a3_idx, a4_idx)
                    else:
                        key = (a4_idx, a3_idx, a2_idx, a1_idx)
                    if key in visited:
                        continue
                    visited.add(key)
                    dihedrals.append(
                        {
                            "idx": dih_idx,
                            "a1": a1_idx + 1,
                            "a2": a2_idx + 1,
                            "a3": a3_idx + 1,
                            "a4": a4_idx + 1,
                        }
                    )
                    dih_idx += 1
        return dihedrals

    def get_impropers(self) -> list[dict]:
        """Get all improper angles as a list of dictionaries.

        An improper is defined by a central atom and three substituents.
        For each atom with 3 or more neighbors, all combinations of 3
        neighbors are generated as impropers.

        Returns:
            List of improper records, each with keys:
                - idx: 1-based improper index
                - a1: 1-based index of center atom
                - a2: 1-based index of second atom
                - a3: 1-based index of third atom
                - a4: 1-based index of fourth atom
        """
        from itertools import combinations

        impropers = []
        imp_idx = 1
        visited = set()

        # For each atom with 3 or more neighbors, generate all combinations of 3 neighbors
        for atom in self._mol.GetAtoms():
            center_idx = atom.GetIdx()
            neighbor_indices = sorted([n.GetIdx() for n in atom.GetNeighbors()])

            if len(neighbor_indices) < 3:
                continue

            # Generate all combinations of 3 neighbors
            for a, b, c in combinations(neighbor_indices, 3):
                # Create canonical key to avoid duplicates
                key = (center_idx, a, b, c)
                if key in visited:
                    continue
                visited.add(key)

                impropers.append(
                    {
                        "idx": imp_idx,
                        "a1": center_idx + 1,
                        "a2": a + 1,
                        "a3": b + 1,
                        "a4": c + 1,
                    }
                )
                imp_idx += 1

        return impropers


# ── Convenience Re-exports from encode.py ──────────────────────────────────────


def compute_hop_keys(mol_reader: MolReader, atom_idx: int) -> tuple[str, str, str, str]:
    """Compute hop3, hop2, hop1, and hop0 keys for an atom using graph-based encoding.

    Args:
        mol_reader: MolReader instance with loaded molecule.
        atom_idx: 1-based atom index.

    Returns:
        Tuple of (hop3_key, hop2_key, hop1_key, hop0_key), all 64-char hex strings.
    """
    rd_atom = mol_reader._mol.GetAtomWithIdx(atom_idx - 1)
    return encode.compute_graph_hop_keys(mol_reader._mol, rd_atom)
