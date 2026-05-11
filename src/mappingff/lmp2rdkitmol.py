"""Element mass table and lmp to RDKit mol converter."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

# Load element masses from mass.txt
_MASSES: dict[str, float] = {}


def _load():
    txt_path = Path(__file__).parent / "mass.txt"
    for line in txt_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            _MASSES[parts[0]] = float(parts[1])


_load()


def mass_to_element(mass: float, tolerance: float = 0.1) -> str | None:
    """Find element symbol by mass within tolerance. Returns None if no match."""
    best, best_diff = None, float("inf")
    for elem, em in _MASSES.items():
        d = abs(mass - em)
        if d <= tolerance and d < best_diff:
            best, best_diff = elem, d
    return best


def lmp_to_rdkit_mol(lmpData, tolerance: float = 0.1) -> Chem.Mol:
    """Convert LammpsData to RDKit Mol.

    1. type_id -> element (mass -> element via massToElement)
    2. Build RWMol with atoms+coords, bonds (connectivity only, no bond order)
    3. DetermineBondOrders() infers bond orders from geometry
    4. Return sanitized Mol
    """
    # 1. type_id -> element
    type_to_elem = {}
    for type_id, mass in lmpData.masses:
        elem = mass_to_element(mass, tolerance)
        if elem is None:
            raise ValueError(f"No element found for mass {mass} (type {type_id})")
        type_to_elem[type_id] = elem

    # 2. Build RWMol
    mol = Chem.RWMol()
    for atom_id, mol_tag, type_id, charge, x, y, z in lmpData.atom_records:
        atom = Chem.Atom(type_to_elem[type_id])
        atom.SetFormalCharge(int(round(charge)))
        mol.AddAtom(atom)

    # Add bonds (connectivity only, DetermineBondOrders will infer order)
    for bondId, bondType, a1, a2 in lmpData.bond_records:
        mol.AddBond(a1 - 1, a2 - 1)  # 0-based index

    # 3. Set coordinates
    conf = Chem.Conformer(len(lmpData.atom_records))
    for i, (atom_id, mol_tag, type_id, charge, x, y, z) in enumerate(lmpData.atom_records):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    # 4. Infer bond orders
    rdDetermineBonds.DetermineBondOrders(mol, charge=0)

    # 5. Sanitize
    Chem.SanitizeMol(mol)
    return mol