"""LAMMPS data file parsing and generation.

This module provides:
    - LammpsData: dataclass representing all sections of a LAMMPS data file
    - parse_lammps(path): Parse a LAMMPS data file and return LammpsData
    - generate_lammps(data, out_path): Write LammpsData to a LAMMPS file

The LAMMPS data file format consists of:
    - Header: atom/bond/angle/dihedral/improper counts and type counts
    - Box dimensions: xlo xhi, ylo yhi, zlo zhi
    - Masses: mass for each atom type
    - Pair Coeffs: epsilon and sigma for each atom type
    - Bond Coeffs: k and r0 for each bond type
    - Angle Coeffs: k and theta0 for each angle type
    - Dihedral Coeffs: OPLS coefficients (k0, k1, k2, k3) for each dihedral type
    - Improper Coeffs: coefficients for each improper type
    - Atoms: atom_id, molecule_tag, type_id, charge, x, y, z
    - Bonds: bond_id, bond_type, atom1, atom2
    - Angles: angle_id, angle_type, atom1, atom2, atom3
    - Dihedrals: dihedral_id, dihedral_type, atom1, atom2, atom3, atom4
    - Impropers: improper_id, improper_type, atom1, atom2, atom3, atom4
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from mappingff.db import MacroMapDB


@dataclass
class LammpsData:
    """Structured representation of a LAMMPS data file.

    This dataclass serves as the canonical internal format for LAMMPS data,
    bridging file I/O (parseLammps/generateLammps) and tool operations
    like parameterize.

    Attributes:
        header_comment: First line description.
        atoms: Number of atoms.
        bonds: Number of bonds.
        angles: Number of angles.
        dihedrals: Number of dihedrals.
        impropers: Number of impropers.
        atom_types: Number of atom types.
        bond_types: Number of bond types.
        angle_types: Number of angle types.
        dihedral_types: Number of dihedral types.
        improper_types: Number of improper types.
        xlo, xhi: X box dimensions.
        ylo, yhi: Y box dimensions.
        zlo, zhi: Z box dimensions.
        masses: List of (type_id, mass).
        pair_coeffs: List of (type_id, epsilon, sigma).
        bond_coeffs: List of (type_id, k, r0).
        angle_coeffs: List of (type_id, k, theta0).
        dihedral_coeffs: List of (type_id, k0, k1, k2, k3).
        improper_coeffs: List of (type_id, ...).
        atom_records: List of (atom_id, mol_tag, type_id, charge, x, y, z).
        bond_records: List of (bond_id, bond_type, a1, a2).
        angle_records: List of (angle_id, angle_type, a1, a2, a3).
        dihedral_records: List of (dih_id, dih_type, a1, a2, a3, a4).
        improper_records: List of (imp_id, imp_type, a1, a2, a3, a4).
    """

    header_comment: str = ""
    atoms: int = 0
    bonds: int = 0
    angles: int = 0
    dihedrals: int = 0
    impropers: int = 0
    atom_types: int = 0
    bond_types: int = 0
    angle_types: int = 0
    dihedral_types: int = 0
    improper_types: int = 0
    xlo: float = 0.0
    xhi: float = 0.0
    ylo: float = 0.0
    yhi: float = 0.0
    zlo: float = 0.0
    zhi: float = 0.0

    # Records as typed lists
    masses: list[tuple[int, float]] = field(default_factory=list)
    pair_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    bond_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    angle_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    dihedral_coeffs: list[tuple[int, float, float, float, float]] = field(
        default_factory=list
    )
    improper_coeffs: list[tuple[int, float, int, int]] = field(default_factory=list)

    atom_records: list[tuple[int, int, int, float, float, float, float]] = field(
        default_factory=list
    )
    bond_records: list[tuple[int, int, int, int]] = field(default_factory=list)
    angle_records: list[tuple[int, int, int, int, int]] = field(default_factory=list)
    dihedral_records: list[tuple[int, int, int, int, int, int]] = field(
        default_factory=list
    )
    improper_records: list[tuple[int, int, int, int, int, int]] = field(
        default_factory=list
    )


# LAMMPS section header patterns
_SECTION_HEADERS = [
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
]


def parse_lammps(path: Path) -> LammpsData:
    """Parse a LAMMPS data file.

    Args:
        path: Path to the LAMMPS data file (.lmp or .lammps).

    Returns:
        LammpsData object containing all parsed sections.
    """
    text = path.read_text()
    lines = text.splitlines()

    data = LammpsData()
    i = 0
    current_section: str | None = None

    # Parse header comment
    if lines:
        data.header_comment = lines[0].strip()
        i = 1

    # Parse counts section
    while i < len(lines):
        line = lines[i]
        if "#" in line:
            line = line[: line.index("#")].strip()
        else:
            line = line.strip()
        if not line:
            i += 1
            continue

        # Section header
        if line in _SECTION_HEADERS:
            current_section = line.lower().replace(" ", "_")
            i += 1
            continue

        if current_section is None:
            # Header section: counts, types, box dimensions
            parts = line.split()
            if len(parts) >= 2:
                try:
                    count = int(parts[0])
                except ValueError:
                    try:
                        if parts[2] == "xlo" and parts[3] == "xhi":
                            data.xlo, data.xhi = float(parts[0]), float(parts[1])
                        elif parts[2] == "ylo" and parts[3] == "yhi":
                            data.ylo, data.yhi = float(parts[0]), float(parts[1])
                        elif parts[2] == "zlo" and parts[3] == "zhi":
                            data.zlo, data.zhi = float(parts[0]), float(parts[1])
                    except (ValueError, IndexError):
                        pass
                    i += 1
                    continue

                key = parts[1]
                if key == "atoms":
                    data.atoms = count
                elif key == "bonds":
                    data.bonds = count
                elif key == "angles":
                    data.angles = count
                elif key == "dihedrals":
                    data.dihedrals = count
                elif key == "impropers":
                    data.impropers = count
                elif len(parts) == 3 and parts[2] == "types":
                    if parts[1] == "atom":
                        data.atom_types = count
                    elif parts[1] == "bond":
                        data.bond_types = count
                    elif parts[1] == "angle":
                        data.angle_types = count
                    elif parts[1] == "dihedral":
                        data.dihedral_types = count
                    elif parts[1] == "improper":
                        data.improper_types = count
        else:
            # Section data
            parts = line.split()
            if parts:
                if current_section == "masses":
                    data.masses.append((int(parts[0]), float(parts[1])))
                elif current_section == "pair_coeffs":
                    data.pair_coeffs.append(
                        (int(parts[0]), float(parts[1]), float(parts[2]))
                    )
                elif current_section == "bond_coeffs":
                    data.bond_coeffs.append(
                        (int(parts[0]), float(parts[1]), float(parts[2]))
                    )
                elif current_section == "angle_coeffs":
                    data.angle_coeffs.append(
                        (int(parts[0]), float(parts[1]), float(parts[2]))
                    )
                elif current_section == "dihedral_coeffs":
                    data.dihedral_coeffs.append(
                        (
                            int(parts[0]),
                            float(parts[1]),
                            float(parts[2]),
                            float(parts[3]),
                            float(parts[4]),
                        )
                    )
                elif current_section == "improper_coeffs":
                    data.improper_coeffs.append(
                        (int(parts[0]), float(parts[1]), int(parts[2]), int(parts[3]))
                    )
                elif current_section == "atoms":
                    data.atom_records.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            float(parts[3]),
                            float(parts[4]),
                            float(parts[5]),
                            float(parts[6]),
                        )
                    )
                elif current_section == "bonds":
                    data.bond_records.append(
                        (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]))
                    )
                elif current_section == "angles":
                    data.angle_records.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            int(parts[3]),
                            int(parts[4]),
                        )
                    )
                elif current_section == "dihedrals":
                    data.dihedral_records.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            int(parts[3]),
                            int(parts[4]),
                            int(parts[5]),
                        )
                    )
                elif current_section == "impropers":
                    data.improper_records.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            int(parts[3]),
                            int(parts[4]),
                            int(parts[5]),
                        )
                    )

        i += 1

    return data


def generate_lammps(data: LammpsData, out_path: Path) -> None:
    """Generate a LAMMPS data file from LammpsData.

    Args:
        data: LammpsData object containing all sections.
        out_path: Path to write the output LAMMPS file.
    """
    lines: list[str] = []

    # Header comment
    lines.append(data.header_comment or "LAMMPS data file")
    lines.append("")

    # Counts
    lines.append(f"{data.atoms:>10} atoms")
    lines.append(f"{data.bonds:>10} bonds")
    lines.append(f"{data.angles:>10} angles")
    lines.append(f"{data.dihedrals:>10} dihedrals")
    lines.append(f"{data.impropers:>10} impropers")
    lines.append("")
    lines.append(f"{data.atom_types:>10} atom types")
    lines.append(f"{data.bond_types:>10} bond types")
    lines.append(f"{data.angle_types:>10} angle types")
    lines.append(f"{data.dihedral_types:>10} dihedral types")
    lines.append(f"{data.improper_types:>10} improper types")
    lines.append("")

    # Box dimensions (6 decimal places)
    lines.append(f"{data.xlo:>12.6f} {data.xhi:>12.6f} xlo xhi")
    lines.append(f"{data.ylo:>12.6f} {data.yhi:>12.6f} ylo yhi")
    lines.append(f"{data.zlo:>12.6f} {data.zhi:>12.6f} zlo zhi")
    lines.append("")

    # Masses (3 decimal places)
    lines.append("Masses")
    lines.append("")
    for type_id, mass in data.masses:
        lines.append(f"{type_id:>10} {mass:>10.3f}")
    lines.append("")

    # Pair Coeffs (epsilon: 3 decimal places, sigma: 7 decimal places)
    lines.append("Pair Coeffs")
    lines.append("")
    for type_id, epsilon, sigma in data.pair_coeffs:
        lines.append(f"{type_id:>10} {epsilon:>14.3f} {sigma:>14.7f}")
    lines.append("")

    # Bond Coeffs (K0, R0: 4 decimal places)
    lines.append("Bond Coeffs")
    lines.append("")
    for type_id, k, r0 in data.bond_coeffs:
        lines.append(f"{type_id:>10} {k:>14.4f} {r0:>14.4f}")
    lines.append("")

    # Angle Coeffs (K0, angle0: 3 decimal places)
    lines.append("Angle Coeffs")
    lines.append("")
    for type_id, k, theta0 in data.angle_coeffs:
        lines.append(f"{type_id:>10} {k:>14.3f} {theta0:>14.3f}")
    lines.append("")

    # Dihedral Coeffs (type_id: integer, V1-V4: 3 decimal places)
    lines.append("Dihedral Coeffs")
    lines.append("")
    for rec in data.dihedral_coeffs:
        type_id, k0, k1, k2, k3 = rec
        line = f"{type_id:>10}" + "".join(f"{x:>14.3f}" for x in (k0, k1, k2, k3))
        lines.append(line)
    lines.append("")

    # Improper Coeffs (V2/2: 3 decimal places, last two fields are integers)
    lines.append("Improper Coeffs")
    lines.append("")
    for rec in data.improper_coeffs:  # type: ignore[assignment]
        # rec[0] is type_id, rec[1:] are coeffs
        # last two fields are integers
        formatted = [f"{rec[0]:>10}"]
        for x in rec[1:-2]:  # type: ignore[misc]
            formatted.append(f"{x:>14.3f}")
        formatted.append(f"{int(rec[-2]):>14}")
        formatted.append(f"{int(rec[-1]):>14}")
        lines.append("".join(formatted))
    lines.append("")

    # Atoms (charge: 6 decimal places, x/y/z: 5 decimal places)
    lines.append("Atoms")
    lines.append("")
    for atom_id, mol_tag, type_id, charge, x, y, z in data.atom_records:
        lines.append(
            f"{atom_id:>10} {mol_tag:>5} {type_id:>5} {charge:>14.6f} "
            f"{x:>14.5f} {y:>14.5f} {z:>14.5f}"
        )
    lines.append("")

    # Bonds
    lines.append("Bonds")
    lines.append("")
    for bond_id, bond_type, a1, a2 in data.bond_records:
        lines.append(f"{bond_id:>10} {bond_type:>10} {a1:>10} {a2:>10}")
    lines.append("")

    # Angles
    lines.append("Angles")
    lines.append("")
    for angle_id, angle_type, a1, a2, a3 in data.angle_records:
        lines.append(f"{angle_id:>10} {angle_type:>10} {a1:>10} {a2:>10} {a3:>10}")
    lines.append("")

    # Dihedrals
    lines.append("Dihedrals")
    lines.append("")
    for dih_id, dih_type, a1, a2, a3, a4 in data.dihedral_records:
        lines.append(f"{dih_id:>10} {dih_type:>10} {a1:>10} {a2:>10} {a3:>10} {a4:>10}")
    lines.append("")

    # Impropers
    lines.append("Impropers")
    lines.append("")
    for imp_id, imp_type, a1, a2, a3, a4 in data.improper_records:
        lines.append(f"{imp_id:>10} {imp_type:>10} {a1:>10} {a2:>10} {a3:>10} {a4:>10}")
    lines.append("")

    out_path.write_text("\n".join(lines) + "\n")


def _solve_weighted_adjustment(
    delta: float,
    indices: list,
    charges: list,
    min_bounds: list,
    max_bounds: list,
    weights: list,
) -> list:
    """Solve weighted charge adjustment with per-atom bounds via iterative fitting.

    Each atom i has current charge q_i, bounds [min_i, max_i], and weight w_i.
    We want to find adjustment adj_i such that:
    - sum(adj_i) = delta
    - min_i - q_i <= adj_i <= max_i - q_i

    Uses iterative approach: distribute delta proportionally by weight,
    but any atom that hits a bound stops participating and its share
    is redistributed to remaining atoms.

    Args:
        delta: Total charge to distribute.
        indices: Data indices for atoms.
        charges: Current charges.
        min_bounds: Per-atom minimum bounds.
        max_bounds: Per-atom maximum bounds.
        weights: Per-atom weights (|charge|).

    Returns:
        List of adjustments per atom.
    """
    adjustments = [0.0] * len(indices)
    remaining_delta = delta
    active_mask = [True] * len(indices)

    while abs(remaining_delta) > 1e-9 and any(active_mask):
        # Find active atoms
        active_indices = [i for i, a in enumerate(active_mask) if a]
        if not active_indices:
            break

        # Calculate total weight of active atoms
        total_weight = sum(weights[i] for i in active_indices)
        if total_weight < 1e-9:
            break

        # Distribute remaining delta proportionally
        overshoot = False
        for i in active_indices:
            w = weights[i]
            proportion = w / total_weight
            ideal_adj = remaining_delta * proportion

            # Check if this adjustment would exceed bounds
            max_possible = max_bounds[i] - charges[i]
            min_possible = min_bounds[i] - charges[i]

            if ideal_adj > max_possible:
                adjustments[i] = max_possible
                active_mask[i] = False
                overshoot = True
            elif ideal_adj < min_possible:
                adjustments[i] = min_possible
                active_mask[i] = False
                overshoot = True
            else:
                adjustments[i] = ideal_adj

        if overshoot:
            # Recalculate remaining delta after hitting bounds
            applied = sum(adjustments)
            remaining_delta = delta - applied
        else:
            remaining_delta = 0.0
            break

    return adjustments


def adjust_total_charge(
    data: LammpsData,
    target_charge: float | None,
    db: MacroMapDB,
    atom_type_params: dict,
) -> tuple[float, float]:
    """Adjust atom charges in LammpsData to achieve target total charge.

    Distributes charge delta using a weighted approach:
    1. First to atoms with charge_list entries > 1 (proportional to |charge|)
    2. Then to sp3 carbons meeting criteria if residual remains

    Args:
        data: LammpsData object to modify in-place.
        target_charge: Desired total charge for the system. If None, no adjustment.
        db: MacroMapDB with atomTypes information.
        atom_type_params: Dict mapping db_lammps_type -> {type, element, mass, ...}.

    Returns:
        Tuple of (charge_after_step1, charge_after_step2).
    """
    current_charge = sum(atom[3] for atom in data.atom_records)

    if target_charge is None:
        return current_charge, current_charge

    delta = target_charge - current_charge
    if abs(delta) < 1e-9:
        return current_charge, current_charge

    # Build output_type -> db_lammps_type mapping
    output_to_db_type = {
        params["type"]: db_type for db_type, params in atom_type_params.items()
    }

    # Build (element, db_lammps_type) -> charge_list lookup
    element_lammps_to_chargelist = {}
    for entry in db.atom_types.values():
        charge_list = entry.get("charge_list", [])
        if len(charge_list) > 1:
            key = (entry["element"], entry["lammps_type"])
            element_lammps_to_chargelist[key] = charge_list

    # Collect multi-entry atoms for step 1
    multi_entry = []
    for data_idx, atom_rec in enumerate(data.atom_records):
        output_type = atom_rec[2]
        if output_type not in output_to_db_type:
            continue
        db_type = output_to_db_type[output_type]
        if db_type is None:
            continue
        element = atom_type_params[db_type]["element"]
        key = (element, db_type)
        if key not in element_lammps_to_chargelist:
            continue
        charge_list = element_lammps_to_chargelist[key]
        current_q = atom_rec[3]
        if abs(current_q) < 1e-9:
            continue
        multi_entry.append(
            (data_idx, current_q, min(charge_list), max(charge_list), abs(current_q))
        )

    charge_after_step1 = current_charge
    if multi_entry:
        indices = [d[0] for d in multi_entry]
        charges = [d[1] for d in multi_entry]
        min_bounds = [d[2] for d in multi_entry]
        max_bounds = [d[3] for d in multi_entry]
        weights = [d[4] for d in multi_entry]

        adjustments = _solve_weighted_adjustment(
            delta, indices, charges, min_bounds, max_bounds, weights
        )

        for i, data_idx in enumerate(indices):
            atom = list(data.atom_records[data_idx])
            atom[3] = charges[i] + adjustments[i]
            data.atom_records[data_idx] = tuple(atom)  # type: ignore[assignment]

        charge_after_step1 = sum(atom[3] for atom in data.atom_records)
        delta = target_charge - charge_after_step1

    # Step 2: evenly distribute residual across all atoms
    charge_after_step2 = charge_after_step1
    if abs(delta) > 1e-9:
        adj_per_atom = delta / len(data.atom_records)
        for i in range(len(data.atom_records)):
            atom = list(data.atom_records[i])
            atom[3] += adj_per_atom
            data.atom_records[i] = tuple(atom)  # type: ignore[assignment]
        charge_after_step2 = sum(atom[3] for atom in data.atom_records)

    return charge_after_step1, charge_after_step2
