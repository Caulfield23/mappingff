#!/usr/bin/env python3
from pathlib import Path


SECTION_NAMES = {
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
    "Velocities",
}


def _parse_atom_line(toks):
    if len(toks) < 6:
        return None

    atom_id = int(toks[0])

    # Support both atom-style formats:
    #  id mol type q x y z  and  id type q x y z
    candidates = []
    if len(toks) >= 7:
        candidates.append((2, 3, 4, 5, 6))
    if len(toks) >= 6:
        candidates.append((1, 2, 3, 4, 5))

    for t_col, q_col, x_col, y_col, z_col in candidates:
        try:
            atom_type = int(toks[t_col])
            charge = float(toks[q_col])
            x = float(toks[x_col])
            y = float(toks[y_col])
            z = float(toks[z_col])
            return atom_id, atom_type, charge, x, y, z
        except Exception:
            continue

    return None


def parse_lammps_sections(lmp_path: Path):
    lines = lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    current = None
    masses = {}
    pair_coeffs = {}
    atoms = {}

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped in SECTION_NAMES:
            current = stripped
            continue

        if stripped.startswith("#"):
            continue

        if current == "Masses":
            toks = stripped.split()
            if len(toks) < 2:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            type_id = int(toks[0])
            masses[type_id] = float(toks[1])

        elif current == "Pair Coeffs":
            toks = stripped.split()
            if len(toks) < 3:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            type_id = int(toks[0])
            pair_coeffs[type_id] = {
                "epsilon": float(toks[1]),
                "sigma": float(toks[2]),
            }

        elif current == "Atoms":
            toks = stripped.split()
            if not toks[0].lstrip("+-").isdigit():
                continue
            parsed = _parse_atom_line(toks)
            if parsed is None:
                continue
            atom_id, atom_type, charge, x, y, z = parsed
            atoms[atom_id] = {
                "lmp_type": atom_type,
                "charge": charge,
                "x": x,
                "y": y,
                "z": z,
            }

    return masses, pair_coeffs, atoms


def parse_lammps_masses(lmp_path: Path):
    masses, _, _ = parse_lammps_sections(lmp_path)
    if not masses:
        raise ValueError(f"No Masses section parsed from: {lmp_path}")
    return masses
