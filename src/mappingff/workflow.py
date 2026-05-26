"""Core workflow functions for mappingff.

This module provides the main workflow functions:
    - build_db: Build parameter database from sample molecules
    - parameterize: Parameterize target molecule and generate LAMMPS file

These functions can be imported and used programmatically, or via the CLI.
"""

from __future__ import annotations

import logging
from pathlib import Path

from mappingff.db import MacroMapDB
from mappingff.encode import get_hop0_subgraph, get_hop3_subgraph
from mappingff.fallback import resolve_atom_type
from mappingff.lmp import LammpsData, adjust_total_charge, generate_lammps, parse_lammps
from mappingff.mol import MolReader, compute_hop_keys


def build_db(samples_dir: Path, db_path: Path, append: bool = False) -> None:
    """Build a parameter database from sample molecules.

    Iterates over subdirectories in samples_dir, each expected to contain
    a .mol (or .pdb) and a lammps data file with the same base name.

    Args:
        samples_dir: Directory containing sample subdirectories.
        db_path: Path to output the SQLite database file.
        append: If False (default), replace existing database.
               If True, merge with existing data.
    """
    log = logging.getLogger("build-db")

    # Find all sample directories
    sample_dirs = sorted([d for d in samples_dir.iterdir() if d.is_dir()])

    # Initialize database
    db = MacroMapDB(db_path)
    db.load(append=append)

    total_atoms = 0
    total_bonds = 0
    total_angles = 0
    total_dihedrals = 0
    total_impropers = 0

    # Collect all items to process: (lmp_path, topo_path or None)
    items: list[tuple[Path, Path | None]] = []
    for sample_dir in sample_dirs:
        topo_path = next(sample_dir.glob("*.mol"), None) or next(
            sample_dir.glob("*.pdb"), None
        )
        lmp_path = (
            next(sample_dir.glob("*.lmp"), None)
            or next(sample_dir.glob("*.lammps"), None)
            or next(sample_dir.glob("*.data"), None)
        )
        if topo_path and lmp_path:
            items.append((lmp_path, topo_path))
    for pattern in ("*.lmp", "*.lammps", "*.data"):
        for lmp_path in samples_dir.glob(pattern):
            items.append((lmp_path, None))

    for lmp_path, topo_path in items:
        seg_name = lmp_path.stem

        file_part = f" + {topo_path.name}" if topo_path else " (lmp only)"
        log.info(f"Processing {seg_name}: {lmp_path.name}{file_part}")

        lmp_data = parse_lammps(lmp_path)
        molReader = MolReader(topo_path if topo_path else lmp_path)

        atoms = molReader.get_atoms()
        bonds = molReader.get_bonds()
        angles = molReader.get_angles()
        dihedrals = molReader.get_dihedrals()
        impropers = molReader.get_impropers()
        total_atoms += len(atoms)
        total_bonds += len(bonds)
        total_angles += len(angles)
        total_dihedrals += len(dihedrals)
        total_impropers += len(impropers)

        # Build atom_idx -> hop0_key mapping for bonded parameter lookup
        atom_hop0_key: dict[int, str] = {}

        # Process each atom
        for atom in atoms:
            atom_idx = atom["idx"]
            hop3_key, hop2_key, hop1_key, hop0_key = compute_hop_keys(
                molReader, atom_idx
            )
            rd_atom = molReader.mol.GetAtomWithIdx(atom_idx - 1)
            hop0_graph = get_hop0_subgraph(molReader.mol, rd_atom)
            hop3_graph = get_hop3_subgraph(molReader.mol, rd_atom)
            atom_hop0_key[atom_idx] = hop0_key

            # Get LAMMPS type from the sample data
            atom_record = next(
                (rec for rec in lmp_data.atom_records if rec[0] == atom_idx),
                None,
            )
            if atom_record is None:
                log.warning(f"Atom {atom_idx} not found in LAMMPS data for {seg_name}")
                continue

            # atom_record format: [atom_id, mol_tag, type_id, charge, x, y, z]
            lammps_type = atom_record[2]

            # Get pair coeffs for this type
            pair_coeff = next(
                (rec for rec in lmp_data.pair_coeffs if rec[0] == lammps_type),
                None,
            )
            mass = next(
                (rec for rec in lmp_data.masses if rec[0] == lammps_type),
                None,
            )

            # Insert atom type at hop3 level (finest classification)
            # lammps_type is auto-generated to ensure uniqueness per hop3_key
            db.insert_atom_type(
                {
                    "hop3_key": hop3_key,
                    "hop2_key": hop2_key,
                    "hop1_key": hop1_key,
                    "hop0_key": hop0_key,
                    "element": atom["symbol"],
                    "hop0_graph": hop0_graph,
                    "hop3_graph": hop3_graph,
                    "mass": mass[1] if mass else 0.0,
                    "sigma": pair_coeff[2] if pair_coeff else 0.0,
                    "epsilon": pair_coeff[1] if pair_coeff else 0.0,
                    "charge": atom_record[3],
                    "source": f"{seg_name}_{atom_idx}",
                },
            )

            # Get the auto-generated lammps_type from database
            inserted_info = db.get_atom_type(hop3_key)
            if inserted_info is None:
                log.warning(
                    f"Atom type for hop3_key {hop3_key[:8]}... not found after insert"
                )
                continue
            assigned_lammps_type = inserted_info["lammps_type"]

            # Insert hop keymap (for external validation only)
            db.insert_hop_key([hop0_key, hop1_key, hop2_key], assigned_lammps_type)

        # Process bonds
        # bond_records: [bond_id, bond_type, a1, a2]
        for bond_rec in lmp_data.bond_records:
            bond_type = bond_rec[1]
            a1 = bond_rec[2]
            a2 = bond_rec[3]

            hop0_key_a = atom_hop0_key[a1]
            hop0_key_b = atom_hop0_key[a2]

            # bond_coeffs: [type_id, K, r0] (bond_style harmonic)
            bond_coeff = next(
                rec for rec in lmp_data.bond_coeffs if rec[0] == bond_type
            )

            # Canonical key (order-independent)
            bond_key = (min(hop0_key_a, hop0_key_b), max(hop0_key_a, hop0_key_b))

            db.insert_bond_param(bond_key, {"k": bond_coeff[1], "r0": bond_coeff[2]})

        # Process angles
        # angle_records: [angle_id, angle_type, a1, a2, a3]
        for angle_rec in lmp_data.angle_records:
            ang_type = angle_rec[1]
            a1 = angle_rec[2]
            a2 = angle_rec[3]
            a3 = angle_rec[4]

            hop0_key_a = atom_hop0_key[a1]
            hop0_key_b = atom_hop0_key[a2]
            hop0_key_c = atom_hop0_key[a3]

            # angle_coeffs: [type_id, K, theta0] (angle_style harmonic)
            angle_coeff = next(
                rec for rec in lmp_data.angle_coeffs if rec[0] == ang_type
            )

            # Canonical key (outer atoms swappable)
            angle_key: tuple[str, str, str] = min(
                (hop0_key_a, hop0_key_b, hop0_key_c),
                (hop0_key_c, hop0_key_b, hop0_key_a),
            )

            db.insert_angle_param(
                angle_key, {"k": angle_coeff[1], "theta0": angle_coeff[2]}
            )

        # Process dihedrals
        # dihedral_records: [dih_id, dih_type, a1, a2, a3, a4]
        for dih_rec in lmp_data.dihedral_records:
            dih_type = dih_rec[1]
            a1 = dih_rec[2]
            a2 = dih_rec[3]
            a3 = dih_rec[4]
            a4 = dih_rec[5]

            hop0_key_a = atom_hop0_key[a1]
            hop0_key_b = atom_hop0_key[a2]
            hop0_key_c = atom_hop0_key[a3]
            hop0_key_d = atom_hop0_key[a4]

            # dih_coeffs: [type_id, K1, K2, K3, K4] (dihedral_style opls)
            dih_coeff = next(
                rec for rec in lmp_data.dihedral_coeffs if rec[0] == dih_type
            )

            # Canonical key (A,B,C,D) and (D,C,B,A) are equivalent, pick lexicographically smaller
            dihedral_key: tuple[str, str, str, str] = min(
                (hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d),
                (hop0_key_d, hop0_key_c, hop0_key_b, hop0_key_a),
            )

            db.insert_dihedral_param(dihedral_key, {"coeffs": dih_coeff[1:]})

        # Process impropers
        # improper_records: [imp_id, imp_type, a1, a2, a3, a4]
        for imp_rec in lmp_data.improper_records:
            imp_type = imp_rec[1]
            a1 = imp_rec[2]
            a2 = imp_rec[3]
            a3 = imp_rec[4]
            a4 = imp_rec[5]

            hop0_key_a = atom_hop0_key[a1]
            hop0_key_b = atom_hop0_key[a2]
            hop0_key_c = atom_hop0_key[a3]
            hop0_key_d = atom_hop0_key[a4]

            # imp_coeffs: [type_id, K, d, n] (improper_style cvff)
            imp_coeff = next(
                rec for rec in lmp_data.improper_coeffs if rec[0] == imp_type
            )

            # Canonical key (center fixed, others sorted)
            others = sorted([hop0_key_b, hop0_key_c, hop0_key_d])
            improper_key: tuple[str, str, str, str] = (
                hop0_key_a,
                others[0],
                others[1],
                others[2],
            )

            db.insert_improper_param(improper_key, {"coeffs": imp_coeff[1:]})

        log.info(
            f"  {seg_name}: {len(atoms)} atoms, {len(bonds)} bonds, "
            f"{len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers"
        )

        # Update metadata for this sample
        db.set_meta(
            seg_name,
            f"{lmp_path} | {topo_path if topo_path else '(lmp only)'}",
        )

    # Update metadata
    db.set_meta("sample_count", str(len(items)))

    # Save database
    db.save()

    log.info(f"Database saved to {db_path}")
    log.info(f"Samples: {len(items)}")
    log.info(
        f"Processed: {total_atoms} atoms, {total_bonds} bonds, "
        f"{total_angles} angles, {total_dihedrals} dihedrals, {total_impropers} impropers"
    )
    log.info(f"Hop3 atom types: {len(db.atom_types)}")
    log.info(f"Hop2 atom types: {len(db.hop2_keymap)}")
    log.info(f"Hop1 atom types: {len(db.hop1_keymap)}")
    log.info(f"Hop0 atom types: {len(db.hop0_keymap)}")
    log.info(f"Bond types: {len(db.bond_params)}")
    log.info(f"Angle types: {len(db.angle_params)}")
    log.info(f"Dihedral types: {len(db.dihedral_params)}")
    log.info(f"Improper types: {len(db.improper_params)}")


def parameterize(
    topo_path: Path,
    db_path: Path,
    out_path: Path | None = None,
    total_charge: float | None = None,
) -> dict:
    """Parameterize a target molecule using the database.

    Args:
        topo_path: Path to target molecule .mol or .pdb file.
        db_path: Path to the SQLite database file.
        out_path: Output LAMMPS file path. If None, uses <mol_file>.data.
        total_charge: Target total charge for the system. If None, no adjustment (default: None).
            If provided, residual charge will be adjusted evenly across all atoms.

    Returns:
        Dictionary with parameterization statistics.
    """
    log = logging.getLogger("parameterize")

    # Load database
    db = MacroMapDB(db_path)
    db.load(append=True)
    log.info(f"Database loaded from {db_path}")
    log.debug(f"  Hop3 atom types: {len(db.atom_types)}")
    log.debug(f"  Hop2 atom types: {len(db.hop2_keymap)}")
    log.debug(f"  Hop1 atom types: {len(db.hop1_keymap)}")
    log.debug(f"  Hop0 atom types: {len(db.hop0_keymap)}")
    log.debug(f"  Bond types: {len(db.bond_params)}")
    log.debug(f"  Angle types: {len(db.angle_params)}")
    log.debug(f"  Dihedral types: {len(db.dihedral_params)}")
    log.debug(f"  Improper types: {len(db.improper_params)}")

    # Parse target molecule
    molReader = MolReader(topo_path)
    atoms = molReader.get_atoms()
    bonds = molReader.get_bonds()
    angles = molReader.get_angles()
    dihedrals = molReader.get_dihedrals()
    impropers = molReader.get_impropers()
    coords = molReader.get_coords()

    log.info(f"Parameterizing molecule: {topo_path.name}")

    # Determine output path
    if out_path is None:
        out_path = topo_path.with_suffix(".data")

    # Resolve atom types
    atom_type_map: dict[int, int] = {}  # atom_idx -> lammps_type (resolved from db)
    atom_type_params: dict[int | None, dict] = (
        {}
    )  # lammps_type -> {type: consecutive_output_type, element, mass, sigma, epsilon, charge}
    atom_hop0_key: dict[int, str | None] = {}  # atom_idx -> hop0_key (resolved from db)
    atom_hop3_key: dict[int, str] = {}  # atom_idx -> hop3_key (for debugging)
    atom_fallback_level: dict[int, str] = (
        {}
    )  # atom_idx -> fallback level (hop2/hop1/hop0)
    next_atom_type_id = 1
    hop3_matches = 0
    hop2_matches = 0
    hop1_matches = 0
    hop0_matches = 0
    no_match = 0

    for atom in atoms:
        atom_idx = atom["idx"]
        hop3_key, hop2_key, hop1_key, hop0_key = compute_hop_keys(molReader, atom_idx)
        atom_hop0_key[atom_idx] = hop0_key
        atom_hop3_key[atom_idx] = hop3_key

        lammps_type, resolved_hop0_key, hop_level, matched_hop3_key = resolve_atom_type(
            hop3_key, hop2_key, hop1_key, hop0_key, db
        )
        atom_type_map[atom_idx] = lammps_type  # type: ignore[assignment]
        atom_hop0_key[atom_idx] = resolved_hop0_key

        if lammps_type is not None and lammps_type not in atom_type_params:
            info = db.atom_types[matched_hop3_key]  # type: ignore[index]
            atom_type_params[lammps_type] = {
                "type": next_atom_type_id,
                "element": info["element"],
                "mass": info["mass"],
                "sigma": info["sigma"],
                "epsilon": info["epsilon"],
                "charge": info["charge"],
            }
            next_atom_type_id += 1
        elif lammps_type is None:
            if None not in atom_type_params:
                atom_type_params[None] = {
                    "type": next_atom_type_id,
                    "element": "unknown",
                    "mass": 0.0,
                    "sigma": 0.0,
                    "epsilon": 0.0,
                    "charge": 0.0,
                }
                next_atom_type_id += 1

        # Track match statistics
        if hop_level == "hop3":
            hop3_matches += 1
            log.debug(
                f"  Atom {atom_idx}: element={atom['symbol']}, lammps_type={lammps_type}"
            )
        elif hop_level == "hop2":
            hop2_matches += 1
            atom_fallback_level[atom_idx] = "hop2"
            log.debug(
                f"  Atom {atom_idx}: element={atom['symbol']}, lammps_type={lammps_type} (hop2 fallback)"
            )
        elif hop_level == "hop1":
            hop1_matches += 1
            atom_fallback_level[atom_idx] = "hop1"
            log.debug(
                f"  Atom {atom_idx}: element={atom['symbol']}, lammps_type={lammps_type} (hop1 fallback)"
            )
        elif hop_level == "hop0":
            hop0_matches += 1
            atom_fallback_level[atom_idx] = "hop0"
            log.debug(
                f"  Atom {atom_idx}: element={atom['symbol']}, lammps_type={lammps_type} (hop0 fallback)"
            )
        else:
            no_match += 1
            atom_fallback_level[atom_idx] = "none"
            log.debug(
                f"  Atom {atom_idx}: element={atom['symbol']}, lammps_type={lammps_type} (NO MATCH)"
            )

    # Look up bond parameters
    bond_type_map: dict[int, int] = {}  # bond_idx -> new_bond_type_id
    bond_type_params: dict[int, tuple] = {}  # new_bond_type_id -> (k, r0)
    next_bond_type_id = 1
    bond_no_match = 0
    bond_no_match_atoms: list[tuple[int, int]] = []

    for bond in bonds:
        bond_idx = bond["idx"]
        a1 = bond["a1"]
        a2 = bond["a2"]
        hop0_key_a = atom_hop0_key.get(a1)
        hop0_key_b = atom_hop0_key.get(a2)

        # Look up bond param - zeroed if hop0_key is None or lookup returns None
        if hop0_key_a is None or hop0_key_b is None:
            bond_no_match += 1
            bond_no_match_atoms.append((a1, a2))
            param_key = (0.0, 0.0)
        else:
            bond_param = db.lookup_bond_param(hop0_key_a, hop0_key_b)
            if bond_param is None:
                bond_no_match += 1
                bond_no_match_atoms.append((a1, a2))
                param_key = (0.0, 0.0)
            else:
                param_key = (bond_param["k"], bond_param["r0"])

        # Create bond type if not seen
        bond_type_id = None
        for btid, (bk, br0) in bond_type_params.items():
            if abs(bk - param_key[0]) < 0.01 and abs(br0 - param_key[1]) < 0.001:
                bond_type_id = btid
                break

        if bond_type_id is None:
            bond_type_id = next_bond_type_id
            next_bond_type_id += 1
            bond_type_params[bond_type_id] = param_key

        bond_type_map[bond_idx] = bond_type_id

    # Look up angle parameters
    angle_type_map: dict[int, int] = {}  # angle_idx -> new_angle_type_id
    angle_type_params: dict[int, tuple] = {}  # new_angle_type_id -> (k, theta0)
    next_angle_type_id = 1
    angle_no_match = 0
    angle_no_match_atoms: list[tuple[int, int, int]] = []

    for angle in angles:
        angle_idx = angle["idx"]
        a1 = angle["a1"]
        a2 = angle["a2"]
        a3 = angle["a3"]
        hop0_key_a = atom_hop0_key.get(a1)
        hop0_key_b = atom_hop0_key.get(a2)
        hop0_key_c = atom_hop0_key.get(a3)

        if hop0_key_a is None or hop0_key_b is None or hop0_key_c is None:
            angle_no_match += 1
            angle_no_match_atoms.append((a1, a2, a3))
            param_key = (0.0, 0.0)
        else:
            angle_param = db.lookup_angle_param(hop0_key_a, hop0_key_b, hop0_key_c)
            if angle_param is None:
                angle_no_match += 1
                angle_no_match_atoms.append((a1, a2, a3))
                param_key = (0.0, 0.0)
            else:
                param_key = (angle_param["k"], angle_param["theta0"])

        # Create angle type if not seen
        angle_type_id = None
        for atid, (ak, ath) in angle_type_params.items():
            if abs(ak - param_key[0]) < 0.01 and abs(ath - param_key[1]) < 0.1:
                angle_type_id = atid
                break

        if angle_type_id is None:
            angle_type_id = next_angle_type_id
            next_angle_type_id += 1
            angle_type_params[angle_type_id] = param_key

        angle_type_map[angle_idx] = angle_type_id

    # Look up dihedral parameters
    dihedral_type_map: dict[int, int] = {}  # dih_idx -> new_dih_type_id
    dihedral_type_params: dict[int, tuple[float, float, float, float]] = (
        {}
    )  # new_dih_type_id -> coeffs tuple
    next_dihedral_type_id = 1
    dihedral_no_match = 0
    dihedral_no_match_atoms: list[tuple[int, int, int, int]] = []

    for dihedral in dihedrals:
        dih_idx = dihedral["idx"]
        a1 = dihedral["a1"]
        a2 = dihedral["a2"]
        a3 = dihedral["a3"]
        a4 = dihedral["a4"]
        hop0_key_a = atom_hop0_key.get(a1)
        hop0_key_b = atom_hop0_key.get(a2)
        hop0_key_c = atom_hop0_key.get(a3)
        hop0_key_d = atom_hop0_key.get(a4)

        if (
            hop0_key_a is None
            or hop0_key_b is None
            or hop0_key_c is None
            or hop0_key_d is None
        ):
            dihedral_no_match += 1
            dihedral_no_match_atoms.append((a1, a2, a3, a4))
            coeffs = (0.0, 0.0, 0.0, 0.0)
        else:
            dihedral_param = db.lookup_dihedral_param(
                hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d
            )
            if dihedral_param is None:
                dihedral_no_match += 1
                dihedral_no_match_atoms.append((a1, a2, a3, a4))
                coeffs = (0.0, 0.0, 0.0, 0.0)
            else:
                coeffs = tuple(dihedral_param.get("coeffs", []))  # type: ignore[assignment]

        # Create dihedral type if not seen
        dihedral_type_id = None
        for dtid, existing_coeffs in dihedral_type_params.items():
            if len(coeffs) == len(existing_coeffs):
                match = all(abs(c - e) < 0.01 for c, e in zip(coeffs, existing_coeffs))
                if match:
                    dihedral_type_id = dtid
                    break

        if dihedral_type_id is None:
            dihedral_type_id = next_dihedral_type_id
            next_dihedral_type_id += 1
            dihedral_type_params[dihedral_type_id] = coeffs

        dihedral_type_map[dih_idx] = dihedral_type_id

    # Look up improper parameters
    improper_type_map: dict[int, int] = {}  # imp_idx -> new_imp_type_id
    improper_type_params: dict[int, tuple[float, int, int]] = (
        {}
    )  # new_imp_type_id -> coeffs tuple
    next_improper_type_id = 1
    improper_no_match = 0
    improper_no_match_atoms: list[tuple[int, int, int, int]] = []

    # Build atom and atom degree index for cvff improper filtering
    atom_name_by_idx = {a["idx"]: a for a in atoms}
    atom_degree_by_idx = {a["idx"]: a["degree"] for a in atoms}

    for improper in impropers:
        imp_idx = improper["idx"]
        a1 = improper["a1"]
        a2 = improper["a2"]
        a3 = improper["a3"]
        a4 = improper["a4"]

        # improper filter: center atom (a1) must be C or N with exactly 3 neighbors
        center_atom = atom_name_by_idx.get(a1)
        if center_atom is None:
            continue
        if center_atom["symbol"] not in ("C", "N"):
            continue
        if atom_degree_by_idx.get(a1, 0) != 3:
            continue

        hop0_key_a = atom_hop0_key.get(a1)
        hop0_key_b = atom_hop0_key.get(a2)
        hop0_key_c = atom_hop0_key.get(a3)
        hop0_key_d = atom_hop0_key.get(a4)

        if (
            hop0_key_a is None
            or hop0_key_b is None
            or hop0_key_c is None
            or hop0_key_d is None
        ):
            improper_no_match += 1
            improper_no_match_atoms.append((a1, a2, a3, a4))
            coeffs = (0.0, -1, 2)  # type: ignore[assignment]
        else:
            improper_param = db.lookup_improper_param(
                hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d
            )
            if improper_param is None:
                improper_no_match += 1
                improper_no_match_atoms.append((a1, a2, a3, a4))
                coeffs = (0.0, -1, 2)  # type: ignore[assignment]
            else:
                coeffs = tuple(improper_param.get("coeffs", []))  # type: ignore[assignment]

        # Create improper type if not seen
        improper_type_id = None
        for itid, existing_coeffs in improper_type_params.items():  # type: ignore[assignment]
            if len(coeffs) == len(existing_coeffs):
                match = all(abs(c - e) < 0.01 for c, e in zip(coeffs, existing_coeffs))
                if match:
                    improper_type_id = itid
                    break

        if improper_type_id is None:
            improper_type_id = next_improper_type_id
            next_improper_type_id += 1
            improper_type_params[improper_type_id] = coeffs  # type: ignore[assignment]

        improper_type_map[imp_idx] = improper_type_id

    # Write all warnings first
    if bond_no_match > 0:
        for a1, a2 in bond_no_match_atoms:
            log.warning(f"  No bond param for atoms {a1}-{a2}")
    if angle_no_match > 0:
        for a1, a2, a3 in angle_no_match_atoms:
            log.warning(f"  No angle param for atoms {a1}-{a2}-{a3}")
    if dihedral_no_match > 0:
        for a1, a2, a3, a4 in dihedral_no_match_atoms:
            log.warning(f"  No dihedral param for atoms {a1}-{a2}-{a3}-{a4}")
    if improper_no_match > 0:
        for a1, a2, a3, a4 in improper_no_match_atoms:
            log.warning(f"  No improper param for atoms {a1}-{a2}-{a3}-{a4}")

    # Write all parameter info lines
    log.info(
        f"Bond parameters: {len(bond_type_map)} bonds, {len(bond_type_params)} unique types, {bond_no_match} no match"
    )
    log.info(
        f"Angle parameters: {len(angle_type_map)} angles, {len(angle_type_params)} unique types, {angle_no_match} no match"
    )
    log.info(
        f"Dihedral parameters: {len(dihedral_type_map)} dihedrals, {len(dihedral_type_params)} unique types, {dihedral_no_match} no match"
    )
    log.info(
        f"Improper parameters: {len(improper_type_map)} impropers, {len(improper_type_params)} unique types, {improper_no_match} no match"
    )

    log.info(
        f"Atom type assignment: {len(atoms)} atoms, {len(atom_type_params)} unique types, {no_match} no match"
    )
    log.info(f"  hop3 matches: {hop3_matches}")
    log.info(f"  hop2 matches: {hop2_matches}")
    log.info(f"  hop1 matches: {hop1_matches}")
    log.info(f"  hop0 matches: {hop0_matches}")
    log.info(f"  no match: {no_match}")

    # Build output data
    # Get box dimensions from coordinates
    all_coords = list(coords.values())
    xvals = [c[0] for c in all_coords]
    yvals = [c[1] for c in all_coords]
    zvals = [c[2] for c in all_coords]
    xlo, xhi = min(xvals) - 5, max(xvals) + 5
    ylo, yhi = min(yvals) - 5, max(yvals) + 5
    zlo, zhi = min(zvals) - 5, max(zvals) + 5

    # Build LammpsData object
    lmp_data = LammpsData()
    lmp_data.header_comment = "LAMMPS data file Generated by mappingff"
    lmp_data.atoms = len(atoms)
    lmp_data.bonds = len(bond_type_map)
    lmp_data.angles = len(angle_type_map)
    lmp_data.dihedrals = len(dihedral_type_map)
    lmp_data.impropers = len(improper_type_map)
    lmp_data.atom_types = len(atom_type_params)
    lmp_data.bond_types = len(bond_type_params)
    lmp_data.angle_types = len(angle_type_params)
    lmp_data.dihedral_types = len(dihedral_type_params)
    lmp_data.improper_types = len(improper_type_params)
    lmp_data.xlo = xlo
    lmp_data.xhi = xhi
    lmp_data.ylo = ylo
    lmp_data.yhi = yhi
    lmp_data.zlo = zlo
    lmp_data.zhi = zhi

    # Masses
    for _, params in sorted(atom_type_params.items(), key=lambda x: x[1]["type"]):
        lmp_data.masses.append((params["type"], params["mass"]))

    # Pair Coeffs
    for _, params in sorted(atom_type_params.items(), key=lambda x: x[1]["type"]):
        lmp_data.pair_coeffs.append(
            (params["type"], params["epsilon"], params["sigma"])
        )

    # Bond Coeffs
    for btid in sorted(bond_type_params.keys()):
        k, r0 = bond_type_params[btid]
        lmp_data.bond_coeffs.append((btid, k, r0))

    # Angle Coeffs
    for atid in sorted(angle_type_params.keys()):
        k, theta0 = angle_type_params[atid]
        lmp_data.angle_coeffs.append((atid, k, theta0))

    # Dihedral Coeffs
    for dtid in sorted(dihedral_type_params.keys()):
        coeffs = dihedral_type_params[dtid]
        lmp_data.dihedral_coeffs.append((dtid,) + coeffs)

    # Improper Coeffs
    for itid in sorted(improper_type_params.keys()):
        coeffs = improper_type_params[itid]  # type: ignore[assignment]
        lmp_data.improper_coeffs.append((itid,) + coeffs)  # type: ignore[arg-type]

    # Atom records
    for atom in atoms:
        atom_idx = atom["idx"]
        lammps_type = atom_type_map[atom_idx]
        output_type = atom_type_params[lammps_type]["type"]
        charge = atom_type_params[lammps_type].get("charge", 0.0)
        x, y, z = coords[atom_idx]
        lmp_data.atom_records.append((atom_idx, 1, output_type, charge, x, y, z))

    # Bond records
    for bond in bonds:
        bond_idx = bond["idx"]
        if bond_idx in bond_type_map:
            bt = bond_type_map[bond_idx]
            a1 = bond["a1"]
            a2 = bond["a2"]
            lmp_data.bond_records.append((bond_idx, bt, a1, a2))

    # Angle records
    for angle in angles:
        angle_idx = angle["idx"]
        if angle_idx in angle_type_map:
            at = angle_type_map[angle_idx]
            a1 = angle["a1"]
            a2 = angle["a2"]
            a3 = angle["a3"]
            lmp_data.angle_records.append((angle_idx, at, a1, a2, a3))

    # Dihedral records
    for dihedral in dihedrals:
        dih_idx = dihedral["idx"]
        if dih_idx in dihedral_type_map:
            dt = dihedral_type_map[dih_idx]
            a1 = dihedral["a1"]
            a2 = dihedral["a2"]
            a3 = dihedral["a3"]
            a4 = dihedral["a4"]
            lmp_data.dihedral_records.append((dih_idx, dt, a1, a2, a3, a4))

    # Improper records
    for improper in impropers:
        imp_idx = improper["idx"]
        if imp_idx in improper_type_map:
            it = improper_type_map[imp_idx]
            a1 = improper["a1"]
            a2 = improper["a2"]
            a3 = improper["a3"]
            a4 = improper["a4"]
            lmp_data.improper_records.append((imp_idx, it, a1, a2, a3, a4))

    # Adjust total charge if needed
    before_charge = sum(atom[3] for atom in lmp_data.atom_records)
    if total_charge is not None:
        after_step1_charge, after_step2_charge = adjust_total_charge(
            lmp_data, total_charge, db, atom_type_params
        )
        log.info(
            f"Total charge: {after_step2_charge:.6f} "
            f"(before: {before_charge:.6f}, "
            f"after step1: {after_step1_charge:.6f}, "
            f"after step2: {after_step2_charge:.6f}, "
            f"target: {total_charge})"
        )
    else:
        log.info(f"Total charge: {before_charge:.6f}")

    # Write LAMMPS file
    generate_lammps(lmp_data, out_path)
    log.info(f"Output written to {out_path}")

    return {
        "atoms": len(atoms),
        "bonds": len(bonds),
        "angles": len(angle_type_map),
        "dihedrals": len(dihedral_type_map),
        "impropers": len(improper_type_map),
        "unique_types": len(atom_type_params),
        "hop3_matches": hop3_matches,
        "hop2_matches": hop2_matches,
        "hop1_matches": hop1_matches,
        "hop0_matches": hop0_matches,
        "no_match": no_match,
    }
