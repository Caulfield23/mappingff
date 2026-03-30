#!/usr/bin/env python3
import argparse
from pathlib import Path

from rdkit import RDLogger

try:
    from .core.atom_match import append_build_log
    from .core.atom_match import build_atom_types
    from .core.atom_match import init_build_log
    from .core.atom_match import load_hop_param_db
    from .core.atom_match import load_input_structure
    from .core.atom_match import write_atom_keytype_map
    from .core.lammps_write import write_lammps_data
    from .core.multi_match import assign_multiatom_params
    from .core.multi_match import build_type_map
    from .core.multi_match import enumerate_terms
    from .core.multi_match import load_multiatom_db
except ImportError:
    from core.atom_match import append_build_log
    from core.atom_match import build_atom_types
    from core.atom_match import init_build_log
    from core.atom_match import load_hop_param_db
    from core.atom_match import load_input_structure
    from core.atom_match import write_atom_keytype_map
    from core.lammps_write import write_lammps_data
    from core.multi_match import assign_multiatom_params
    from core.multi_match import build_type_map
    from core.multi_match import enumerate_terms
    from core.multi_match import load_multiatom_db


INTERNAL_HOP_DEPTH = 2
INTERNAL_FALLBACK_HOPS = (1, 0)
DB_FINAL_ENV = "final_env_keymap.csv"
DB_HOP2_ENV = "hop2_env_keymap.csv"
DB_HOP1_ENV = "hop1_env_keymap.csv"
DB_HOP0_ENV = "hop0_env_keymap.csv"
DB_MULTIATOM = "multiatom_master_keytype.csv"


def _resolve_db_paths(db_dir: Path):
    base = db_dir.expanduser().resolve()
    return {
        "final_env": base / DB_FINAL_ENV,
        "hop2_env": base / DB_HOP2_ENV,
        "hop1_env": base / DB_HOP1_ENV,
        "hop0_env": base / DB_HOP0_ENV,
        "multiatom": base / DB_MULTIATOM,
    }


def main():
    RDLogger.DisableLog("rdApp.warning")

    parser = argparse.ArgumentParser(
        description=(
            "Read a full molecule (.mol only), match parameters using final_env_keymap "
            "and multiatom_master databases, and write a complete LAMMPS data file."
        )
    )
    parser.add_argument(
        "--structure", required=True, type=Path, help="Input full-molecule .mol file"
    )
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=Path("database"),
        help="Database directory containing fixed internal CSV names",
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="Output LAMMPS data file path"
    )
    parser.add_argument(
        "--box-padding", type=float, default=20.0, help="Simulation box padding (default: 20.0)"
    )
    parser.add_argument(
        "--molecule-id", type=int, default=1, help="Molecule id in the Atoms section (default: 1)"
    )
    parser.add_argument(
        "--build-log",
        type=Path,
        default=None,
        help="Build-log output path (default: <out_dir>/build.log)",
    )
    parser.add_argument(
        "--atom-keytype-map",
        type=Path,
        default=None,
        help="Output map of atom index to key_type (default: <out_dir>/atom_index_key_types.csv)",
    )
    args = parser.parse_args()
    db_paths = _resolve_db_paths(args.db_dir)

    build_log_path = (
        args.build_log if args.build_log is not None else args.out.parent / "build.log"
    )
    atom_keytype_map_path = (
        args.atom_keytype_map
        if args.atom_keytype_map is not None
        else args.out.parent / "atom_index_key_types.csv"
    )

    init_build_log(
        build_log_path,
        [
            "build log",
            f"structure: {args.structure}",
            f"db_dir: {args.db_dir.expanduser().resolve()}",
            f"hop_depth: {INTERNAL_HOP_DEPTH}",
            "fallback_hops: " + ",".join(str(x) for x in INTERNAL_FALLBACK_HOPS),
        ],
    )

    try:
        mol = load_input_structure(args.structure)
        if mol is None:
            raise ValueError(f"RDKit failed to read structure file: {args.structure}")
        if mol.GetNumConformers() == 0:
            raise ValueError("Input structure has no 3D conformer coordinates")

        hop2_env_to_atom_param = load_hop_param_db(db_paths["hop2_env"])
        hop1_env_to_atom_param = load_hop_param_db(db_paths["hop1_env"])
        hop0_env_to_atom_param = load_hop_param_db(db_paths["hop0_env"])
        idx_rev, idx_imp, idx_rev_inv, idx_imp_center_inv = load_multiatom_db(
            db_paths["multiatom"]
        )

        missing_log_path = build_log_path
        atom_records, atom_type_rows, fallback_hit_counter = build_atom_types(
            mol,
            hop2_env_to_atom_param,
            hop1_env_to_atom_param,
            hop0_env_to_atom_param,
            INTERNAL_HOP_DEPTH,
            fallback_hops=INTERNAL_FALLBACK_HOPS,
            missing_log_path=missing_log_path,
            structure_path=args.structure,
            build_log_path=build_log_path,
        )
        write_atom_keytype_map(atom_keytype_map_path, atom_records)
        append_build_log(
            build_log_path,
            [f"[ATOM] key_type_map={atom_keytype_map_path}"],
        )
        atom_key_sets = {r["atom_id"]: set(r["global_key_ids"]) for r in atom_records}

        bonds, angles, dihedrals, impropers = enumerate_terms(mol)

        bond_records, bond_missing, bond_amb = assign_multiatom_params(
            "bond",
            bonds,
            atom_key_sets,
            idx_rev,
            idx_imp,
            idx_rev_inv,
            idx_imp_center_inv,
            strict_missing=False,
            build_log_path=build_log_path,
        )
        angle_records, angle_missing, angle_amb = assign_multiatom_params(
            "angle",
            angles,
            atom_key_sets,
            idx_rev,
            idx_imp,
            idx_rev_inv,
            idx_imp_center_inv,
            strict_missing=False,
            build_log_path=build_log_path,
        )
        dihedral_records, dihedral_missing, dihedral_amb = assign_multiatom_params(
            "dihedral",
            dihedrals,
            atom_key_sets,
            idx_rev,
            idx_imp,
            idx_rev_inv,
            idx_imp_center_inv,
            strict_missing=False,
            build_log_path=build_log_path,
        )
        improper_records, improper_missing, improper_amb = assign_multiatom_params(
            "improper",
            impropers,
            atom_key_sets,
            idx_rev,
            idx_imp,
            idx_rev_inv,
            idx_imp_center_inv,
            strict_missing=False,
            build_log_path=build_log_path,
        )

        append_build_log(
            build_log_path,
            [
                f"[BOND] matched={len(bond_records)} missing={len(bond_missing)} ambiguous={len(bond_amb)}",
                f"[ANGLE] matched={len(angle_records)} missing={len(angle_missing)} ambiguous={len(angle_amb)}",
                f"[DIHEDRAL] matched={len(dihedral_records)} missing={len(dihedral_missing)} ambiguous={len(dihedral_amb)}",
                f"[IMPROPER] matched={len(improper_records)} missing={len(improper_missing)} ambiguous={len(improper_amb)}",
            ],
        )

        bond_type_rows = build_type_map(bond_records)
        angle_type_rows = build_type_map(angle_records)
        dihedral_type_rows = build_type_map(dihedral_records)
        improper_type_rows = build_type_map(improper_records)

        write_lammps_data(
            out_path=args.out,
            mol=mol,
            atom_records=atom_records,
            atom_type_rows=atom_type_rows,
            bond_records=bond_records,
            angle_records=angle_records,
            dihedral_records=dihedral_records,
            improper_records=improper_records,
            bond_type_rows=bond_type_rows,
            angle_type_rows=angle_type_rows,
            dihedral_type_rows=dihedral_type_rows,
            improper_type_rows=improper_type_rows,
            molecule_id=args.molecule_id,
            box_padding=args.box_padding,
        )
    except Exception as exc:
        append_build_log(build_log_path, [f"[ERROR] {type(exc).__name__}: {exc}"])
        raise

    print("Done:")
    print(f"- output: {args.out}")
    print(f"- atoms: {len(atom_records)}")
    if fallback_hit_counter:
        detail = ", ".join(
            f"hop{hop}:{count}" for hop, count in fallback_hit_counter.items()
        )
        print(f"- atom env_key fallback hits: {detail}")
    else:
        print("- atom env_key fallback hits: 0")
    print(f"- bonds: {len(bond_records)}")
    print(f"- angles: {len(angle_records)}")
    print(f"- dihedrals: {len(dihedral_records)}")
    print(f"- impropers (matched): {len(improper_records)}")
    print(f"- impropers (removed due to no candidate match): {len(improper_missing)}")


if __name__ == "__main__":
    main()
