"""Parameterization pipeline that maps atoms/bonded terms and writes LAMMPS output."""

from pathlib import Path

from rdkit import RDLogger

from rdkit import Chem
from macromapff.domain.atom_match import build_atom_match
from macromapff.domain.bonded_match import assign_bonded_params
from macromapff.domain.bonded_match import build_type_map
from macromapff.domain.env_key_match import INTERNAL_HOP_DEPTH
from macromapff.domain.term_enumeration import enumerate_terms
from macromapff.io.input import load_bonded_db
from macromapff.io.input import load_hop_param_db
from macromapff.io.input import load_input_structure
from macromapff.io.log import append_build_log
from macromapff.io.log import init_build_log
from macromapff.io.log import log_bonded_match
from macromapff.io.log import write_missing_env_log
from macromapff.io.output import write_atom_keytype_map
from macromapff.io.output import write_lammps_data as _write_lammps_data

INTERNAL_FALLBACK_HOPS = (1, 0)
INTERNAL_BOX_PADDING = 20.0
INTERNAL_MOLECULE_ID = 1
INTERNAL_BUILD_LOG_NAME = "parameterize.log"
INTERNAL_ATOM_KEYTYPE_MAP_NAME = "atom_index_key_types.csv"
DB_HOP2_ENV = "hop_env/hop2_KeyMap.csv"
DB_HOP1_ENV = "hop_env/hop1_KeyMap.csv"
DB_HOP0_ENV = "hop_env/hop0_KeyMap.csv"
DB_BONDED = "Global_BondedTerms.csv"


def build_atom_match_with_logs(
    mol: Chem.Mol,
    hop2_env_to_atom_param: dict,
    hop1_env_to_atom_param: dict,
    hop0_env_to_atom_param: dict,
    structure_path: Path,
    build_log_path: Path,
):
    """Assign atom-level match results and record fallback/missing diagnostics."""
    atom_records, local_type_params, fallback_hit_counter, missing = build_atom_match(
        mol,
        hop2_env_to_atom_param,
        hop1_env_to_atom_param,
        hop0_env_to_atom_param,
        fallback_hops=INTERNAL_FALLBACK_HOPS,
    )

    if missing:
        write_missing_env_log(
            log_path=build_log_path,
            structure_path=structure_path,
            missing=missing,
            hop_depth=INTERNAL_HOP_DEPTH,
            fallback_hops=INTERNAL_FALLBACK_HOPS,
        )

    append_build_log(
        build_log_path,
        [
            "[ATOM] build_atom_match finished",
            f"[ATOM] total_atoms={len(atom_records)}",
            f"[ATOM] missing_atoms={len(missing)}",
            "[ATOM] fallback_hits="
            + ", ".join(f"hop{h}:{c}" for h, c in sorted(fallback_hit_counter.items())),
        ],
    )

    return atom_records, local_type_params, dict(sorted(fallback_hit_counter.items()))


def assign_bonded_params_with_logs(
    kind: str,
    terms,
    atom_key_sets,
    matcher_indexes,
    build_log_path: Path,
):
    """Assign one interaction kind's multi-atom parameters and emit logs."""
    (
        idx_rev_patterns,
        idx_imp_patterns,
        idx_rev_inverted,
        idx_imp_center_inverted,
    ) = matcher_indexes

    (
        records,
        missing,
        ambiguous,
        cache_hit,
        cache_miss,
        cache_size,
    ) = assign_bonded_params(
        kind,
        terms,
        atom_key_sets,
        idx_rev_patterns,
        idx_imp_patterns,
        idx_rev_inverted,
        idx_imp_center_inverted,
    )

    log_bonded_match(
        kind,
        missing,
        ambiguous,
        cache_hit,
        cache_miss,
        cache_size,
        build_log_path=build_log_path,
    )

    return records, missing, ambiguous

def parameterize_lammps(db_dir: Path, structure: Path, out: Path):
    """Run atom matching, bonded matching, and LAMMPS file generation."""
    db_root = db_dir.expanduser().resolve()
    RDLogger.DisableLog("rdApp.warning")
    hop2_env = db_root / DB_HOP2_ENV
    hop1_env = db_root / DB_HOP1_ENV
    hop0_env = db_root / DB_HOP0_ENV
    bonded = db_root / DB_BONDED

    build_log_path = out.parent / INTERNAL_BUILD_LOG_NAME
    atom_keytype_map_path = out.parent / INTERNAL_ATOM_KEYTYPE_MAP_NAME

    init_build_log(
        build_log_path,
        [
            "build log",
            f"structure: {structure}",
            f"db_dir: {db_root}",
            f"hop_depth: {INTERNAL_HOP_DEPTH}",
            "fallback_hops: " + ",".join(str(x) for x in INTERNAL_FALLBACK_HOPS),
        ],
    )

    try:
        mol = load_input_structure(structure)
        if mol is None:
            raise ValueError(f"RDKit failed to read structure file: {structure}")
        if mol.GetNumConformers() == 0:
            raise ValueError("Input structure has no 3D conformer coordinates")
        hop2_env_to_atom_param = load_hop_param_db(hop2_env)
        hop1_env_to_atom_param = load_hop_param_db(hop1_env)
        hop0_env_to_atom_param = load_hop_param_db(hop0_env)

        matcher_indexes = load_bonded_db(bonded)

        atom_records, atom_type_rows, fallback_hit_counter = build_atom_match_with_logs(
            mol,
            hop2_env_to_atom_param,
            hop1_env_to_atom_param,
            hop0_env_to_atom_param,
            structure_path=structure,
            build_log_path=build_log_path,
        )
        write_atom_keytype_map(atom_keytype_map_path, atom_records)
        append_build_log(
            build_log_path,
            [f"[ATOM] key_type_map={atom_keytype_map_path}"],
        )
        atom_key_sets = {r["atom_id"]: set(r["global_key_ids"]) for r in atom_records}

        bonds, angles, dihedrals, impropers = enumerate_terms(mol)

        bond_records, bond_missing, bond_amb = assign_bonded_params_with_logs(
            "bond",
            bonds,
            atom_key_sets,
            matcher_indexes,
            build_log_path=build_log_path,
        )
        angle_records, angle_missing, angle_amb = assign_bonded_params_with_logs(
            "angle",
            angles,
            atom_key_sets,
            matcher_indexes,
            build_log_path=build_log_path,
        )
        dihedral_records, dihedral_missing, dihedral_amb = assign_bonded_params_with_logs(
            "dihedral",
            dihedrals,
            atom_key_sets,
            matcher_indexes,
            build_log_path=build_log_path,
        )
        improper_records, improper_missing, improper_amb = assign_bonded_params_with_logs(
            "improper",
            impropers,
            atom_key_sets,
            matcher_indexes,
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

        _write_lammps_data(
            out_path=out,
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
            molecule_id=INTERNAL_MOLECULE_ID,
            box_padding=INTERNAL_BOX_PADDING,
        )
    except Exception as exc:
        append_build_log(build_log_path, [f"[ERROR] {type(exc).__name__}: {exc}"])
        raise

    return {
        "atoms": len(atom_records),
        "bonds": len(bond_records),
        "angles": len(angle_records),
        "dihedrals": len(dihedral_records),
        "impropers": len(improper_records),
        "impropers_missing": len(improper_missing),
        "fallback_hits": fallback_hit_counter,
    }


