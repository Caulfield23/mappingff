from pathlib import Path

from rdkit import RDLogger

from rdkit import Chem
from macromapff.domain import assign_multiatom_params_core
from macromapff.domain import build_atom_types_core
from macromapff.domain import build_type_map
from macromapff.domain import enumerate_terms
from macromapff.io import append_build_log
from macromapff.io import init_build_log
from macromapff.io import load_hop_param_db
from macromapff.io import load_input_structure
from macromapff.io import load_multiatom_db
from macromapff.io import log_multiatom_match
from macromapff.io import write_atom_keytype_map
from macromapff.io import write_lammps_data as _write_lammps_data
from macromapff.io import write_missing_env_log

INTERNAL_HOP_DEPTH = 2
INTERNAL_FALLBACK_HOPS = (1, 0)
DB_FINAL_ENV = "final_env_keymap.csv"
DB_HOP2_ENV = "hop_env/hop2_env_keymap.csv"
DB_HOP1_ENV = "hop_env/hop1_env_keymap.csv"
DB_HOP0_ENV = "hop_env/hop0_env_keymap.csv"
DB_MULTIATOM = "multiatom_master_keytype.csv"


def build_atom_types(
    mol: Chem.Mol,
    hop2_env_to_atom_param: dict,
    hop1_env_to_atom_param: dict,
    hop0_env_to_atom_param: dict,
    hop_depth: int,
    fallback_hops=(),
    missing_log_path: Path = None,
    structure_path: Path = None,
    build_log_path: Path = None,
):
    atom_records, local_type_params, fallback_hit_counter, missing = build_atom_types_core(
        mol,
        hop2_env_to_atom_param,
        hop1_env_to_atom_param,
        hop0_env_to_atom_param,
        hop_depth,
        fallback_hops=fallback_hops,
    )

    if missing and missing_log_path is not None:
        write_missing_env_log(
            log_path=missing_log_path,
            structure_path=structure_path or Path("unknown"),
            missing=missing,
            hop_depth=hop_depth,
            fallback_hops=fallback_hops,
        )

    append_build_log(
        build_log_path,
        [
            "[ATOM] build_atom_types finished",
            f"[ATOM] total_atoms={len(atom_records)}",
            f"[ATOM] missing_atoms={len(missing)}",
            "[ATOM] fallback_hits="
            + ", ".join(f"hop{h}:{c}" for h, c in sorted(fallback_hit_counter.items())),
        ],
    )

    return atom_records, local_type_params, dict(sorted(fallback_hit_counter.items()))


def assign_multiatom_params(
    kind: str,
    terms,
    atom_key_sets,
    idx_rev_patterns,
    idx_imp_patterns,
    idx_rev_inverted,
    idx_imp_center_inverted,
    strict_missing=True,
    build_log_path: Path = None,
):
    (
        records,
        missing,
        ambiguous,
        cache_hit,
        cache_miss,
        cache_size,
    ) = assign_multiatom_params_core(
        kind,
        terms,
        atom_key_sets,
        idx_rev_patterns,
        idx_imp_patterns,
        idx_rev_inverted,
        idx_imp_center_inverted,
    )

    log_multiatom_match(
        kind,
        missing,
        ambiguous,
        cache_hit,
        cache_miss,
        cache_size,
        build_log_path=build_log_path,
        strict_missing=strict_missing,
    )

    return records, missing, ambiguous

def _resolve_db_paths(db_dir: Path):
    base = db_dir.expanduser().resolve()
    return {
        "final_env": base / DB_FINAL_ENV,
        "hop2_env": base / DB_HOP2_ENV,
        "hop1_env": base / DB_HOP1_ENV,
        "hop0_env": base / DB_HOP0_ENV,
        "multiatom": base / DB_MULTIATOM,
    }


class LammpsGenerator:
    def __init__(
        self,
        db_dir: Path,
        box_padding: float = 20.0,
        molecule_id: int = 1,
        build_log: Path | None = None,
        atom_keytype_map: Path | None = None,
    ) -> None:
        self.db_dir = db_dir
        self.box_padding = box_padding
        self.molecule_id = molecule_id
        self.build_log = build_log
        self.atom_keytype_map = atom_keytype_map

    def generate(self, structure: Path, out: Path):
        RDLogger.DisableLog("rdApp.warning")
        db_paths = _resolve_db_paths(self.db_dir)

        build_log_path = (
            self.build_log if self.build_log is not None else out.parent / "build.log"
        )
        atom_keytype_map_path = (
            self.atom_keytype_map
            if self.atom_keytype_map is not None
            else out.parent / "atom_index_key_types.csv"
        )

        init_build_log(
            build_log_path,
            [
                "build log",
                f"structure: {structure}",
                f"db_dir: {self.db_dir.expanduser().resolve()}",
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
                molecule_id=self.molecule_id,
                box_padding=self.box_padding,
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


