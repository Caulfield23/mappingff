#!/usr/bin/env python3
"""Build atom-level environment mappings from structure and LAMMPS data."""

from pathlib import Path

from rdkit import Chem

from macromapff.domain import EnvFeatureBuilder
from macromapff.io import load_structure_any
from macromapff.io import parse_lammps_data
from macromapff.io import write_atom_env_csv


INTERNAL_HOP_DEPTH = 2


def build_mapping(
    structure_path: Path,
    out_dir: Path,
    module: str,
    lmp_path: Path = None,
    hop_depth: int = 2,
):
    """Build one module's atom-env CSV from structure and LAMMPS data."""
    if lmp_path is None:
        raise ValueError("--lmp is required for the pure LAMMPS route")

    lmp_atoms, _ = parse_lammps_data(lmp_path)

    def _build_with_loaded_mol(mol: Chem.Mol, source_path: Path):
        """Construct atom-env rows from an already loaded RDKit molecule."""
        if mol.GetNumAtoms() != len(lmp_atoms):
            raise ValueError(
                f"Atom count mismatch: structure={mol.GetNumAtoms()} vs lammps_atoms={len(lmp_atoms)}."
                f" structure={source_path}"
            )

        atom_rows = []
        env_builder = EnvFeatureBuilder()
        atom_context = env_builder.precompute_atom_context(mol)

        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom_id = atom_idx + 1
            lmp_atom = lmp_atoms.get(atom_id)
            if lmp_atom is None:
                raise ValueError(f"Missing atom_id={atom_id} in LAMMPS Atoms section")

            type_id = lmp_atom["lmp_type"]
            if lmp_atom["atomic_num"] != atom.GetAtomicNum():
                raise ValueError(
                    f"Atomic element mismatch: atom_index={atom_id}, structure={source_path}, "
                    f"structure Z={atom.GetAtomicNum()} vs lammps_mass Z={lmp_atom['atomic_num']}"
                )

            charge = lmp_atom["charge"]
            sigma = lmp_atom["sigma"]
            epsilon = lmp_atom["epsilon"]

            env_key, _ = env_builder.make_env_key(
                mol,
                atom,
                hop_depth=hop_depth,
                atom_ctx=atom_context.get(atom_idx),
            )
            row = {
                "module": module,
                "atom_index": atom_idx + 1,
                "atom_name": atom.GetSymbol(),
                "opls_type_id": type_id,
                "opls_type_name": lmp_atom["type_name"],
                "charge": charge,
                "sigma": sigma,
                "epsilon": epsilon,
                "env_key": env_key,
            }
            atom_rows.append(row)
        return write_atom_env_csv(out_dir, module, atom_rows)

    mol = load_structure_any(structure_path)
    try:
        return _build_with_loaded_mol(mol, structure_path)
    except Exception as primary_exc:
        if structure_path.suffix.lower() != ".mol":
            raise

        pdb_fallback = structure_path.with_suffix(".pdb")
        if not pdb_fallback.exists():
            raise ValueError(
                f"Failed to build from mol and no fallback pdb found: {pdb_fallback}\n"
                f"Original error: {primary_exc}"
            )

        try:
            mol_pdb = load_structure_any(pdb_fallback)
            result = _build_with_loaded_mol(mol_pdb, pdb_fallback)
            print(
                f"[INFO] Building from mol failed; switched to pdb automatically: {pdb_fallback}",
                flush=True,
            )
            return result
        except Exception as pdb_exc:
            raise ValueError(
                "Both mol and pdb build attempts failed.\n"
                f"- mol: {structure_path}\n"
                f"  error: {primary_exc}\n"
                f"- pdb: {pdb_fallback}\n"
                f"  error: {pdb_exc}"
            )


