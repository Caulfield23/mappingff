#!/usr/bin/env python3
"""Build atom-level environment mappings from structure and LAMMPS data.

Outputs:
- ``{module}_atom_env.csv``
"""

from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from macromapff.domain import EnvFeatureBuilder
from macromapff.io import load_structure_any
from macromapff.io import parse_lammps_data
from macromapff.io import write_atom_env_csv


INTERNAL_HOP_DEPTH = 2


def _compile_smarts_patterns():
    pattern_defs = [
        ("amide_n", "[NX3][CX3](=[OX1])[#6]"),
        ("amide_c", "[CX3](=[OX1])[NX3]"),
        ("carboxylate", "[CX3](=O)[O-]"),
        ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
        ("ester", "[CX3](=O)[OX2][#6]"),
        ("ether_o", "[OD2]([#6])[#6]"),
        ("hydroxyl_o", "[OX2H]"),
        ("carbonyl_c", "[CX3]=[OX1]"),
        ("carbonyl_o", "[OX1]=[CX3]"),
        ("amine_n", "[NX3;H2,H1;!$(NC=O)]"),
        ("quaternary_n", "[NX4+]"),
        ("nitrile_c", "[CX2]#N"),
        ("nitrile_n", "N#[CX2]"),
        ("sulfone_s", "S(=O)(=O)"),
        ("sulfoxide_s", "[SX3](=O)"),
        ("phosphate_p", "P(=O)(O)O"),
        ("aromatic_atom", "[a]"),
        ("silicon_atom", "[Si]"),
        ("siloxane_o", "[O;X2]-[Si]-[O;X2]"),
    ]

    compiled = []
    for name, smarts in pattern_defs:
        q = Chem.MolFromSmarts(smarts)
        if q is not None:
            compiled.append((name, q))
    return compiled


SMARTS_PATTERNS = _compile_smarts_patterns()




def bond_stereo_code(bond: Chem.Bond) -> str:
    stereo = bond.GetStereo()
    if stereo == Chem.rdchem.BondStereo.STEREOE:
        return "E"
    if stereo == Chem.rdchem.BondStereo.STEREOZ:
        return "Z"
    if stereo == Chem.rdchem.BondStereo.STEREOANY:
        return "ANY"
    return "NONE"


def atom_smarts_hits_map(mol: Chem.Mol):
    hit_map = {idx: [] for idx in range(mol.GetNumAtoms())}
    for name, query in SMARTS_PATTERNS:
        for match in mol.GetSubstructMatches(query, uniquify=True):
            for atom_idx in match:
                hit_map[atom_idx].append(name)
    return {idx: sorted(set(tags)) for idx, tags in hit_map.items()}


def atom_bond_order_hist(atom: Chem.Atom):
    env_builder = EnvFeatureBuilder()
    hist = Counter()
    for bond in atom.GetBonds():
        hist[env_builder.bond_type_code(bond)] += 1
    return dict(sorted(hist.items(), key=lambda x: x[0]))


def get_ring_sizes(atom: Chem.Atom, min_size: int = 3, max_size: int = 8):
    sizes = []
    for size in range(min_size, max_size + 1):
        if atom.IsInRingSize(size):
            sizes.append(size)
    return sizes


def second_shell_element_counts(atom: Chem.Atom):
    center_idx = atom.GetIdx()
    first_shell = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
    counts = Counter()
    for nbr in atom.GetNeighbors():
        for nbr2 in nbr.GetNeighbors():
            idx2 = nbr2.GetIdx()
            if idx2 == center_idx or idx2 in first_shell:
                continue
            counts[nbr2.GetAtomicNum()] += 1
    return dict(sorted(counts.items(), key=lambda x: x[0]))


def _path_bond_signature(mol: Chem.Mol, src_idx: int, dst_idx: int):
    path = Chem.GetShortestPath(mol, src_idx, dst_idx)
    if len(path) <= 1:
        return ""
    env_builder = EnvFeatureBuilder()
    sig = []
    for left, right in zip(path[:-1], path[1:]):
        bond = mol.GetBondBetweenAtoms(left, right)
        sig.append(env_builder.bond_type_code(bond))
    return "-".join(sig)


def atom_env_fragment_smiles(mol: Chem.Mol, atom_idx: int, radius: int):
    env_bonds = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
    atom_ids = {atom_idx}
    for bidx in env_bonds:
        bond = mol.GetBondWithIdx(bidx)
        atom_ids.add(bond.GetBeginAtomIdx())
        atom_ids.add(bond.GetEndAtomIdx())

    if not atom_ids:
        atom = mol.GetAtomWithIdx(atom_idx)
        return f"{atom.GetAtomicNum()}"

    frag = Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(atom_ids),
        bondsToUse=list(env_bonds),
        isomericSmiles=True,
        canonical=True,
    )
    return frag


def build_mapping(
    structure_path: Path,
    out_dir: Path,
    module: str,
    lmp_path: Path = None,
    hop_depth: int = 2,
):
    if lmp_path is None:
        raise ValueError("--lmp is required for the pure LAMMPS route")

    lmp_atoms, _ = parse_lammps_data(lmp_path)

    def _build_with_loaded_mol(mol: Chem.Mol, source_path: Path):
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


