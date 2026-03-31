#!/usr/bin/env python3
"""Build atom-level environment mappings from structure and LAMMPS data.

Outputs:
- ``{module}_atom_env.csv``
"""

import csv
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from macromapff.pipeline.core.lammps_parse import parse_lammps_sections


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


ENV_KEY_PRIORITY = [
    "z",
    "formal_charge",
    "aromatic",
    "hybridization",
    "degree",
    "total_hs",
    "in_ring",
    "ring_count",
    "neighbor_sig",
    "bond_kinds",
]


def ordered_env_key_obj(obj: dict):
    if not isinstance(obj, dict):
        return obj

    out = {}
    for key in ENV_KEY_PRIORITY:
        if key in obj:
            out[key] = obj[key]

    for key in sorted(k for k in obj.keys() if k not in out):
        out[key] = obj[key]

    return out


def infer_atomic_num_from_mass(mass: float, max_z: int = 36, tol: float = 0.2):
    ptable = Chem.GetPeriodicTable()
    best = None
    for z in range(1, max_z + 1):
        w = ptable.GetAtomicWeight(z)
        d = abs(w - mass)
        if best is None or d < best[0]:
            best = (d, z)
    if best is None or best[0] > tol:
        raise ValueError(
            f"Failed to infer element from mass {mass} (min delta {best[0] if best else 'N/A'})"
        )
    return best[1]


def parse_lammps_data(lmp_path: Path):
    masses, pair_coeffs, atoms = parse_lammps_sections(lmp_path)

    if not masses:
        raise ValueError(f"No Masses section parsed from: {lmp_path}")
    if not pair_coeffs:
        raise ValueError(f"No Pair Coeffs section parsed from: {lmp_path}")
    if not atoms:
        raise ValueError(f"No Atoms section parsed from: {lmp_path}")

    type_info = {}
    ptable = Chem.GetPeriodicTable()
    for type_id, mass in masses.items():
        atomic_num = infer_atomic_num_from_mass(mass)
        type_info[type_id] = {
            "mass": mass,
            "atomic_num": atomic_num,
            "type_name": f"LMP_{type_id}",
            "element": ptable.GetElementSymbol(atomic_num),
        }

    for atom_id, info in atoms.items():
        lmp_type = info["lmp_type"]
        if lmp_type not in pair_coeffs:
            raise ValueError(f"LAMMPS atom type {lmp_type} is missing Pair Coeffs")
        if lmp_type not in type_info:
            raise ValueError(f"LAMMPS atom type {lmp_type} is missing Masses")
        info["sigma"] = pair_coeffs[lmp_type]["sigma"]
        info["epsilon"] = pair_coeffs[lmp_type]["epsilon"]
        info["atomic_num"] = type_info[lmp_type]["atomic_num"]
        info["type_name"] = type_info[lmp_type]["type_name"]

    return atoms, type_info


def load_structure(structure_path: Path):
    suffix = structure_path.suffix.lower()

    def _sanitize_or_raise(mol, source_path: Path, label: str):
        if mol is None:
            raise ValueError(f"RDKit failed to read {label} file: {source_path}")
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            raise ValueError(f"RDKit read {label} but sanitization failed: {source_path} ({exc})")
        return mol

    def _load_pdb(pdb_path: Path):
        if not pdb_path.exists():
            return None
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
        if mol is None:
            return None
        return _sanitize_or_raise(mol, pdb_path, "pdb")

    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(structure_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol

        mol = Chem.MolFromMolFile(str(structure_path), removeHs=False, sanitize=False)
        if mol is not None:
            return _sanitize_or_raise(mol, structure_path, "mol")

        pdb_fallback = structure_path.with_suffix(".pdb")
        mol = _load_pdb(pdb_fallback)
        if mol is not None:
            print(
                f"[INFO] Failed to read mol; automatically fell back to pdb input: {pdb_fallback}",
                flush=True,
            )
            return mol

        raise ValueError(
            f"RDKit failed to read mol file: {structure_path}; fallback to same-name pdb also failed: {pdb_fallback}"
        )

    if suffix == ".pdb":
        mol = _load_pdb(structure_path)
        if mol is None:
            raise ValueError(f"RDKit failed to read pdb file: {structure_path}")
        return mol

    if suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(structure_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol

        mol = Chem.MolFromMol2File(str(structure_path), removeHs=False, sanitize=False)
        if mol is not None:
            return _sanitize_or_raise(mol, structure_path, "mol2")

        return mol

    raise ValueError(
        f"Unsupported structure format: {structure_path}. Supported formats: .mol/.mol2/.pdb"
    )


def bond_type_code(bond: Chem.Bond) -> str:
    bt = bond.GetBondType()
    if bt == Chem.rdchem.BondType.SINGLE:
        return "S"
    if bt == Chem.rdchem.BondType.DOUBLE:
        return "D"
    if bt == Chem.rdchem.BondType.TRIPLE:
        return "T"
    if bt == Chem.rdchem.BondType.AROMATIC:
        return "A"
    return str(bt)


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
    hist = Counter()
    for bond in atom.GetBonds():
        hist[bond_type_code(bond)] += 1
    return dict(sorted(hist.items(), key=lambda x: x[0]))


def precompute_atom_context(mol: Chem.Mol):
    dist_matrix = Chem.GetDistanceMatrix(mol)

    context = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        context[idx] = {
            "dist_matrix": dist_matrix,
        }

    return context


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
    sig = []
    for left, right in zip(path[:-1], path[1:]):
        bond = mol.GetBondBetweenAtoms(left, right)
        sig.append(bond_type_code(bond))
    return "-".join(sig)


def hop_shell_signatures(
    mol: Chem.Mol, atom_idx: int, max_hop: int = 2, dist_matrix=None
):
    max_hop = min(int(max_hop), 2)
    dist = dist_matrix if dist_matrix is not None else Chem.GetDistanceMatrix(mol)
    shells = {hop: [] for hop in range(1, max_hop + 1)}
    for idx in range(mol.GetNumAtoms()):
        hop = int(dist[atom_idx][idx])
        if hop < 1 or hop > max_hop:
            continue
        atom = mol.GetAtomWithIdx(idx)
        token = {
            "z": atom.GetAtomicNum(),
            "fc": atom.GetFormalCharge(),
            "ar": int(atom.GetIsAromatic()),
            "deg": atom.GetDegree(),
            "h": atom.GetTotalNumHs(includeNeighbors=True),
            "ring": int(atom.IsInRing()),
        }
        shells[hop].append(token)

    out = {}
    for hop in range(1, max_hop + 1):
        out[f"hop{hop}_shell"] = sorted(
            shells[hop],
            key=lambda x: (
                x["z"],
                x["fc"],
                x["ar"],
                x["deg"],
                x["h"],
                x["ring"],
            ),
        )
    return out


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


def make_env_key(
    mol: Chem.Mol, atom: Chem.Atom, hop_depth: int = 2, atom_ctx: dict | None = None
):
    atom_idx = atom.GetIdx()
    nbr_sigs = []
    bond_kinds = []

    for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        bcode = bond_type_code(bond)
        bond_kinds.append(bcode)
        nbr_sigs.append(f"{nbr.GetAtomicNum()}:{bcode}:{int(nbr.GetIsAromatic())}")

    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumAtomRings(atom_idx)
    ctx = atom_ctx or {}
    dist_matrix = ctx.get("dist_matrix")

    features = {
        "z": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "aromatic": int(atom.GetIsAromatic()),
        "hybridization": str(atom.GetHybridization()),
        "in_ring": int(atom.IsInRing()),
        "ring_count": ring_count,
        "degree": atom.GetDegree(),
        "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
        "neighbor_sig": sorted(nbr_sigs),
        "bond_kinds": sorted(bond_kinds),
    }
    features.update(
        hop_shell_signatures(mol, atom_idx, max_hop=hop_depth, dist_matrix=dist_matrix)
    )

    key = json.dumps(
        ordered_env_key_obj(features),
        ensure_ascii=False,
        sort_keys=False,
        separators=(",", ":"),
    )
    return key, features


def write_outputs(out_dir: Path, module: str, atom_rows):
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"{module}_atom_env.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "module",
                "atom_index",
                "atom_name",
                "opls_type_id",
                "opls_type_name",
                "charge",
                "sigma",
                "epsilon",
                "env_key",
            ],
        )
        writer.writeheader()
        for row in atom_rows:
            writer.writerow(row)

    return csv_path


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
        atom_context = precompute_atom_context(mol)

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

            env_key, _ = make_env_key(
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
        return write_outputs(out_dir, module, atom_rows)

    mol = load_structure(structure_path)
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
            mol_pdb = load_structure(pdb_fallback)
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


