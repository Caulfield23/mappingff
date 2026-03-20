#!/usr/bin/env python3
"""
脚本说明
--------
基于模块结构文件（.mol/.mol2）与对应的 LAMMPS data（.lammps.lmp），
提取每个原子的局部化学环境（env_key），并建立 env_key -> lmp_type 的映射库。

输出文件
--------
- {module}_atom_env.csv
- {module}_env_to_opls.json
- {module}_env_to_opls.sqlite
- {module}_keymap.csv

用法示例
--------
python scripts/build_envkey_mapping.py \
    --mol /path/to/segment1.mol \
    --lmp /path/to/segment1.lammps.lmp \
    --module segment1 \
    --outdir /path/to/outputs/segment1_envdb
"""

import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


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
            f"无法根据质量 {mass} 推断元素（最小差值 {best[0] if best else 'N/A'}）"
        )
    return best[1]


def parse_lammps_data(lmp_path: Path):
    lines = lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    section_names = {
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

    current = None
    masses = {}
    pair_coeffs = {}
    atoms = {}

    def _parse_atom_line(toks):
        if len(toks) < 6:
            return None

        atom_id = int(toks[0])

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

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped in section_names:
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
            epsilon = float(toks[1])
            sigma = float(toks[2])
            pair_coeffs[type_id] = {"epsilon": epsilon, "sigma": sigma}

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

    if not masses:
        raise ValueError(f"在 {lmp_path} 中未解析到 Masses")
    if not pair_coeffs:
        raise ValueError(f"在 {lmp_path} 中未解析到 Pair Coeffs")
    if not atoms:
        raise ValueError(f"在 {lmp_path} 中未解析到 Atoms")

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
            raise ValueError(f"LAMMPS atom type {lmp_type} 缺少 Pair Coeffs")
        if lmp_type not in type_info:
            raise ValueError(f"LAMMPS atom type {lmp_type} 缺少 Masses")
        info["sigma"] = pair_coeffs[lmp_type]["sigma"]
        info["epsilon"] = pair_coeffs[lmp_type]["epsilon"]
        info["atomic_num"] = type_info[lmp_type]["atomic_num"]
        info["type_name"] = type_info[lmp_type]["type_name"]

    return atoms, type_info


def load_structure(structure_path: Path):
    suffix = structure_path.suffix.lower()

    def _sanitize_or_raise(mol, source_path: Path, label: str):
        if mol is None:
            raise ValueError(f"RDKit 无法读取 {label} 文件: {source_path}")
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            raise ValueError(f"RDKit 读取到 {label} 但清洗失败: {source_path} ({exc})")
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
                f"[INFO] mol 读取失败，已自动回退为 pdb 输入: {pdb_fallback}",
                flush=True,
            )
            return mol

        raise ValueError(
            f"RDKit 无法读取 mol 文件: {structure_path}；且未能从同名 pdb 回退: {pdb_fallback}"
        )

    if suffix == ".pdb":
        mol = _load_pdb(structure_path)
        if mol is None:
            raise ValueError(f"RDKit 无法读取 pdb 文件: {structure_path}")
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
        f"不支持的结构文件格式: {structure_path}. 目前支持 .mol/.mol2/.pdb"
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


def third_shell_element_counts(mol: Chem.Mol, atom_idx: int, dist_matrix=None):
    dist = dist_matrix if dist_matrix is not None else Chem.GetDistanceMatrix(mol)
    counts = Counter()
    for idx in range(mol.GetNumAtoms()):
        if int(dist[atom_idx][idx]) == 3:
            counts[mol.GetAtomWithIdx(idx).GetAtomicNum()] += 1
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
    mol: Chem.Mol, atom_idx: int, max_hop: int = 3, dist_matrix=None
):
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
    mol: Chem.Mol, atom: Chem.Atom, hop_depth: int = 3, atom_ctx: dict | None = None
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
    key_hash = hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]
    return key, key_hash, features


def write_outputs(out_dir: Path, module: str, atom_rows, env_map_counter):
    out_dir.mkdir(parents=True, exist_ok=True)

    key_meta = {}
    for env_key, counter in env_map_counter.items():
        candidates = sorted(counter.items(), key=lambda x: (-x[1], x[0][0]))
        n_candidates_raw = len(candidates)
        name_counts = Counter()
        for (type_id, type_name), count in counter.items():
            name_counts[type_name] += int(count)
        canonical_candidates = sorted(name_counts.items(), key=lambda x: (-x[1], x[0]))
        n_candidates_name = len(canonical_candidates)
        primary_raw = candidates[0][0] if candidates else None
        primary_name = canonical_candidates[0][0] if canonical_candidates else None
        key_meta[env_key] = {
            "n_candidates": n_candidates_name,
            "n_candidates_raw": n_candidates_raw,
            "n_candidates_name": n_candidates_name,
            "primary_raw": primary_raw,
            "primary_name": primary_name,
        }

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
                "env_key_hash",
                "n_candidates",
                "n_candidates_raw",
                "n_candidates_name",
                "is_primary_name",
                "is_primary_raw",
                "env_key",
            ],
        )
        writer.writeheader()
        for row in atom_rows:
            meta = key_meta[row["env_key"]]
            row["n_candidates"] = meta["n_candidates"]
            row["n_candidates_raw"] = meta["n_candidates_raw"]
            row["n_candidates_name"] = meta["n_candidates_name"]
            row["is_primary_name"] = int(row["opls_type_name"] == meta["primary_name"])
            row["is_primary_raw"] = int(
                meta["primary_raw"] is not None
                and (row["opls_type_id"], row["opls_type_name"]) == meta["primary_raw"]
            )
            writer.writerow(row)

    json_path = out_dir / f"{module}_env_to_opls.json"
    keymap_csv_path = out_dir / f"{module}_keymap.csv"
    env_items = []
    keymap_rows = []
    for env_key, counter in env_map_counter.items():
        candidates = []
        for (type_id, type_name), count in sorted(
            counter.items(), key=lambda x: (-x[1], x[0][0])
        ):
            candidates.append(
                {
                    "opls_type_id": type_id,
                    "opls_type_name": type_name,
                    "count": count,
                }
            )
        key_hash = hashlib.sha1(env_key.encode("utf-8")).hexdigest()[:16]
        name_counts = Counter()
        for (type_id, type_name), count in counter.items():
            name_counts[type_name] += int(count)

        canonical_candidates = [
            {"opls_type_name": name, "count": cnt}
            for name, cnt in sorted(name_counts.items(), key=lambda x: (-x[1], x[0]))
        ]

        n_candidates_raw = len(candidates)
        n_candidates_name = len(canonical_candidates)
        primary_raw = candidates[0] if candidates else None
        primary_name = canonical_candidates[0] if canonical_candidates else None

        env_items.append(
            {
                "env_key_hash": key_hash,
                "env_key": env_key,
                "n_candidates": n_candidates_name,
                "n_candidates_raw": n_candidates_raw,
                "n_candidates_name": n_candidates_name,
                "candidates": candidates,
                "canonical_candidates": canonical_candidates,
                "primary": primary_name,
                "primary_raw": primary_raw,
            }
        )

        for cand in candidates:
            keymap_rows.append(
                {
                    "env_key_hash": key_hash,
                    "opls_type_id": cand["opls_type_id"],
                    "opls_type_name": cand["opls_type_name"],
                    "count": cand["count"],
                    "n_candidates": n_candidates_name,
                    "n_candidates_raw": n_candidates_raw,
                    "n_candidates_name": n_candidates_name,
                    "is_primary_name": int(
                        primary_name is not None
                        and cand["opls_type_name"] == primary_name["opls_type_name"]
                    ),
                    "is_primary_raw": int(
                        primary_raw is not None
                        and cand["opls_type_id"] == primary_raw["opls_type_id"]
                        and cand["opls_type_name"] == primary_raw["opls_type_name"]
                    ),
                    "env_key": env_key,
                }
            )

    json_data = {
        "module": module,
        "n_atoms": len(atom_rows),
        "n_unique_env_keys": len(env_items),
        "env_to_opls": sorted(env_items, key=lambda x: x["env_key_hash"]),
    }
    json_path.write_text(
        json.dumps(json_data, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    db_path = out_dir / f"{module}_env_to_opls.sqlite"
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS atom_env")
    cur.execute("DROP TABLE IF EXISTS env_to_opls")
    cur.execute(
        """
        CREATE TABLE atom_env (
            module TEXT,
            atom_index INTEGER,
            atom_name TEXT,
            opls_type_id INTEGER,
            opls_type_name TEXT,
            charge REAL,
            sigma REAL,
            epsilon REAL,
            env_key_hash TEXT,
            n_candidates INTEGER,
            n_candidates_raw INTEGER,
            n_candidates_name INTEGER,
            is_primary_name INTEGER,
            is_primary_raw INTEGER,
            env_key TEXT
        )
        """
    )
    cur.execute(
        """
        CREATE TABLE env_to_opls (
            env_key_hash TEXT,
            env_key TEXT,
            opls_type_id INTEGER,
            opls_type_name TEXT,
            count INTEGER,
            n_candidates INTEGER,
            n_candidates_raw INTEGER,
            n_candidates_name INTEGER,
            is_primary_name INTEGER,
            is_primary_raw INTEGER,
            PRIMARY KEY (env_key, opls_type_id, opls_type_name)
        )
        """
    )

    cur.executemany(
        """
        INSERT INTO atom_env (
            module, atom_index, atom_name, opls_type_id, opls_type_name,
            charge, sigma, epsilon, env_key_hash,
            n_candidates, n_candidates_raw, n_candidates_name, is_primary_name, is_primary_raw,
            env_key
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                r["module"],
                r["atom_index"],
                r["atom_name"],
                r["opls_type_id"],
                r["opls_type_name"],
                r["charge"],
                r["sigma"],
                r["epsilon"],
                r["env_key_hash"],
                r["n_candidates"],
                r["n_candidates_raw"],
                r["n_candidates_name"],
                r["is_primary_name"],
                r["is_primary_raw"],
                r["env_key"],
            )
            for r in atom_rows
        ],
    )

    env_rows = []
    for row in keymap_rows:
        env_rows.append(
            (
                row["env_key_hash"],
                row["env_key"],
                row["opls_type_id"],
                row["opls_type_name"],
                row["count"],
                row["n_candidates"],
                row["n_candidates_raw"],
                row["n_candidates_name"],
                row["is_primary_name"],
                row["is_primary_raw"],
            )
        )

    cur.executemany(
        """
        INSERT INTO env_to_opls (
            env_key_hash, env_key, opls_type_id, opls_type_name, count,
            n_candidates, n_candidates_raw, n_candidates_name, is_primary_name, is_primary_raw
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        env_rows,
    )

    with keymap_csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "env_key_hash",
                "opls_type_id",
                "opls_type_name",
                "count",
                "n_candidates",
                "n_candidates_raw",
                "n_candidates_name",
                "is_primary_name",
                "is_primary_raw",
                "env_key",
            ],
        )
        writer.writeheader()
        for row in keymap_rows:
            writer.writerow(row)

    conn.commit()
    conn.close()

    return csv_path, json_path, db_path, keymap_csv_path


def build_mapping(
    structure_path: Path,
    out_dir: Path,
    module: str,
    lmp_path: Path = None,
    hop_depth: int = 3,
):
    if lmp_path is None:
        raise ValueError("纯 LAMMPS 路线必须提供 --lmp")

    lmp_atoms, _ = parse_lammps_data(lmp_path)

    def _build_with_loaded_mol(mol: Chem.Mol, source_path: Path):
        if mol.GetNumAtoms() != len(lmp_atoms):
            raise ValueError(
                f"原子数不一致: structure={mol.GetNumAtoms()} vs lammps_atoms={len(lmp_atoms)}。"
                f" structure={source_path}"
            )

        atom_rows = []
        env_map_counter = defaultdict(Counter)
        atom_context = precompute_atom_context(mol)

        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom_id = atom_idx + 1
            lmp_atom = lmp_atoms.get(atom_id)
            if lmp_atom is None:
                raise ValueError(f"LAMMPS Atoms 中缺少 atom_id={atom_id}")

            type_id = lmp_atom["lmp_type"]
            if lmp_atom["atomic_num"] != atom.GetAtomicNum():
                raise ValueError(
                    f"原子元素不一致: atom_index={atom_id}, structure={source_path}, "
                    f"structure Z={atom.GetAtomicNum()} vs lammps_mass Z={lmp_atom['atomic_num']}"
                )

            charge = lmp_atom["charge"]
            sigma = lmp_atom["sigma"]
            epsilon = lmp_atom["epsilon"]

            env_key, env_hash, _ = make_env_key(
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
                "env_key_hash": env_hash,
                "env_key": env_key,
            }
            atom_rows.append(row)
            env_map_counter[env_key][(type_id, lmp_atom["type_name"])] += 1

        return write_outputs(out_dir, module, atom_rows, env_map_counter)

    mol = load_structure(structure_path)
    try:
        return _build_with_loaded_mol(mol, structure_path)
    except Exception as primary_exc:
        if structure_path.suffix.lower() != ".mol":
            raise

        pdb_fallback = structure_path.with_suffix(".pdb")
        if not pdb_fallback.exists():
            raise ValueError(
                f"使用 mol 建库失败，且未找到可回退 pdb: {pdb_fallback}\n"
                f"原始错误: {primary_exc}"
            )

        try:
            mol_pdb = load_structure(pdb_fallback)
            result = _build_with_loaded_mol(mol_pdb, pdb_fallback)
            print(
                f"[INFO] 使用 mol 建库失败，已自动切换 pdb: {pdb_fallback}",
                flush=True,
            )
            return result
        except Exception as pdb_exc:
            raise ValueError(
                "mol 与 pdb 建库均失败。\n"
                f"- mol: {structure_path}\n"
                f"  错误: {primary_exc}\n"
                f"- pdb: {pdb_fallback}\n"
                f"  错误: {pdb_exc}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="为模块分子生成 env_key，并基于 LAMMPS data 文件建立 env_key -> lmp_type 映射数据库。",
        epilog=(
            "示例:\n"
            "  python scripts/build_envkey_mapping.py "
            "--mol /path/to/segment1.mol "
            "--lmp /path/to/segment1.lammps.lmp "
            "--module segment1 --outdir /path/to/outputs/segment1_envdb"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--mol",
        required=True,
        type=Path,
        help="模块结构文件路径（支持 .mol/.mol2/.pdb；若 .mol 失败会自动尝试同名 .pdb）",
    )
    parser.add_argument(
        "--lmp",
        required=True,
        type=Path,
        help="同模块的 .lammps.lmp 文件（纯 LAMMPS 路线）",
    )
    parser.add_argument(
        "--hop-depth",
        required=False,
        type=int,
        default=3,
        help="环境 hop 深度，默认 3；更高更精细但键空间更稀疏。",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="输出目录")
    parser.add_argument("--module", required=True, help="模块名，例如 segment1")

    args = parser.parse_args()

    csv_path, json_path, db_path, keymap_csv_path = build_mapping(
        structure_path=args.mol,
        out_dir=args.outdir,
        module=args.module,
        lmp_path=args.lmp,
        hop_depth=args.hop_depth,
    )

    print("完成：")
    print(f"- atom env 表: {csv_path}")
    print(f"- env->opls JSON: {json_path}")
    print(f"- env->opls SQLite: {db_path}")
    print(f"- module keymap CSV: {keymap_csv_path}")


if __name__ == "__main__":
    main()
