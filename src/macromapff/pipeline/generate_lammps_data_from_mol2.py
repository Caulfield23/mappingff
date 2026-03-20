#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from itertools import combinations
from itertools import product
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

from build_envkey_mapping import make_env_key, precompute_atom_context


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

SHELL_ITEM_KEY_PRIORITY = ["z", "fc", "ar", "deg", "h", "ring"]


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


def _stable_sort_list(values):
    return sorted(
        values,
        key=lambda x: json.dumps(
            x, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ),
    )


def _normalize_env_obj(obj, parent_key: str = ""):
    if isinstance(obj, dict):
        if parent_key == "":
            ordered = ordered_env_key_obj(obj)
            out = {}
            for key, val in ordered.items():
                out[key] = _normalize_env_obj(val, parent_key=key)
            return out

        if all(k in SHELL_ITEM_KEY_PRIORITY for k in obj.keys()):
            out = {}
            for key in SHELL_ITEM_KEY_PRIORITY:
                if key in obj:
                    out[key] = _normalize_env_obj(obj[key], parent_key=key)
            for key in sorted(k for k in obj.keys() if k not in out):
                out[key] = _normalize_env_obj(obj[key], parent_key=key)
            return out

        out = {}
        for key in sorted(obj.keys()):
            out[key] = _normalize_env_obj(obj[key], parent_key=key)
        return out

    if isinstance(obj, list):
        norm = [_normalize_env_obj(v, parent_key=parent_key) for v in obj]
        if parent_key in {"neighbor_sig", "bond_kinds"} or (
            parent_key.startswith("hop") and parent_key.endswith("_shell")
        ):
            return _stable_sort_list(norm)
        return norm

    return obj


def load_input_structure(structure_path: Path) -> Chem.Mol:
    suffix = structure_path.suffix.lower()
    if suffix != ".mol":
        raise ValueError(f"当前仅支持 .mol 输入，收到: {structure_path}")

    mol = Chem.MolFromMolFile(
        str(structure_path), removeHs=False, sanitize=True, strictParsing=False
    )
    if mol is None:
        mol = Chem.MolFromMolFile(
            str(structure_path), removeHs=False, sanitize=False, strictParsing=False
        )
        if mol is None:
            raise ValueError(f"RDKit 无法读取 .mol 文件: {structure_path}")
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ),
            )
    return mol


def canonicalize_env_key(env_key_raw: str) -> str:
    try:
        obj = json.loads(env_key_raw)
        obj = _normalize_env_obj(obj)
        return json.dumps(
            obj,
            ensure_ascii=False,
            sort_keys=False,
            separators=(",", ":"),
        )
    except Exception:
        return (env_key_raw or "").strip()


def load_hop_param_db(hop_csv: Path):
    if not hop_csv.exists():
        raise FileNotFoundError(f"找不到 hop 参数数据库: {hop_csv}")

    env_to_atom_param = {}
    with hop_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            env_key = canonicalize_env_key(row["env_key"])
            source_ids = [x for x in str(row.get("source_key_ids", "")).split(";") if x]
            key_ids = sorted({int(x) for x in source_ids}) if source_ids else [-1]
            env_to_atom_param[env_key] = {
                "key_ids": key_ids,
                "charge": float(row["charge_mean"]),
                "sigma": float(row["sigma_mean"]),
                "epsilon": float(row["epsilon_mean"]),
                "mass": float(row["mass_mean"]),
                "match_count": int(float(row.get("match_count", "1") or 1)),
            }

    if not env_to_atom_param:
        raise ValueError(f"hop 参数数据库为空: {hop_csv}")
    return env_to_atom_param


def _parse_json(raw: str):
    try:
        return json.loads(raw)
    except Exception:
        return None


def _norm_rev_tuple(values):
    vals = tuple(values)
    rev = tuple(reversed(vals))
    return vals if vals <= rev else rev


def _parse_fallback_hops(raw: str):
    if raw is None:
        return tuple()
    items = []
    for token in str(raw).split(","):
        token = token.strip()
        if not token:
            continue
        hop = int(token)
        if hop < 0:
            raise ValueError(f"fallback hop 必须 >= 0，收到: {hop}")
        if hop not in items:
            items.append(hop)
    return tuple(items)


def _env_key_at_hop_or_none(env_key_raw: str, hop: int):
    obj = _parse_json(env_key_raw)
    if not isinstance(obj, dict):
        return None

    reduced = dict(obj)
    for key in list(reduced.keys()):
        if key.startswith("hop") and key.endswith("_shell"):
            try:
                num = int(key[3:-6])
            except Exception:
                continue
            if num > hop:
                reduced.pop(key, None)

    return canonicalize_env_key(json.dumps(reduced, ensure_ascii=False))


def write_missing_env_log(
    log_path: Path,
    structure_path: Path,
    missing,
    hop_depth: int,
    fallback_hops,
):
    lines = [
        "[ATOM][MISSING] report",
        f"[ATOM][MISSING] structure={structure_path}",
        f"[ATOM][MISSING] hop_depth={hop_depth}",
        f"[ATOM][MISSING] fallback_hops={','.join(str(x) for x in fallback_hops)}",
        f"[ATOM][MISSING] missing_count={len(missing)}",
    ]
    for item in missing:
        lines.append(
            f"[ATOM][MISSING] atom_index={item['atom_index']} symbol={item['symbol']}"
        )
        lines.append(f"[ATOM][MISSING] env_key_exact={item['env_key_exact']}")
        if item.get("fallback_candidates"):
            for hop, key in item["fallback_candidates"]:
                lines.append(f"[ATOM][MISSING] env_key_hop{hop}={key}")
    append_build_log(log_path, lines)


def append_build_log(log_path: Path, lines):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


def write_atom_keytype_map(out_path: Path, atom_records):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["atom_index", "chosen_key_type", "all_key_types"])
        for rec in atom_records:
            all_key_types = sorted({int(x) for x in rec.get("global_key_ids", [])})
            writer.writerow(
                [
                    int(rec["atom_id"]),
                    int(rec.get("global_key_id", 0)),
                    ";".join(str(x) for x in all_key_types),
                ]
            )


def init_build_log(log_path: Path, lines):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


def load_multiatom_db(multiatom_csv: Path):
    if not multiatom_csv.exists():
        raise FileNotFoundError(f"找不到多原子参数数据库: {multiatom_csv}")

    reversible_kinds = {"bond", "angle", "dihedral"}
    idx_rev_patterns = defaultdict(list)
    idx_imp_patterns = []

    def _to_allowed_tuple(raw_list):
        allowed = []
        for item in raw_list:
            if isinstance(item, list):
                allowed.append(frozenset(int(x) for x in item))
            else:
                allowed.append(frozenset([int(item)]))
        return tuple(allowed)

    with multiatom_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            kind = row["interaction_kind"].strip().lower()
            key_tuple_raw = _parse_json(row["key_type_tuple"])
            coeff_sets_raw = _parse_json(row["coeff_param_sets"])
            if not isinstance(key_tuple_raw, list) or not isinstance(
                coeff_sets_raw, list
            ):
                continue
            allowed = _to_allowed_tuple(key_tuple_raw)
            coeff_sets = {
                tuple(str(v) for v in one)
                for one in coeff_sets_raw
                if isinstance(one, list)
            }
            if not coeff_sets:
                continue

            if kind in reversible_kinds:
                idx_rev_patterns[kind].append(
                    {
                        "allowed": allowed,
                        "coeff_candidates": set(coeff_sets),
                    }
                )
            elif kind == "improper":
                if len(allowed) != 4:
                    continue
                idx_imp_patterns.append(
                    {
                        "center_allowed": allowed[0],
                        "other_allowed": (allowed[1], allowed[2], allowed[3]),
                        "coeff_candidates": set(coeff_sets),
                    }
                )

    return idx_rev_patterns, idx_imp_patterns


def _tuple_matches_allowed(key_tuple, allowed):
    if len(key_tuple) != len(allowed):
        return False
    return all(key_tuple[i] in allowed[i] for i in range(len(key_tuple)))


def _improper_matches_pattern(key_tuple, pattern):
    if len(key_tuple) != 4:
        return False
    if key_tuple[0] not in pattern["center_allowed"]:
        return False

    vals = [key_tuple[1], key_tuple[2], key_tuple[3]]
    slots = pattern["other_allowed"]

    for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        if all(vals[perm[i]] in slots[i] for i in range(3)):
            return True
    return False


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
    atom_context = precompute_atom_context(mol)

    atom_records = []
    missing = []
    fallback_hit_counter = defaultdict(int)

    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        env_key_raw, _, _ = make_env_key(
            mol,
            atom,
            hop_depth=hop_depth,
            atom_ctx=atom_context.get(atom_idx),
        )
        env_key = canonicalize_env_key(env_key_raw)
        key_at_hop2 = _env_key_at_hop_or_none(env_key, 2)
        params = hop2_env_to_atom_param.get(key_at_hop2) if key_at_hop2 else None
        fallback_candidates = []

        matched_hop = 2
        if params is None:
            for fb_hop in fallback_hops:
                if fb_hop == 1:
                    fb_key = _env_key_at_hop_or_none(env_key, 1)
                    if fb_key:
                        fallback_candidates.append((1, fb_key))
                        params = hop1_env_to_atom_param.get(fb_key)
                elif fb_hop == 0:
                    fb_key = _env_key_at_hop_or_none(env_key, 0)
                    if fb_key:
                        fallback_candidates.append((0, fb_key))
                        params = hop0_env_to_atom_param.get(fb_key)
                else:
                    continue

                if params is not None:
                    matched_hop = fb_hop
                    fallback_hit_counter[fb_hop] += 1
                    break

        if params is None:
            missing.append(
                {
                    "atom_index": atom_idx + 1,
                    "symbol": atom.GetSymbol(),
                    "env_key_exact": env_key,
                    "fallback_candidates": fallback_candidates,
                }
            )
            atom_records.append(
                {
                    "atom_id": atom_idx + 1,
                    "global_key_id": 0,
                    "global_key_ids": [0],
                    "charge": 0.0,
                    "sigma": 0.0,
                    "epsilon": 0.0,
                    "mass": 0.0,
                    "matched_hop_depth": -1,
                }
            )
            continue

        atom_records.append(
            {
                "atom_id": atom_idx + 1,
                "global_key_id": params["key_ids"][0],
                "global_key_ids": sorted(set(params["key_ids"])),
                "charge": params["charge"],
                "sigma": params["sigma"],
                "epsilon": params["epsilon"],
                "mass": params["mass"],
                "matched_hop_depth": matched_hop,
            }
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
            "[ATOM] build_atom_types 完成",
            f"[ATOM] total_atoms={len(atom_records)}",
            f"[ATOM] missing_atoms={len(missing)}",
            "[ATOM] fallback_hits="
            + ", ".join(f"hop{h}:{c}" for h, c in sorted(fallback_hit_counter.items())),
        ],
    )

    used_global = sorted({a["global_key_id"] for a in atom_records})
    g2l = {gid: i + 1 for i, gid in enumerate(used_global)}

    for a in atom_records:
        a["atom_type"] = g2l[a["global_key_id"]]

    local_type_params = []
    for gid in used_global:
        any_atom = next(a for a in atom_records if a["global_key_id"] == gid)
        local_type_params.append(
            {
                "local_type": g2l[gid],
                "global_key_id": gid,
                "mass": any_atom["mass"],
                "sigma": any_atom["sigma"],
                "epsilon": any_atom["epsilon"],
            }
        )
    local_type_params.sort(key=lambda x: x["local_type"])

    return atom_records, local_type_params, dict(sorted(fallback_hit_counter.items()))


def enumerate_terms(mol: Chem.Mol):
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx() + 1
        j = bond.GetEndAtomIdx() + 1
        bonds.append((i, j))

    angles = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        for i, k in combinations(sorted(nbs), 2):
            angles.append((i, center + 1, k))

    dihedral_set = set()
    for bond in mol.GetBonds():
        j = bond.GetBeginAtomIdx() + 1
        k = bond.GetEndAtomIdx() + 1
        j_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(j - 1).GetNeighbors()
            if x.GetIdx() + 1 != k
        ]
        k_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(k - 1).GetNeighbors()
            if x.GetIdx() + 1 != j
        ]
        for i in j_neighbors:
            for l in k_neighbors:
                t = (i, j, k, l)
                dihedral_set.add(_norm_rev_tuple(t))
    dihedrals = sorted(dihedral_set)

    impropers = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        if len(nbs) < 3:
            continue
        for a, b, c in combinations(sorted(nbs), 3):
            impropers.append((center + 1, a, b, c))

    return bonds, angles, dihedrals, impropers


def assign_multiatom_params(
    kind: str,
    terms,
    atom_key_sets,
    idx_rev_patterns,
    idx_imp_patterns,
    strict_missing=True,
    build_log_path: Path = None,
):
    records = []
    missing = []
    ambiguous = []

    for atom_ids in terms:
        key_options = [sorted(atom_key_sets[x]) for x in atom_ids]
        if any(len(options) == 0 for options in key_options):
            missing.append((atom_ids, tuple()))
            continue

        matched_tuples = {}
        for key_tuple in product(*key_options):
            if kind in {"bond", "angle", "dihedral"}:
                rev = tuple(reversed(key_tuple))
                coeff_candidates = set()
                for pattern in idx_rev_patterns[kind]:
                    if _tuple_matches_allowed(
                        key_tuple, pattern["allowed"]
                    ) or _tuple_matches_allowed(rev, pattern["allowed"]):
                        coeff_candidates.update(pattern["coeff_candidates"])
            else:
                coeff_candidates = set()
                for pattern in idx_imp_patterns:
                    if _improper_matches_pattern(key_tuple, pattern):
                        coeff_candidates.update(pattern["coeff_candidates"])
            if coeff_candidates:
                matched_tuples[key_tuple] = sorted(coeff_candidates)

        if not matched_tuples:
            combo_sample = []
            for i, combo in enumerate(product(*key_options)):
                combo_sample.append(combo)
                if i >= 19:
                    break
            missing.append((atom_ids, tuple(combo_sample)))
            records.append(
                {
                    "atom_ids": atom_ids,
                    "key_type_tuple": (
                        tuple(combo_sample[0]) if combo_sample else tuple()
                    ),
                    "coeff": None,
                }
            )
            continue

        sorted_matches = sorted(matched_tuples.items(), key=lambda x: x[0])
        chosen_key_tuple, chosen_coeff_list = sorted_matches[0]
        coeff = chosen_coeff_list[0]

        unique_coeffs = sorted(
            {tuple(c) for values in matched_tuples.values() for c in values}
        )
        if len(unique_coeffs) > 1:
            ambiguous.append(
                {
                    "atom_ids": atom_ids,
                    "matched_key_tuples": [k for k, _ in sorted_matches],
                    "matched_coeffs": unique_coeffs,
                    "chosen_key_tuple": chosen_key_tuple,
                    "chosen_coeff": coeff,
                }
            )

        records.append(
            {"atom_ids": atom_ids, "key_type_tuple": chosen_key_tuple, "coeff": coeff}
        )

    if missing:
        sample = ", ".join(f"atom_ids={x[0]} candidates={x[1]}" for x in missing[:10])
        append_build_log(
            build_log_path,
            [
                (
                    f"[{kind.upper()}][ERROR] missing_count={len(missing)}"
                    if strict_missing
                    else f"[{kind.upper()}][WARN] missing_count={len(missing)}"
                ),
                (
                    f"[{kind.upper()}][ERROR] sample_missing={sample}"
                    if strict_missing
                    else f"[{kind.upper()}][WARN] sample_missing={sample}"
                ),
            ],
        )
        if strict_missing:
            raise ValueError(
                f"{kind} 有 {len(missing)} 条连接在参数表中找不到。示例 key_type_tuple: {sample}"
            )

    if ambiguous:
        lines = [f"[{kind.upper()}][AMBIGUOUS] count={len(ambiguous)}"]
        for item in ambiguous[:200]:
            lines.append(
                f"[{kind.upper()}][AMBIGUOUS] atom_ids={item['atom_ids']} "
                f"chosen_key_tuple={item['chosen_key_tuple']} chosen_coeff={item['chosen_coeff']} "
                f"all_key_tuples={item['matched_key_tuples']} all_coeffs={item['matched_coeffs']}"
            )
        append_build_log(build_log_path, lines)

    return records, missing, ambiguous


def build_type_map(records):
    coeff_to_type = {}
    for rec in records:
        coeff = rec["coeff"]
        if coeff is None:
            continue
        if coeff not in coeff_to_type:
            coeff_to_type[coeff] = len(coeff_to_type) + 1

    for rec in records:
        if rec["coeff"] is None:
            rec["type_id"] = 0
        else:
            rec["type_id"] = coeff_to_type[rec["coeff"]]

    type_rows = sorted(
        [(tid, coeff) for coeff, tid in coeff_to_type.items()], key=lambda x: x[0]
    )
    return type_rows


def _fmt_box(lo, hi):
    return f"{lo:12.6f} {hi:12.6f}"


def write_lammps_data(
    out_path: Path,
    mol: Chem.Mol,
    atom_records,
    atom_type_rows,
    bond_records,
    angle_records,
    dihedral_records,
    improper_records,
    bond_type_rows,
    angle_type_rows,
    dihedral_type_rows,
    improper_type_rows,
    molecule_id: int,
    box_padding: float,
):
    conf = mol.GetConformer()
    xs, ys, zs = [], [], []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(float(p.x))
        ys.append(float(p.y))
        zs.append(float(p.z))

    xlo, xhi = min(xs) - box_padding, max(xs) + box_padding
    ylo, yhi = min(ys) - box_padding, max(ys) + box_padding
    zlo, zhi = min(zs) - box_padding, max(zs) + box_padding

    atom_by_id = {a["atom_id"]: a for a in atom_records}

    improper_has_zero = any(int(r.get("type_id", 0)) == 0 for r in improper_records)
    if improper_has_zero:
        improper_records_out = []
        for rec in improper_records:
            rec_new = dict(rec)
            rec_new["type_id"] = max(1, int(rec_new.get("type_id", 0)))
            improper_records_out.append(rec_new)
        improper_type_rows_out = improper_type_rows
    else:
        improper_records_out = improper_records
        improper_type_rows_out = improper_type_rows

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("LAMMPS data file generated by generate_lammps_data_from_mol2.py\n\n")
        f.write(f"{len(atom_records):8d} atoms\n")
        f.write(f"{len(bond_records):8d} bonds\n")
        f.write(f"{len(angle_records):8d} angles\n")
        f.write(f"{len(dihedral_records):8d} dihedrals\n")
        f.write(f"{len(improper_records_out):8d} impropers\n\n")

        f.write(f"{len(atom_type_rows):8d} atom types\n")
        f.write(f"{len(bond_type_rows):8d} bond types\n")
        f.write(f"{len(angle_type_rows):8d} angle types\n")
        f.write(f"{len(dihedral_type_rows):8d} dihedral types\n")
        f.write(f"{len(improper_type_rows_out):8d} improper types\n\n")

        f.write(f"{_fmt_box(xlo, xhi)} xlo xhi\n")
        f.write(f"{_fmt_box(ylo, yhi)} ylo yhi\n")
        f.write(f"{_fmt_box(zlo, zhi)} zlo zhi\n\n")

        f.write("Masses\n\n")
        for row in atom_type_rows:
            f.write(f"{row['local_type']:8d} {row['mass']:10.6f}\n")
        f.write("\n")

        f.write("Pair Coeffs\n\n")
        for row in atom_type_rows:
            f.write(
                f"{row['local_type']:8d} {row['epsilon']:10.6f} {row['sigma']:10.6f}\n"
            )
        f.write("\n")

        f.write("Bond Coeffs\n\n")
        for tid, coeff in bond_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.4f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Angle Coeffs\n\n")
        for tid, coeff in angle_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.3f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Dihedral Coeffs\n\n")
        for tid, coeff in dihedral_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.3f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Improper Coeffs\n\n")
        for tid, coeff in improper_type_rows_out:
            f.write(
                f"{tid:8d} {float(coeff[0]):10.3f} {int(float(coeff[1])):7d} {int(float(coeff[2])):7d}\n"
            )
        f.write("\n")

        f.write("Atoms\n\n")
        for atom_id in range(1, mol.GetNumAtoms() + 1):
            p = conf.GetAtomPosition(atom_id - 1)
            rec = atom_by_id[atom_id]
            f.write(
                f"{atom_id:8d} {molecule_id:6d} {rec['atom_type']:6d} "
                f"{rec['charge']:10.6f} {float(p.x):10.5f} {float(p.y):10.5f} {float(p.z):10.5f}\n"
            )
        f.write("\n")

        f.write("Bonds\n\n")
        for idx, rec in enumerate(bond_records, start=1):
            i, j = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d}\n")
        f.write("\n")

        f.write("Angles\n\n")
        for idx, rec in enumerate(angle_records, start=1):
            i, j, k = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d}\n")
        f.write("\n")

        f.write("Dihedrals\n\n")
        for idx, rec in enumerate(dihedral_records, start=1):
            i, j, k, l = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d} {l:6d}\n")
        f.write("\n")

        f.write("Impropers\n\n")
        for idx, rec in enumerate(improper_records_out, start=1):
            i, j, k, l = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d} {l:6d}\n")


def main():
    RDLogger.DisableLog("rdApp.warning")

    parser = argparse.ArgumentParser(
        description=(
            "读取完整分子（仅 .mol），按 final_env_keymap + multiatom_master 数据库匹配参数，"
            "写出完整 LAMMPS data 文件。"
        )
    )
    parser.add_argument(
        "--structure", required=True, type=Path, help="输入完整分子 .mol 文件"
    )
    parser.add_argument(
        "--final-env-csv",
        type=Path,
        default=Path("outputs/final_env_keymap.csv"),
        help="原子参数数据库 CSV",
    )
    parser.add_argument(
        "--hop2-env-csv",
        type=Path,
        default=Path("outputs/hop2_env_keymap.csv"),
        help="hop2 原子参数数据库 CSV",
    )
    parser.add_argument(
        "--hop1-env-csv",
        type=Path,
        default=Path("outputs/hop1_env_keymap.csv"),
        help="hop1 原子参数数据库 CSV",
    )
    parser.add_argument(
        "--hop0-env-csv",
        type=Path,
        default=Path("outputs/hop0_env_keymap.csv"),
        help="hop0 兜底原子参数数据库 CSV",
    )
    parser.add_argument(
        "--multiatom-csv",
        type=Path,
        default=Path("outputs/multiatom_master_keytype.csv"),
        help="多原子参数数据库 CSV",
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="输出 LAMMPS data 文件路径"
    )
    parser.add_argument(
        "--hop-depth", type=int, default=4, help="环境特征 hop 深度，默认 4"
    )
    parser.add_argument(
        "--fallback-hops",
        type=str,
        default="1,0",
        help="env_key 未命中时按顺序降级匹配的 hop 深度列表（逗号分隔），默认 1,0",
    )
    parser.add_argument(
        "--box-padding", type=float, default=20.0, help="边界盒 padding，默认 20.0"
    )
    parser.add_argument(
        "--molecule-id", type=int, default=1, help="Atoms 段中的 molecule id，默认 1"
    )
    parser.add_argument(
        "--build-log",
        type=Path,
        default=None,
        help="构建过程日志输出路径（默认: <out_dir>/build.log）",
    )
    parser.add_argument(
        "--atom-keytype-map",
        type=Path,
        default=None,
        help="输出每个原子索引对应 key_type 的映射文件（默认: <out_dir>/atom_index_key_types.csv）",
    )
    args = parser.parse_args()

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
            f"hop_depth: {args.hop_depth}",
            f"fallback_hops: {args.fallback_hops}",
        ],
    )

    try:
        mol = load_input_structure(args.structure)
        if mol is None:
            raise ValueError(f"RDKit 无法读取结构文件: {args.structure}")
        if mol.GetNumConformers() == 0:
            raise ValueError("输入结构没有 3D 坐标（conformer）")

        hop2_env_to_atom_param = load_hop_param_db(args.hop2_env_csv)
        hop1_env_to_atom_param = load_hop_param_db(args.hop1_env_csv)
        hop0_env_to_atom_param = load_hop_param_db(args.hop0_env_csv)
        idx_rev, idx_imp = load_multiatom_db(args.multiatom_csv)

        fallback_hops = _parse_fallback_hops(args.fallback_hops)
        missing_log_path = build_log_path
        atom_records, atom_type_rows, fallback_hit_counter = build_atom_types(
            mol,
            hop2_env_to_atom_param,
            hop1_env_to_atom_param,
            hop0_env_to_atom_param,
            args.hop_depth,
            fallback_hops=fallback_hops,
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
            strict_missing=False,
            build_log_path=build_log_path,
        )
        angle_records, angle_missing, angle_amb = assign_multiatom_params(
            "angle",
            angles,
            atom_key_sets,
            idx_rev,
            idx_imp,
            strict_missing=False,
            build_log_path=build_log_path,
        )
        dihedral_records, dihedral_missing, dihedral_amb = assign_multiatom_params(
            "dihedral",
            dihedrals,
            atom_key_sets,
            idx_rev,
            idx_imp,
            strict_missing=False,
            build_log_path=build_log_path,
        )
        improper_records, improper_missing, improper_amb = assign_multiatom_params(
            "improper",
            impropers,
            atom_key_sets,
            idx_rev,
            idx_imp,
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

    print("完成：")
    print(f"- 输出: {args.out}")
    print(f"- atoms: {len(atom_records)}")
    if fallback_hit_counter:
        detail = ", ".join(
            f"hop{hop}:{count}" for hop, count in fallback_hit_counter.items()
        )
        print(f"- atoms env_key 回退命中: {detail}")
    else:
        print("- atoms env_key 回退命中: 0")
    print(f"- bonds: {len(bond_records)}")
    print(f"- angles: {len(angle_records)}")
    print(f"- dihedrals: {len(dihedral_records)}")
    print(f"- impropers (匹配后): {len(improper_records)}")
    print(f"- impropers (候选未命中已删除): {len(improper_missing)}")


if __name__ == "__main__":
    main()
