#!/usr/bin/env python3
import csv
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem

from macromapff.pipeline.env_build import make_env_key, precompute_atom_context
from macromapff.pipeline.core.env import json_cell


ENV_INDEX_COLUMNS = [
    "z",
    "formal_charge",
    "aromatic",
    "hybridization",
    "degree",
    "total_hs",
    "in_ring",
    "ring_count",
    "hop1_shell",
    "hop2_shell",
    "neighbor_sig",
    "bond_kinds",
]


def load_input_structure(structure_path: Path) -> Chem.Mol:
    suffix = structure_path.suffix.lower()
    if suffix != ".mol":
        raise ValueError(f"Only .mol input is supported; got: {structure_path}")

    mol = Chem.MolFromMolFile(
        str(structure_path), removeHs=False, sanitize=True, strictParsing=False
    )
    if mol is None:
        mol = Chem.MolFromMolFile(
            str(structure_path), removeHs=False, sanitize=False, strictParsing=False
        )
        if mol is None:
            raise ValueError(f"RDKit failed to read .mol file: {structure_path}")
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


def _env_index_key_from_obj(obj: dict):
    key = []
    for col in ENV_INDEX_COLUMNS:
        key.append(json_cell(obj.get(col, "")))
    return tuple(key)


def _env_key_from_obj(obj: dict):
    norm = {}
    for col in ENV_INDEX_COLUMNS:
        val = obj.get(col, "")
        if val == "" or val is None:
            continue
        norm[col] = val
    return json.dumps(norm, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _env_features_at_hop(features: dict, hop: int):
    reduced = dict(features)
    for key in list(reduced.keys()):
        if key.startswith("hop") and key.endswith("_shell"):
            try:
                num = int(key[3:-6])
            except Exception:
                continue
            if num > hop:
                reduced[key] = ""
    return reduced


def load_hop_param_db(hop_csv: Path):
    if not hop_csv.exists():
        raise FileNotFoundError(f"Hop parameter database not found: {hop_csv}")

    structured_index = {}
    env_index = {}
    with hop_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            source_ids = [x for x in str(row.get("source_key_ids", "")).split(";") if x]
            key_ids = sorted({int(x) for x in source_ids}) if source_ids else [-1]
            payload = {
                "key_ids": key_ids,
                "charge": float(row["charge_mean"]),
                "sigma": float(row["sigma_mean"]),
                "epsilon": float(row["epsilon_mean"]),
                "mass": float(row["mass_mean"]),
            }
            idx_key = tuple(str(row.get(c, "") or "") for c in ENV_INDEX_COLUMNS)
            structured_index[idx_key] = payload

            env_key = str(row.get("env_key", "") or "").strip()
            if env_key:
                env_index[env_key] = payload

    if not structured_index:
        raise ValueError(f"Hop parameter database is empty: {hop_csv}")
    return {
        "structured": structured_index,
        "env": env_index,
    }


def parse_fallback_hops(raw: str):
    if raw is None:
        return tuple()
    items = []
    for token in str(raw).split(","):
        token = token.strip()
        if not token:
            continue
        hop = int(token)
        if hop < 0:
            raise ValueError(f"Fallback hop must be >= 0; got: {hop}")
        if hop not in items:
            items.append(hop)
    return tuple(items)


def append_build_log(log_path: Path, lines):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


def init_build_log(log_path: Path, lines):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


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
    struct_key_cache = {}

    def _lookup_params(db, features: dict, hop: int):
        cache_key = (id(features), hop)
        reduced = struct_key_cache.get(cache_key)
        if reduced is None:
            reduced = _env_features_at_hop(features, hop)
            struct_key_cache[cache_key] = reduced

        env_key = _env_key_from_obj(reduced)
        env_hit = db["env"].get(env_key)
        if env_hit is not None:
            return env_hit

        idx_key = _env_index_key_from_obj(reduced)
        return db["structured"].get(idx_key)

    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        env_key_raw, env_features = make_env_key(
            mol,
            atom,
            hop_depth=hop_depth,
            atom_ctx=atom_context.get(atom_idx),
        )
        params = _lookup_params(hop2_env_to_atom_param, env_features, 2)
        fallback_candidates = []

        matched_hop = 2
        if params is None:
            for fb_hop in fallback_hops:
                if fb_hop == 1:
                    params = _lookup_params(hop1_env_to_atom_param, env_features, 1)
                    fallback_candidates.append((1, "[structured-index]"))
                elif fb_hop == 0:
                    params = _lookup_params(hop0_env_to_atom_param, env_features, 0)
                    fallback_candidates.append((0, "[structured-index]"))
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
                    "env_key_exact": env_key_raw,
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
            "[ATOM] build_atom_types finished",
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
