#!/usr/bin/env python3
"""
脚本说明
--------
提取 LAMMPS data 文件中的多原子项（bond/angle/dihedral/improper），并结合
atom_env 构建“观测到的”env_key 组合 -> LMP 原子类型组合映射表。

核心特性
--------
- 支持四类多原子项：bond, angle, dihedral, improper
- 对 bond/angle/dihedral 做方向规范化（正反等价）
- 仅保留 LAMMPS 中真实出现的组合（无候选展开、无稀疏空值）

输出文件
--------
- {out_prefix}.csv：展开后的组合表
- {out_prefix}.json：同内容 JSON + env 组合摘要
- {out_prefix}.sqlite：`expanded_multiatom_map` 与 `env_tuple_summary` 两张表

用法示例
--------
python scripts/extract_multiatom_terms.py \
  --lmp /path/to/segment1.lammps.lmp \
  --atom-env-csv /path/to/segment1_atom_env.csv \
  --out-prefix /path/to/outputs/segment1_multiatom_map
"""

import argparse
import csv
import hashlib
import json
import sqlite3
from pathlib import Path


def _parse_int_tokens(line: str):
    body = line.split("#", 1)[0].strip()
    if not body:
        return []
    toks = body.split()
    if not toks:
        return []
    return toks


def parse_lammps_topology_and_coeffs(lmp_path: Path):
    lines = lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    coeff_sections = {
        "Bond Coeffs": "bond",
        "Angle Coeffs": "angle",
        "Dihedral Coeffs": "dihedral",
        "Improper Coeffs": "improper",
    }
    topo_sections = {
        "Bonds": ("bond", 2),
        "Angles": ("angle", 3),
        "Dihedrals": ("dihedral", 4),
        "Impropers": ("improper", 4),
    }
    known_sections = (
        set(coeff_sections)
        | set(topo_sections)
        | {
            "Masses",
            "Pair Coeffs",
            "Atoms",
            "Velocities",
        }
    )

    coeffs = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
    terms = {"bond": [], "angle": [], "dihedral": [], "improper": []}

    current = None
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped in known_sections:
            current = stripped
            continue
        if stripped.startswith("#"):
            continue

        if current in coeff_sections:
            toks = _parse_int_tokens(stripped)
            if len(toks) < 2:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            kind = coeff_sections[current]
            type_id = int(toks[0])
            params = toks[1:]
            coeffs[kind][type_id] = params
            continue

        if current in topo_sections:
            kind, n_atoms = topo_sections[current]
            toks = _parse_int_tokens(stripped)
            if len(toks) < 2 + n_atoms:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            term_id = int(toks[0])
            term_type = int(toks[1])
            atom_ids = [int(x) for x in toks[2 : 2 + n_atoms]]
            terms[kind].append(
                {
                    "term_id": term_id,
                    "term_type": term_type,
                    "atom_ids": atom_ids,
                }
            )

    return coeffs, terms


def load_atom_env(atom_env_csv: Path):
    atom_map = {}
    module_name = ""
    with atom_env_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            atom_index = int(row["atom_index"])
            atom_map[atom_index] = {
                "env_key": row["env_key"],
                "lmp_type": int(row["opls_type_id"]),
            }
            if not module_name:
                module_name = row.get("module", "")
    if not atom_map:
        raise ValueError(f"atom_env 为空或不可解析: {atom_env_csv}")
    return module_name, atom_map


def canonicalize_tuple(kind: str, env_tuple: tuple, type_tuple: tuple):
    if kind in {"bond", "angle", "dihedral"}:
        fwd = (env_tuple, type_tuple)
        rev = (tuple(reversed(env_tuple)), tuple(reversed(type_tuple)))
        return fwd if fwd <= rev else rev
    if kind == "improper" and len(env_tuple) == 4 and len(type_tuple) == 4:
        center_env = env_tuple[0]
        center_type = type_tuple[0]
        others = sorted(zip(env_tuple[1:], type_tuple[1:]))
        env_new = (center_env,) + tuple(x[0] for x in others)
        type_new = (center_type,) + tuple(x[1] for x in others)
        return env_new, type_new
    return env_tuple, type_tuple


def env_tuple_hash(env_tuple: tuple):
    payload = "|".join(env_tuple)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


def build_observed_mapping(coeffs, terms, atom_map):
    observed_records = {}
    env_summary = {}

    for kind in ["bond", "angle", "dihedral", "improper"]:
        for term in terms[kind]:
            atom_ids = term["atom_ids"]
            term_type = term["term_type"]
            coeff_params = coeffs[kind].get(term_type, [])

            try:
                env_tuple_raw = tuple(atom_map[a]["env_key"] for a in atom_ids)
                obs_types_raw = tuple(atom_map[a]["lmp_type"] for a in atom_ids)
            except KeyError as exc:
                raise ValueError(
                    f"atom_env 缺少原子索引 {exc.args[0]}，无法构建多原子映射。"
                )

            can_env_obs, can_types_obs = canonicalize_tuple(
                kind, env_tuple_raw, obs_types_raw
            )

            summary_key = (kind, can_env_obs)
            if summary_key not in env_summary:
                env_summary[summary_key] = {
                    "interaction_kind": kind,
                    "n_atoms": len(can_env_obs),
                    "env_key_hash": env_tuple_hash(can_env_obs),
                    "n_terms": 0,
                    "observed_type_tuples": set(),
                    "observed_term_types": set(),
                }
            env_summary[summary_key]["n_terms"] += 1
            env_summary[summary_key]["observed_type_tuples"].add(can_types_obs)
            env_summary[summary_key]["observed_term_types"].add(term_type)

            rec_key = (kind, can_env_obs, can_types_obs)
            if rec_key not in observed_records:
                observed_records[rec_key] = {
                    "interaction_kind": kind,
                    "n_atoms": len(can_env_obs),
                    "env_key_hash": env_tuple_hash(can_env_obs),
                    "lmp_type_tuple": list(can_types_obs),
                    "term_count": 0,
                    "source_term_types": set(),
                    "coeff_param_sets": set(),
                }

            rec = observed_records[rec_key]
            rec["term_count"] += 1
            rec["source_term_types"].add(term_type)
            rec["coeff_param_sets"].add(tuple(coeff_params))

    observed_list = []
    for _, rec in sorted(
        observed_records.items(),
        key=lambda x: (
            x[1]["interaction_kind"],
            x[1]["env_key_hash"],
            tuple(x[1]["lmp_type_tuple"]),
        ),
    ):
        observed_list.append(
            {
                "interaction_kind": rec["interaction_kind"],
                "n_atoms": rec["n_atoms"],
                "env_key_hash": rec["env_key_hash"],
                "lmp_type_tuple": rec["lmp_type_tuple"],
                "term_count": int(rec["term_count"]),
                "source_term_types": sorted(rec["source_term_types"]),
                "coeff_param_sets": [list(x) for x in sorted(rec["coeff_param_sets"])],
            }
        )

    summary_list = []
    for _, s in sorted(
        env_summary.items(),
        key=lambda x: (x[1]["interaction_kind"], x[1]["env_key_hash"]),
    ):
        summary_list.append(
            {
                "interaction_kind": s["interaction_kind"],
                "n_atoms": s["n_atoms"],
                "env_key_hash": s["env_key_hash"],
                "n_terms": int(s["n_terms"]),
                "n_observed_type_tuples": len(s["observed_type_tuples"]),
                "observed_type_tuples": [
                    list(x) for x in sorted(s["observed_type_tuples"])
                ],
                "observed_term_types": sorted(s["observed_term_types"]),
            }
        )

    return observed_list, summary_list


def write_csv(csv_path: Path, observed_list):
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "interaction_kind",
                "n_atoms",
                "env_key_hash",
                "lmp_type_tuple",
                "term_count",
                "source_term_types",
                "coeff_param_sets",
            ],
        )
        writer.writeheader()
        for row in observed_list:
            writer.writerow(
                {
                    "interaction_kind": row["interaction_kind"],
                    "n_atoms": row["n_atoms"],
                    "env_key_hash": row["env_key_hash"],
                    "lmp_type_tuple": json.dumps(
                        row["lmp_type_tuple"], ensure_ascii=False, separators=(",", ":")
                    ),
                    "term_count": row["term_count"],
                    "source_term_types": json.dumps(
                        row["source_term_types"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                    "coeff_param_sets": json.dumps(
                        row["coeff_param_sets"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                }
            )


def write_json(json_path: Path, module_name: str, observed_list, summary_list):
    json_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "module": module_name,
        "n_observed_rows": len(observed_list),
        "n_env_tuple_groups": len(summary_list),
        "observed_multiatom_map": observed_list,
        "env_tuple_summary": summary_list,
    }
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def write_sqlite(db_path: Path, observed_list, summary_list):
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS observed_multiatom_map")
    cur.execute("DROP TABLE IF EXISTS env_tuple_summary")

    cur.execute(
        """
        CREATE TABLE observed_multiatom_map (
            interaction_kind TEXT,
            n_atoms INTEGER,
            env_key_hash TEXT,
            lmp_type_tuple TEXT,
            term_count INTEGER,
            source_term_types TEXT,
            coeff_param_sets TEXT,
            PRIMARY KEY (interaction_kind, env_key_hash, lmp_type_tuple)
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE env_tuple_summary (
            interaction_kind TEXT,
            n_atoms INTEGER,
            env_key_hash TEXT,
            n_terms INTEGER,
            n_observed_type_tuples INTEGER,
            observed_type_tuples TEXT,
            observed_term_types TEXT,
            PRIMARY KEY (interaction_kind, env_key_hash)
        )
        """
    )

    cur.executemany(
        """
        INSERT INTO observed_multiatom_map (
            interaction_kind, n_atoms, env_key_hash, lmp_type_tuple,
            term_count, source_term_types, coeff_param_sets
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                row["interaction_kind"],
                int(row["n_atoms"]),
                row["env_key_hash"],
                json.dumps(
                    row["lmp_type_tuple"], ensure_ascii=False, separators=(",", ":")
                ),
                int(row["term_count"]),
                json.dumps(
                    row["source_term_types"], ensure_ascii=False, separators=(",", ":")
                ),
                json.dumps(
                    row["coeff_param_sets"], ensure_ascii=False, separators=(",", ":")
                ),
            )
            for row in observed_list
        ],
    )

    cur.executemany(
        """
        INSERT INTO env_tuple_summary (
            interaction_kind, n_atoms, env_key_hash, n_terms,
            n_observed_type_tuples,
            observed_type_tuples, observed_term_types
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                row["interaction_kind"],
                int(row["n_atoms"]),
                row["env_key_hash"],
                int(row["n_terms"]),
                int(row["n_observed_type_tuples"]),
                json.dumps(
                    row["observed_type_tuples"],
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                json.dumps(
                    row["observed_term_types"],
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
            )
            for row in summary_list
        ],
    )

    conn.commit()
    conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="提取 bond/angle/dihedral/improper，并建立仅观测组合的 env_key_hash -> LMP 类型组合映射表。",
        epilog=(
            "示例:\n"
            "  python scripts/extract_multiatom_terms.py "
            "--lmp /path/to/segment1.lammps.lmp "
            "--atom-env-csv /path/to/segment1_atom_env.csv "
            "--out-prefix /path/to/outputs/segment1_multiatom_map"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--lmp", required=True, type=Path, help="模块 LAMMPS data 文件")
    parser.add_argument(
        "--atom-env-csv",
        required=True,
        type=Path,
        help="build_envkey_mapping 生成的 atom_env.csv",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        type=Path,
        help="输出前缀（不带后缀），将生成 .csv/.json/.sqlite",
    )

    args = parser.parse_args()

    coeffs, terms = parse_lammps_topology_and_coeffs(args.lmp)
    module_name, atom_map = load_atom_env(args.atom_env_csv)

    observed_list, summary_list = build_observed_mapping(
        coeffs=coeffs,
        terms=terms,
        atom_map=atom_map,
    )

    out_csv = args.out_prefix.with_suffix(".csv")
    out_json = args.out_prefix.with_suffix(".json")
    out_db = args.out_prefix.with_suffix(".sqlite")

    write_csv(out_csv, observed_list)
    write_json(out_json, module_name, observed_list, summary_list)
    write_sqlite(out_db, observed_list, summary_list)

    print("完成：")
    print(f"- observed rows: {len(observed_list)}")
    print(f"- env tuple groups: {len(summary_list)}")
    print(f"- CSV: {out_csv}")
    print(f"- JSON: {out_json}")
    print(f"- SQLite: {out_db}")


if __name__ == "__main__":
    main()
