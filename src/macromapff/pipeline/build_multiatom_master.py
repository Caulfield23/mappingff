#!/usr/bin/env python3
"""
脚本说明
--------
把各 segment 的 multiatom observed 表合并成总表，并把 lmp_type_tuple 转换为 key_type 索引（key_id）。

合并规则
--------
- 先通过 final_env_keymap_type_stats.csv 构建映射：(module_name, opls_type_id) -> key_id
- 将每条 multiatom 的 lmp_type_tuple（如 [88,87,110]）映射为 key_type_tuple（如 [12,34,56]）
- 按以下键合并：
    interaction_kind + n_atoms + key_type_tuple
- term_count 求和；source_term_types 去重汇总
- coeff 采用众数（按 term_count 加权计数），并写成单一 coeff_param_sets

注意
----
- 当前不做取反/对称归一化（例如 bond 的 [1,2] 与 [2,1]），与你的要求一致。

输出
----
- {out_prefix}.csv
- {out_prefix}.json
- {out_prefix}.sqlite

用法示例
--------
python scripts/build_multiatom_master.py \
  --type-stats-csv outputs/final_env_keymap_type_stats.csv \
  --multiatom-spec 'segment1::outputs/segment1_multiatom_observed.csv' \
  --multiatom-spec 'segment2::outputs/segment2_multiatom_observed.csv' \
  --multiatom-spec 'segment3::outputs/segment3_multiatom_observed.csv' \
  --out-prefix outputs/multiatom_master_keytype
"""

import argparse
import csv
import json
import sqlite3
from collections import defaultdict
from pathlib import Path


class DSU:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def parse_multiatom_spec(spec: str):
    parts = spec.split("::")
    if len(parts) != 2:
        raise ValueError(
            f"--multiatom-spec 格式错误: {spec}\n" f"应为: module_name::multiatom_csv"
        )
    module_name, csv_path = parts
    return module_name.strip(), Path(csv_path).expanduser()


def _parse_json_field(raw: str):
    try:
        return json.loads(raw)
    except Exception:
        return []


def _canonicalize_key_type_tuple(interaction_kind: str, key_type_tuple):
    if interaction_kind != "improper" or len(key_type_tuple) != 4:
        return key_type_tuple

    center = key_type_tuple[0]
    others = sorted(
        key_type_tuple[1:],
        key=lambda x: json.dumps(x, ensure_ascii=False, separators=(",", ":")),
    )
    return [center] + others


def load_type_to_keyid(type_stats_csv: Path):
    if not type_stats_csv.exists():
        raise FileNotFoundError(f"找不到 type stats 文件: {type_stats_csv}")

    mapping = {}
    with type_stats_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key_id = int(row["key_id"])
            module_name = row["module_name"]
            opls_type_id = int(row["opls_type_id"])
            mapping[(module_name, opls_type_id)] = key_id

    if not mapping:
        raise ValueError(f"type stats 为空: {type_stats_csv}")
    return mapping


def load_hop0_key_classes(hop0_csv: Path):
    if not hop0_csv.exists():
        raise FileNotFoundError(f"找不到 hop0 csv: {hop0_csv}")

    dsu = DSU()
    seen_ids = set()
    with hop0_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = str(row.get("source_key_ids", "")).strip()
            if not raw:
                continue
            ids = [int(x) for x in raw.split(";") if x]
            if not ids:
                continue
            for key_id in ids:
                seen_ids.add(key_id)
                dsu.find(key_id)
            head = ids[0]
            for key_id in ids[1:]:
                dsu.union(head, key_id)

    groups = defaultdict(set)
    for key_id in seen_ids:
        groups[dsu.find(key_id)].add(key_id)

    key_to_class = {}
    for _, members in groups.items():
        cls = sorted(members)
        for key_id in members:
            key_to_class[key_id] = cls
    return key_to_class


def build_master(type_to_keyid, key_to_class, multiatom_specs):
    merged = defaultdict(
        lambda: {
            "term_count": 0,
            "source_term_types": set(),
            "source_modules": set(),
            "coeff_counter": defaultdict(int),
        }
    )

    missing_type_refs = []

    for module_name, multiatom_csv in multiatom_specs:
        if not multiatom_csv.exists():
            raise FileNotFoundError(f"找不到 multiatom csv: {multiatom_csv}")

        with multiatom_csv.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                interaction_kind = row["interaction_kind"]
                n_atoms = int(row["n_atoms"])
                lmp_type_tuple = _parse_json_field(row["lmp_type_tuple"])
                coeff_param_sets = _parse_json_field(row.get("coeff_param_sets", "[]"))
                term_count = int(row.get("term_count", 0))
                source_term_types = _parse_json_field(
                    row.get("source_term_types", "[]")
                )

                key_type_tuple = []
                missing = False
                for type_id in lmp_type_tuple:
                    lookup_key = (module_name, int(type_id))
                    if lookup_key not in type_to_keyid:
                        missing = True
                        missing_type_refs.append(
                            {
                                "module_name": module_name,
                                "opls_type_id": int(type_id),
                                "interaction_kind": interaction_kind,
                                "n_atoms": n_atoms,
                                "lmp_type_tuple": lmp_type_tuple,
                            }
                        )
                        break
                    key_id = type_to_keyid[lookup_key]
                    key_type_tuple.append(key_to_class.get(key_id, [key_id]))

                if missing:
                    continue

                key_type_tuple = _canonicalize_key_type_tuple(
                    interaction_kind, key_type_tuple
                )

                key = (
                    interaction_kind,
                    n_atoms,
                    json.dumps(key_type_tuple, separators=(",", ":")),
                )
                node = merged[key]
                node["term_count"] += term_count
                node["source_modules"].add(module_name)
                for term_type in source_term_types:
                    node["source_term_types"].add(int(term_type))

                if isinstance(coeff_param_sets, list):
                    for coeff in coeff_param_sets:
                        if isinstance(coeff, list):
                            coeff_tuple = tuple(str(x) for x in coeff)
                            node["coeff_counter"][coeff_tuple] += max(term_count, 1)

    rows = []
    for (interaction_kind, n_atoms, key_type_tuple_json), payload in merged.items():
        coeff_counter = payload["coeff_counter"]
        if coeff_counter:
            mode_coeff = sorted(
                coeff_counter.items(),
                key=lambda kv: (-kv[1], kv[0]),
            )[
                0
            ][0]
            coeff_param_sets = json.dumps([list(mode_coeff)], ensure_ascii=False)
        else:
            coeff_param_sets = json.dumps([], ensure_ascii=False)

        rows.append(
            {
                "interaction_kind": interaction_kind,
                "n_atoms": n_atoms,
                "key_type_tuple": key_type_tuple_json,
                "term_count": payload["term_count"],
                "source_modules": ";".join(sorted(payload["source_modules"])),
                "source_term_types": json.dumps(
                    sorted(payload["source_term_types"]), ensure_ascii=False
                ),
                "coeff_param_sets": coeff_param_sets,
            }
        )

    rows.sort(
        key=lambda r: (
            r["interaction_kind"],
            r["key_type_tuple"],
            r["coeff_param_sets"],
        )
    )
    return rows, missing_type_refs


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "interaction_kind",
        "n_atoms",
        "key_type_tuple",
        "term_count",
        "source_modules",
        "source_term_types",
        "coeff_param_sets",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, rows, missing_type_refs):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "n_rows": len(rows),
        "rows": rows,
        "n_missing_type_refs": len(missing_type_refs),
        "missing_type_refs": missing_type_refs[:200],
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def write_sqlite(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    conn = sqlite3.connect(str(path))
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE multiatom_master (
            interaction_kind TEXT,
            n_atoms INTEGER,
            key_type_tuple TEXT,
            term_count INTEGER,
            source_modules TEXT,
            source_term_types TEXT,
            coeff_param_sets TEXT
        )
        """
    )

    for row in rows:
        cur.execute(
            """
            INSERT INTO multiatom_master VALUES (?,?,?,?,?,?,?)
            """,
            (
                row["interaction_kind"],
                row["n_atoms"],
                row["key_type_tuple"],
                row["term_count"],
                row["source_modules"],
                row["source_term_types"],
                row["coeff_param_sets"],
            ),
        )

    conn.commit()
    conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="合并 multiatom observed 为 key_type 索引总表"
    )
    parser.add_argument(
        "--type-stats-csv",
        required=True,
        help="final_env_keymap_type_stats.csv 路径",
    )
    parser.add_argument(
        "--multiatom-spec",
        action="append",
        required=True,
        help="模块输入，格式: module_name::multiatom_csv，可重复传入",
    )
    parser.add_argument(
        "--hop0-env-csv",
        required=True,
        help="hop0_env_keymap.csv 路径（用于 key_id 等价类）",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="输出前缀（不带后缀），将生成 .csv/.json/.sqlite",
    )
    args = parser.parse_args()

    type_stats_csv = Path(args.type_stats_csv).expanduser()
    hop0_env_csv = Path(args.hop0_env_csv).expanduser()
    multiatom_specs = [parse_multiatom_spec(s) for s in args.multiatom_spec]
    out_prefix = Path(args.out_prefix).expanduser()

    type_to_keyid = load_type_to_keyid(type_stats_csv)
    key_to_class = load_hop0_key_classes(hop0_env_csv)
    rows, missing_type_refs = build_master(type_to_keyid, key_to_class, multiatom_specs)

    out_csv = out_prefix.with_suffix(".csv")
    out_json = out_prefix.with_suffix(".json")
    out_sqlite = out_prefix.with_suffix(".sqlite")

    write_csv(out_csv, rows)
    write_json(out_json, rows, missing_type_refs)
    write_sqlite(out_sqlite, rows)

    print("完成：")
    print(f"- CSV: {out_csv}")
    print(f"- JSON: {out_json}")
    print(f"- SQLite: {out_sqlite}")
    print(f"- merged rows: {len(rows)}")
    print(f"- missing type refs: {len(missing_type_refs)}")


if __name__ == "__main__":
    main()
