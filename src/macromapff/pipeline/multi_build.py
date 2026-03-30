#!/usr/bin/env python3
"""Build a merged multi-atom master keymap from observed module data.

This script maps per-module ``lmp_type_tuple`` values into global ``key_type_tuple``
indices and merges rows by interaction kind and normalized key tuple.
"""

import argparse
import csv
import json
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
            f"Invalid --multiatom-spec format: {spec}\nExpected: module_name::multiatom_csv"
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
        raise FileNotFoundError(f"Type-stats file not found: {type_stats_csv}")

    mapping = {}
    with type_stats_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key_id = int(row["key_id"])
            module_name = row["module_name"]
            opls_type_id = int(row["opls_type_id"])
            mapping[(module_name, opls_type_id)] = key_id

    if not mapping:
        raise ValueError(f"Type-stats file is empty: {type_stats_csv}")
    return mapping


def load_hop0_key_classes(hop0_csv: Path):
    if not hop0_csv.exists():
        raise FileNotFoundError(f"hop0 CSV not found: {hop0_csv}")

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
            raise FileNotFoundError(f"multiatom CSV not found: {multiatom_csv}")

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
        "key_type_tuple",
        "coeff_param_sets",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Merge observed multiatom data into a global key_type-indexed master table."
    )
    parser.add_argument(
        "--type-stats-csv",
        required=True,
        help="Path to final_env_keymap_type_stats.csv",
    )
    parser.add_argument(
        "--multiatom-spec",
        action="append",
        required=True,
        help="Module input in the form module_name::multiatom_csv (repeatable).",
    )
    parser.add_argument(
        "--hop0-env-csv",
        required=True,
        help="Path to hop0_env_keymap.csv (used for key_id equivalence classes).",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix (without extension); writes .csv.",
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

    write_csv(out_csv, rows)

    print("Done:")
    print(f"- CSV: {out_csv}")
    print(f"- merged rows: {len(rows)}")
    print(f"- missing type refs: {len(missing_type_refs)}")


if __name__ == "__main__":
    main()
