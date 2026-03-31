#!/usr/bin/env python3
"""Build a merged multi-atom master keymap from observed module data.

This script maps per-module ``lmp_type_tuple`` values into global ``key_type_tuple``
indices and merges rows by interaction kind and normalized key tuple.
"""

import csv
import json
from collections import defaultdict
from pathlib import Path


INTERACTION_ORDER = {"bond": 0, "angle": 1, "dihedral": 2, "improper": 3}


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


def load_type_to_keyid(final_env_csv: Path):
    if not final_env_csv.exists():
        raise FileNotFoundError(f"Final env CSV not found: {final_env_csv}")

    mapping = {}
    with final_env_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key_id = int(row["key_id"])
            raw = str(row.get("global_type_ids", "") or "")
            for global_type in [x.strip() for x in raw.split(";") if x.strip()]:
                if global_type not in mapping or key_id < mapping[global_type]:
                    mapping[global_type] = key_id

    if not mapping:
        raise ValueError(f"Final env CSV has no global_type_ids: {final_env_csv}")
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
                    lookup_key = f"{module_name}_{int(type_id)}"
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

                # Merge only when coefficient parameter sets are exactly the same.
                # This prevents collapsing records that share hop0-equivalent key classes
                # but correspond to different force-field parameters.
                normalized_coeff_sets = []
                if isinstance(coeff_param_sets, list):
                    for coeff in coeff_param_sets:
                        if isinstance(coeff, list):
                            normalized_coeff_sets.append(tuple(str(x) for x in coeff))
                coeff_signature = json.dumps(
                    [list(x) for x in sorted(set(normalized_coeff_sets))],
                    ensure_ascii=False,
                    separators=(",", ":"),
                )

                key = (
                    interaction_kind,
                    n_atoms,
                    json.dumps(key_type_tuple, separators=(",", ":")),
                    coeff_signature,
                )
                node = merged[key]
                node["term_count"] += term_count
                node["source_modules"].add(module_name)
                for term_type in source_term_types:
                    node["source_term_types"].add(int(term_type))

    rows = []
    for (interaction_kind, n_atoms, key_type_tuple_json, coeff_signature), payload in merged.items():
        coeff_param_sets = coeff_signature

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
            INTERACTION_ORDER.get(r["interaction_kind"], 99),
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


def append_merge_conflict_log(log_path: Path, rows):
    # Same key_type_tuple with multiple coeff sets indicates parameter split kept intentionally.
    grouped = defaultdict(set)
    for row in rows:
        base_key = (row["interaction_kind"], row["key_type_tuple"])
        grouped[base_key].add(row["coeff_param_sets"])

    conflicts = [
        (kind, key_type_tuple, sorted(coeffs))
        for (kind, key_type_tuple), coeffs in grouped.items()
        if len(coeffs) > 1
    ]
    conflicts.sort(key=lambda x: (INTERACTION_ORDER.get(x[0], 99), x[1]))

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write("\n")
        f.write("=== Multiatom merge check ===\n")
        f.write(
            "rule: merge only when interaction_kind + key_type_tuple + coeff_param_sets are all identical\n"
        )
        f.write(f"rows in multiatom master: {len(rows)}\n")
        f.write(f"key_type groups with multiple coeff sets: {len(conflicts)}\n")

        if conflicts:
            f.write("--- examples (up to 20) ---\n")
            for idx, (kind, key_type_tuple, coeffs) in enumerate(conflicts[:20], start=1):
                f.write(
                    f"{idx}. kind={kind} key_type_tuple={key_type_tuple} coeff_set_count={len(coeffs)}\n"
                )
                for coeff in coeffs[:3]:
                    f.write(f"    coeff={coeff}\n")


def build_multiatom_master(
    final_env_csv: Path,
    hop0_env_csv: Path,
    multiatom_specs,
    out_prefix: Path,
    log_file: Path | None = None,
):
    type_to_keyid = load_type_to_keyid(final_env_csv)
    key_to_class = load_hop0_key_classes(hop0_env_csv)
    rows, missing_type_refs = build_master(type_to_keyid, key_to_class, multiatom_specs)

    out_prefix = Path(out_prefix).expanduser()
    out_csv = out_prefix.with_suffix(".csv")
    write_csv(out_csv, rows)

    if log_file is not None:
        append_merge_conflict_log(Path(log_file).expanduser(), rows)

    return out_csv, rows, missing_type_refs


