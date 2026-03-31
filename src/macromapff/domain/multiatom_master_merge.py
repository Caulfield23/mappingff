#!/usr/bin/env python3
import json
from collections import defaultdict

INTERACTION_ORDER = {"bond": 0, "angle": 1, "dihedral": 2, "improper": 3}


def _parse_json_field(raw: str):
    """Parse JSON field content and return empty list on failure."""
    try:
        return json.loads(raw)
    except Exception:
        return []


def _canonicalize_key_type_tuple(interaction_kind: str, key_type_tuple):
    """Canonicalize key-type tuple, normalizing improper permutations."""
    if interaction_kind != "improper" or len(key_type_tuple) != 4:
        return key_type_tuple

    center = key_type_tuple[0]
    others = sorted(
        key_type_tuple[1:],
        key=lambda x: json.dumps(x, ensure_ascii=False, separators=(",", ":")),
    )
    return [center] + others


def build_master(type_to_keyid, key_to_class, multiatom_specs):
    """Merge module-level observed terms into master key-type mappings."""
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
            import csv

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
