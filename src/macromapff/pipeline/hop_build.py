#!/usr/bin/env python3
import csv
import json
from collections import defaultdict
from pathlib import Path

from macromapff.pipeline.core.env import ENV_SPLIT_COLUMNS, split_env_key_columns


def _row_index_key_at_hop(row: dict, hop: int):
    key_vals = []
    for col in ENV_SPLIT_COLUMNS:
        val = str(row.get(col, "") or "")
        if col.startswith("hop") and col.endswith("_shell"):
            try:
                num = int(col[3:-6])
            except Exception:
                num = -1
            if num > hop:
                val = ""
        key_vals.append(val)
    return tuple(key_vals)


def _env_key_from_split_cols(split_cols: dict):
    obj = {}
    for col in ENV_SPLIT_COLUMNS:
        raw = str(split_cols.get(col, "") or "")
        if raw == "":
            continue
        try:
            obj[col] = json.loads(raw)
        except Exception:
            obj[col] = raw
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _new_stats():
    return {
        "n": 0,
        "charge_sum": 0.0,
        "sigma_sum": 0.0,
        "epsilon_sum": 0.0,
        "mass_sum": 0.0,
        "source_key_ids": set(),
    }

def _add_stats(node, key_id, charge, sigma, epsilon, mass):
    node["n"] += 1
    node["charge_sum"] += charge
    node["sigma_sum"] += sigma
    node["epsilon_sum"] += epsilon
    node["mass_sum"] += mass
    node["source_key_ids"].add(str(key_id))


def build_hop_map(final_env_csv: Path, hop: int):
    out = defaultdict(_new_stats)

    with final_env_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx_key = _row_index_key_at_hop(row, hop)
            node = out[idx_key]
            _add_stats(
                node,
                key_id=row["key_id"],
                charge=float(row["charge_mean"]),
                sigma=float(row["sigma_mean"]),
                epsilon=float(row["epsilon_mean"]),
                mass=float(row["mass_mean"]),
            )

    rows = []
    for idx_key, node in sorted(out.items(), key=lambda x: x[0]):
        n = node["n"]
        split_cols = {ENV_SPLIT_COLUMNS[i]: idx_key[i] for i in range(len(ENV_SPLIT_COLUMNS))}
        rows.append(
            {
                "source_key_ids": ";".join(sorted(node["source_key_ids"], key=int)),
                "charge_mean": node["charge_sum"] / n,
                "sigma_mean": node["sigma_sum"] / n,
                "epsilon_mean": node["epsilon_sum"] / n,
                "mass_mean": node["mass_sum"] / n,
                "env_key": _env_key_from_split_cols(split_cols),
                **split_cols,
            }
        )
    return rows


def write_hop_csv(rows, out_csv: Path):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "source_key_ids",
                "charge_mean",
                "sigma_mean",
                "epsilon_mean",
                "mass_mean",
                "env_key",
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
                "hop1_shell",
                "hop2_shell",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "source_key_ids": row["source_key_ids"],
                    "charge_mean": f"{row['charge_mean']:.8f}",
                    "sigma_mean": f"{row['sigma_mean']:.8f}",
                    "epsilon_mean": f"{row['epsilon_mean']:.8f}",
                    "mass_mean": f"{row['mass_mean']:.8f}",
                    "env_key": row["env_key"],
                    "z": row["z"],
                    "formal_charge": row["formal_charge"],
                    "aromatic": row["aromatic"],
                    "hybridization": row["hybridization"],
                    "degree": row["degree"],
                    "total_hs": row["total_hs"],
                    "in_ring": row["in_ring"],
                    "ring_count": row["ring_count"],
                    "neighbor_sig": row["neighbor_sig"],
                    "bond_kinds": row["bond_kinds"],
                    "hop1_shell": row["hop1_shell"],
                    "hop2_shell": row["hop2_shell"],
                }
            )


def build_hop_databases(
    final_env_csv: Path,
    hop2_out: Path,
    hop1_out: Path,
    hop0_out: Path,
):
    hop2_rows = build_hop_map(final_env_csv, hop=2)
    hop1_rows = build_hop_map(final_env_csv, hop=1)
    hop0_rows = build_hop_map(final_env_csv, hop=0)

    write_hop_csv(hop2_rows, hop2_out)
    write_hop_csv(hop1_rows, hop1_out)
    write_hop_csv(hop0_rows, hop0_out)
    return hop2_rows, hop1_rows, hop0_rows


