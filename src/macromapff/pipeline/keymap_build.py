#!/usr/bin/env python3
"""Build a final merged atom-environment keymap across modules.

Inputs:
- Repeated ``--module-spec`` values in the form
    ``module_name::atom_env_csv::lammps_data``

Outputs:
- ``{out_prefix}.csv``: merged final keymap
- ``{out_prefix}_type_stats.csv``: per-type statistics per key
- ``{out_prefix}.log``: merge diagnostics
"""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

try:
    from .core.env import (
        ENV_SPLIT_COLUMNS,
        canonicalize_env_key,
        split_env_key_columns,
    )
    from .core.lammps_parse import parse_lammps_masses
except ImportError:
    from core.env import ENV_SPLIT_COLUMNS, canonicalize_env_key, split_env_key_columns
    from core.lammps_parse import parse_lammps_masses


def parse_module_spec(spec: str):
    parts = spec.split("::")
    if len(parts) != 3:
        raise ValueError(
            f"Invalid --module-spec format: {spec}\n"
            f"Expected: module_name::atom_env_csv::lammps_data"
        )
    module, atom_env_csv, lmp_data = parts
    return module.strip(), Path(atom_env_csv).expanduser(), Path(lmp_data).expanduser()


def _new_stats():
    return {
        "n": 0,
        "charge_sum": 0.0,
        "charge_sum2": 0.0,
        "sigma_sum": 0.0,
        "sigma_sum2": 0.0,
        "epsilon_sum": 0.0,
        "epsilon_sum2": 0.0,
        "mass_sum": 0.0,
        "mass_sum2": 0.0,
    }


def _add_stats(stats, charge, sigma, epsilon, mass):
    stats["n"] += 1
    stats["charge_sum"] += charge
    stats["charge_sum2"] += charge * charge
    stats["sigma_sum"] += sigma
    stats["sigma_sum2"] += sigma * sigma
    stats["epsilon_sum"] += epsilon
    stats["epsilon_sum2"] += epsilon * epsilon
    stats["mass_sum"] += mass
    stats["mass_sum2"] += mass * mass


def _mean_std(sum_v, sum2_v, n):
    if n <= 0:
        return None, None
    mean = sum_v / n
    var = max(0.0, sum2_v / n - mean * mean)
    return mean, math.sqrt(var)


def _fmt_float(v):
    if v is None:
        return ""
    return f"{v:.8f}"


def _round6(v):
    if v is None:
        return None
    return round(float(v), 6)


def stable_env_hash(canonical_env_key: str):
    return hashlib.sha256(canonical_env_key.encode("utf-8")).hexdigest()[:16]


def build_final_map(module_specs):
    merged = {}

    for module_name, atom_env_csv, lmp_data in module_specs:
        if not atom_env_csv.exists():
            raise FileNotFoundError(f"atom_env.csv not found: {atom_env_csv}")
        if not lmp_data.exists():
            raise FileNotFoundError(f"LAMMPS data file not found: {lmp_data}")

        mass_map = parse_lammps_masses(lmp_data)

        with atom_env_csv.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                env_key = canonicalize_env_key(row["env_key"], deep_normalize=False)
                source_env_key_hash = row.get("env_key_hash") or stable_env_hash(env_key)
                opls_type_id = int(row["opls_type_id"])
                opls_type_name = row["opls_type_name"]
                charge = float(row["charge"])
                sigma = float(row["sigma"])
                epsilon = float(row["epsilon"])

                if opls_type_id not in mass_map:
                    raise ValueError(
                        f"Missing mass for type={opls_type_id} in {lmp_data}; cannot merge module {module_name}"
                    )
                mass = float(mass_map[opls_type_id])

                if env_key not in merged:
                    merged[env_key] = {
                        "env_key": env_key,
                        "source_hashes": set(),
                        "modules": set(),
                        "stats": _new_stats(),
                        "type_stats": defaultdict(
                            lambda: {
                                "type_name": "",
                                "modules": set(),
                                "stats": _new_stats(),
                            }
                        ),
                    }

                node = merged[env_key]
                node["source_hashes"].add(source_env_key_hash)
                node["modules"].add(module_name)
                _add_stats(node["stats"], charge, sigma, epsilon, mass)

                type_key = (module_name, opls_type_id, opls_type_name)
                tnode = node["type_stats"][type_key]
                tnode["type_name"] = opls_type_name
                tnode["modules"].add(module_name)
                _add_stats(tnode["stats"], charge, sigma, epsilon, mass)

    return merged


def finalize_records(merged):
    sorted_items = sorted(merged.items(), key=lambda x: stable_env_hash(x[0]))
    final_rows = []
    type_rows = []

    for idx, (canonical_env_key, node) in enumerate(sorted_items, start=1):
        env_key = node["env_key"]
        env_key_hash = stable_env_hash(canonical_env_key)
        modules = sorted(node["modules"])

        type_counts = Counter()
        for type_key, payload in node["type_stats"].items():
            n = payload["stats"]["n"]
            type_counts[type_key] = n

        canonical_type_key, canonical_count = sorted(
            type_counts.items(), key=lambda x: (-x[1], x[0])
        )[0]
        canonical_module, canonical_type_id_raw, canonical_type_name_raw = (
            canonical_type_key
        )
        canonical_type_id = f"{canonical_module}:{canonical_type_id_raw}"
        canonical_type_name = f"{canonical_module}:{canonical_type_name_raw}"

        total_n = node["stats"]["n"]
        charge_mean, charge_std = _mean_std(
            node["stats"]["charge_sum"], node["stats"]["charge_sum2"], total_n
        )
        sigma_mean, sigma_std = _mean_std(
            node["stats"]["sigma_sum"], node["stats"]["sigma_sum2"], total_n
        )
        epsilon_mean, epsilon_std = _mean_std(
            node["stats"]["epsilon_sum"], node["stats"]["epsilon_sum2"], total_n
        )
        mass_mean, mass_std = _mean_std(
            node["stats"]["mass_sum"], node["stats"]["mass_sum2"], total_n
        )

        lmp_type_keys = sorted(
            node["type_stats"].keys(), key=lambda x: (x[0], x[1], x[2])
        )
        lmp_type_ids = [f"{m}:{t}" for m, t, _ in lmp_type_keys]
        lmp_type_names = [f"{m}:{n}" for m, _, n in lmp_type_keys]

        final_rows.append(
            {
                "key_id": idx,
                "env_key_hash": env_key_hash,
                "env_key": env_key,
                **split_env_key_columns(env_key),
                "n_source_hashes": len(node["source_hashes"]),
                "source_hashes": ";".join(sorted(node["source_hashes"])),
                "n_modules": len(modules),
                "modules": ";".join(modules),
                "n_rows": total_n,
                "n_lmp_types": len(lmp_type_ids),
                "lmp_type_ids": ";".join(lmp_type_ids),
                "lmp_type_names": ";".join(lmp_type_names),
                "canonical_lmp_type_id": canonical_type_id,
                "canonical_lmp_type_name": canonical_type_name,
                "canonical_count": canonical_count,
                "charge_mean": _round6(charge_mean),
                "charge_std": charge_std,
                "sigma_mean": _round6(sigma_mean),
                "sigma_std": sigma_std,
                "epsilon_mean": _round6(epsilon_mean),
                "epsilon_std": epsilon_std,
                "mass_mean": _round6(mass_mean),
                "mass_std": mass_std,
            }
        )

        for module_name, type_id, type_name in lmp_type_keys:
            payload = node["type_stats"][(module_name, type_id, type_name)]
            n = payload["stats"]["n"]
            c_mean, c_std = _mean_std(
                payload["stats"]["charge_sum"], payload["stats"]["charge_sum2"], n
            )
            s_mean, s_std = _mean_std(
                payload["stats"]["sigma_sum"], payload["stats"]["sigma_sum2"], n
            )
            e_mean, e_std = _mean_std(
                payload["stats"]["epsilon_sum"], payload["stats"]["epsilon_sum2"], n
            )
            m_mean, m_std = _mean_std(
                payload["stats"]["mass_sum"], payload["stats"]["mass_sum2"], n
            )
            type_rows.append(
                {
                    "key_id": idx,
                    "env_key_hash": env_key_hash,
                    "module_name": module_name,
                    "opls_type_id": type_id,
                    "opls_type_name": type_name,
                    "modules": ";".join(sorted(payload["modules"])),
                    "n_rows": n,
                    "charge_mean": _round6(c_mean),
                    "charge_std": c_std,
                    "sigma_mean": _round6(s_mean),
                    "sigma_std": s_std,
                    "epsilon_mean": _round6(e_mean),
                    "epsilon_std": e_std,
                    "mass_mean": _round6(m_mean),
                    "mass_std": m_std,
                }
            )

    return final_rows, type_rows


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_log(path: Path, module_specs, final_rows, type_rows):
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append("=== Final keymap merge log ===")
    lines.append(f"modules: {len(module_specs)}")
    for module_name, atom_env_csv, lmp_data in module_specs:
        lines.append(f"  - {module_name}: atom_env={atom_env_csv}, lmp={lmp_data}")

    lines.append("")
    lines.append(f"total keys: {len(final_rows)}")

    multi_type_keys = [r for r in final_rows if r["n_lmp_types"] > 1]
    lines.append(f"keys with multiple LMP types: {len(multi_type_keys)}")
    hash_collision_keys = [r for r in final_rows if r.get("n_source_hashes", 1) > 1]
    lines.append(
        f"merged-from-multiple-source-hashes (same env_key normalized): {len(hash_collision_keys)}"
    )

    lines.append("")
    lines.append("--- key-level stats (mean ± std) ---")
    for row in final_rows:
        lines.append(
            f"key_id={row['key_id']} hash={row['env_key_hash']} "
            f"types={row['n_lmp_types']} rows={row['n_rows']} "
            f"charge={_fmt_float(row['charge_mean'])}±{_fmt_float(row['charge_std'])} "
            f"sigma={_fmt_float(row['sigma_mean'])}±{_fmt_float(row['sigma_std'])} "
            f"epsilon={_fmt_float(row['epsilon_mean'])}±{_fmt_float(row['epsilon_std'])} "
            f"mass={_fmt_float(row['mass_mean'])}±{_fmt_float(row['mass_std'])}"
        )

    lines.append("")
    lines.append("--- multi-type key details ---")
    type_rows_by_key = defaultdict(list)
    for tr in type_rows:
        type_rows_by_key[tr["key_id"]].append(tr)

    for row in multi_type_keys:
        lines.append(
            f"key_id={row['key_id']} hash={row['env_key_hash']} "
            f"lmp_type_ids=[{row['lmp_type_ids']}]"
        )
        for tr in sorted(
            type_rows_by_key[row["key_id"]],
            key=lambda x: (x["module_name"], x["opls_type_id"]),
        ):
            lines.append(
                f"    type={tr['module_name']}:{tr['opls_type_id']}({tr['opls_type_name']}), n={tr['n_rows']}, "
                f"charge={_fmt_float(tr['charge_mean'])}±{_fmt_float(tr['charge_std'])}, "
                f"sigma={_fmt_float(tr['sigma_mean'])}±{_fmt_float(tr['sigma_std'])}, "
                f"epsilon={_fmt_float(tr['epsilon_mean'])}±{_fmt_float(tr['epsilon_std'])}, "
                f"mass={_fmt_float(tr['mass_mean'])}±{_fmt_float(tr['mass_std'])}, "
                f"modules={tr['modules']}"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Merge atom_env outputs from multiple modules into a final keymap with statistics."
    )
    parser.add_argument(
        "--module-spec",
        action="append",
        required=True,
        help="Module input in the form module_name::atom_env_csv::lammps_data (repeatable).",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix (without extension); writes .csv and .log files.",
    )
    args = parser.parse_args()

    module_specs = [parse_module_spec(spec) for spec in args.module_spec]
    out_prefix = Path(args.out_prefix).expanduser()

    merged = build_final_map(module_specs)
    final_rows, type_rows = finalize_records(merged)

    final_csv = out_prefix.with_suffix(".csv")
    final_log = out_prefix.with_suffix(".log")
    type_csv = out_prefix.parent / f"{out_prefix.name}_type_stats.csv"

    write_csv(
        final_csv,
        final_rows,
        [
            "key_id",
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
            "charge_mean",
            "sigma_mean",
            "epsilon_mean",
            "mass_mean",
            "hop1_shell",
            "hop2_shell",
        ],
    )

    write_csv(
        type_csv,
        type_rows,
        [
            "key_id",
            "module_name",
            "opls_type_id",
        ],
    )

    write_log(final_log, module_specs, final_rows, type_rows)

    print("Done:")
    print(f"- final keymap CSV: {final_csv}")
    print(f"- type stats CSV: {type_csv}")
    print(f"- merge log: {final_log}")
    print(f"- total keys: {len(final_rows)}")


if __name__ == "__main__":
    main()
