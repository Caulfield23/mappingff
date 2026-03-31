#!/usr/bin/env python3
import math
from collections import Counter, defaultdict

from macromapff.domain.env import canonicalize_env_key, split_env_key_columns


def _new_stats():
    """Create an empty accumulator for running mean/std statistics."""
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
    """Accumulate one row's scalar values into running statistics."""
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
    """Compute mean and standard deviation from running sums."""
    if n <= 0:
        return None, None
    mean = sum_v / n
    var = max(0.0, sum2_v / n - mean * mean)
    return mean, math.sqrt(var)


def _round6(v):
    """Round numeric output to 6 decimals while preserving None."""
    if v is None:
        return None
    return round(float(v), 6)


def build_final_map(module_specs):
    """Merge all module atom-env rows into canonical env-key groups."""
    merged = {}

    for module_name, atom_env_csv in module_specs:
        if not atom_env_csv.exists():
            raise FileNotFoundError(f"atom_env.csv not found: {atom_env_csv}")

        with atom_env_csv.open("r", encoding="utf-8") as f:
            import csv

            reader = csv.DictReader(f)
            for row in reader:
                env_key = canonicalize_env_key(row["env_key"], deep_normalize=False)
                opls_type_id = int(row["opls_type_id"])
                opls_type_name = row["opls_type_name"]
                charge = float(row["charge"])
                sigma = float(row["sigma"])
                epsilon = float(row["epsilon"])
                mass_raw = str(row.get("mass", "") or "").strip()
                if not mass_raw:
                    raise ValueError(
                        f"Missing mass column value in atom_env row for module {module_name}: {atom_env_csv}"
                    )
                mass = float(mass_raw)

                if env_key not in merged:
                    merged[env_key] = {
                        "env_key": env_key,
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
                node["modules"].add(module_name)
                _add_stats(node["stats"], charge, sigma, epsilon, mass)

                type_key = (module_name, opls_type_id, opls_type_name)
                tnode = node["type_stats"][type_key]
                tnode["type_name"] = opls_type_name
                tnode["modules"].add(module_name)
                _add_stats(tnode["stats"], charge, sigma, epsilon, mass)

    return merged


def finalize_records(merged):
    """Finalize merged groups into export-ready final and type-level rows."""
    sorted_items = sorted(merged.items(), key=lambda x: x[0])
    final_rows = []
    type_rows = []

    for idx, (canonical_env_key, node) in enumerate(sorted_items, start=1):
        env_key = node["env_key"]
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
        global_type_ids = [f"{m}_{t}" for m, t, _ in lmp_type_keys]

        final_rows.append(
            {
                "key_id": idx,
                "global_type_ids": ";".join(global_type_ids),
                "env_key": env_key,
                **split_env_key_columns(env_key),
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
