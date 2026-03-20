#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


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


def _env_key_at_hop_or_none(env_key_raw: str, hop: int):
    try:
        obj = json.loads(env_key_raw)
    except Exception:
        return None
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
            full_key = canonicalize_env_key(row["env_key"])
            reduced_key = _env_key_at_hop_or_none(full_key, hop)
            if not reduced_key:
                continue
            node = out[reduced_key]
            _add_stats(
                node,
                key_id=row["key_id"],
                charge=float(row["charge_mean"]),
                sigma=float(row["sigma_mean"]),
                epsilon=float(row["epsilon_mean"]),
                mass=float(row["mass_mean"]),
            )

    rows = []
    for env_key, node in sorted(out.items(), key=lambda x: x[0]):
        n = node["n"]
        rows.append(
            {
                "source_key_ids": ";".join(sorted(node["source_key_ids"], key=int)),
                "match_count": n,
                "charge_mean": node["charge_sum"] / n,
                "sigma_mean": node["sigma_sum"] / n,
                "epsilon_mean": node["epsilon_sum"] / n,
                "mass_mean": node["mass_sum"] / n,
                "env_key": env_key,
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
                "match_count",
                "charge_mean",
                "sigma_mean",
                "epsilon_mean",
                "mass_mean",
                "env_key",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "source_key_ids": row["source_key_ids"],
                    "match_count": row["match_count"],
                    "charge_mean": f"{row['charge_mean']:.8f}",
                    "sigma_mean": f"{row['sigma_mean']:.8f}",
                    "epsilon_mean": f"{row['epsilon_mean']:.8f}",
                    "mass_mean": f"{row['mass_mean']:.8f}",
                    "env_key": row["env_key"],
                }
            )


def main():
    parser = argparse.ArgumentParser(
        description="从 final_env_keymap.csv 聚合构建 hop2、hop1 与 hop0（去掉所有 hop shell）env_key 库。"
    )
    parser.add_argument(
        "--final-env-csv",
        type=Path,
        default=Path("outputs/final_env_keymap.csv"),
        help="输入 final_env_keymap.csv",
    )
    parser.add_argument(
        "--hop2-out",
        type=Path,
        default=Path("outputs/hop2_env_keymap.csv"),
        help="输出 hop2 库 CSV",
    )
    parser.add_argument(
        "--hop1-out",
        type=Path,
        default=Path("outputs/hop1_env_keymap.csv"),
        help="输出 hop1 库 CSV",
    )
    parser.add_argument(
        "--hop0-out",
        type=Path,
        default=Path("outputs/hop0_env_keymap.csv"),
        help="输出 hop0 兜底库 CSV（去掉所有 hop shell）",
    )
    args = parser.parse_args()

    if not args.final_env_csv.exists():
        raise FileNotFoundError(f"找不到输入数据库: {args.final_env_csv}")

    hop2_rows = build_hop_map(args.final_env_csv, hop=2)
    hop1_rows = build_hop_map(args.final_env_csv, hop=1)
    hop0_rows = build_hop_map(args.final_env_csv, hop=0)

    write_hop_csv(hop2_rows, args.hop2_out)
    write_hop_csv(hop1_rows, args.hop1_out)
    write_hop_csv(hop0_rows, args.hop0_out)

    print(f"[DONE] hop2: {args.hop2_out} ({len(hop2_rows)} keys)")
    print(f"[DONE] hop1: {args.hop1_out} ({len(hop1_rows)} keys)")
    print(f"[DONE] hop0: {args.hop0_out} ({len(hop0_rows)} keys)")


if __name__ == "__main__":
    main()
