#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

from macromapff.domain import INTERACTION_ORDER


def append_build_log(log_path: Path, lines):
    """Append plain-text lines to build log file if enabled."""
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


def init_build_log(log_path: Path, lines):
    """Initialize or overwrite build log with header lines."""
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
    """Write structured report for atoms missing env-key matches."""
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


def log_multiatom_match(
    kind: str,
    missing,
    ambiguous,
    cache_hit: int,
    cache_miss: int,
    cache_size: int,
    build_log_path: Path = None,
):
    """Log multi-atom matching misses/ambiguities and cache statistics."""
    if missing:
        sample = ", ".join(f"atom_ids={x[0]} candidates={x[1]}" for x in missing[:10])
        append_build_log(
            build_log_path,
            [
                f"[{kind.upper()}][WARN] missing_count={len(missing)}",
                f"[{kind.upper()}][WARN] sample_missing={sample}",
            ],
        )

    if ambiguous:
        lines = [f"[{kind.upper()}][AMBIGUOUS] count={len(ambiguous)}"]
        for item in ambiguous[:200]:
            lines.append(
                f"[{kind.upper()}][AMBIGUOUS] atom_ids={item['atom_ids']} "
                f"chosen_key_tuple={item['chosen_key_tuple']} chosen_coeff={item['chosen_coeff']} "
                f"all_key_tuples={item['matched_key_tuples']} all_coeffs={item['matched_coeffs']}"
            )
        append_build_log(build_log_path, lines)

    total_cache_req = cache_hit + cache_miss
    hit_rate = (cache_hit / total_cache_req) if total_cache_req else 0.0
    append_build_log(
        build_log_path,
        [
            f"[{kind.upper()}][CACHE] hit={cache_hit} miss={cache_miss} size={cache_size} hit_rate={hit_rate:.4f}"
        ],
    )


def _fmt_float(v):
    """Format optional float values for human-readable logs."""
    if v is None:
        return ""
    return f"{v:.8f}"


def write_keymap_log(path: Path, module_specs, final_rows, type_rows):
    """Write detailed final keymap merge diagnostics report."""
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append("=== Final keymap merge log ===")
    lines.append(f"modules: {len(module_specs)}")
    for module_name, atom_env_csv in module_specs:
        lines.append(f"  - {module_name}: atom_env={atom_env_csv}")

    lines.append("")
    lines.append(f"total keys: {len(final_rows)}")

    multi_type_keys = [r for r in final_rows if r["n_lmp_types"] > 1]
    lines.append(f"keys with multiple LMP types: {len(multi_type_keys)}")

    lines.append("")
    lines.append("notes: types = number of distinct LMP atom types merged into this key")
    lines.append("notes: rows = number of atom rows merged for this key across all modules")
    lines.append("--- key-level stats (mean ± std) ---")
    for row in final_rows:
        lines.append(
            f"key_id={row['key_id']} "
            f"types={row['n_lmp_types']} rows={row['n_rows']} "
            f"charge={_fmt_float(row['charge_mean'])}±{_fmt_float(row['charge_std'])} "
            f"sigma={_fmt_float(row['sigma_mean'])}±{_fmt_float(row['sigma_std'])} "
            f"epsilon={_fmt_float(row['epsilon_mean'])}±{_fmt_float(row['epsilon_std'])} "
            f"mass={_fmt_float(row['mass_mean'])}±{_fmt_float(row['mass_std'])}"
        )

    lines.append("")
    lines.append("notes: multi-type key details is atom-type merge detail only (not bond/angle/dihedral/improper)")
    lines.append("--- multi-type key details ---")

    type_rows_by_key = defaultdict(list)
    for tr in type_rows:
        type_rows_by_key[tr["key_id"]].append(tr)

    for row in multi_type_keys:
        lines.append(f"key_id={row['key_id']} lmp_type_ids=[{row['lmp_type_ids']}]")
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


def append_merge_conflict_log(log_path: Path, rows):
    """Append grouped multi-atom coeff conflict summary to log file."""
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