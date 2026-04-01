#!/usr/bin/env python3
"""Build a merged global bonded keymap from observed module data."""

from pathlib import Path

from macromapff.domain.bonded_global_merge import build_global_bonded_map
from macromapff.io.input import load_hop0_key_classes
from macromapff.io.input import load_type_to_keyid
from macromapff.io.log import append_merge_conflict_log
from macromapff.io.output import write_global_bonded_csv


def build_global_bonded(
    final_env_csv: Path,
    hop0_env_csv: Path,
    bonded_specs,
    out_prefix: Path,
    log_file: Path,
):
    """Merge module observations and export global bonded CSV (plus log)."""
    type_to_keyid = load_type_to_keyid(final_env_csv)
    key_to_class = load_hop0_key_classes(hop0_env_csv)
    rows, missing_type_refs = build_global_bonded_map(
        type_to_keyid,
        key_to_class,
        bonded_specs,
    )

    out_prefix = Path(out_prefix).expanduser()
    out_csv = out_prefix.with_suffix(".csv")
    write_global_bonded_csv(out_csv, rows)

    append_merge_conflict_log(Path(log_file).expanduser(), rows)

    return out_csv, rows, missing_type_refs
