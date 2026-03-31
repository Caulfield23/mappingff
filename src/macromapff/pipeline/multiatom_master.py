#!/usr/bin/env python3
"""Build a merged multi-atom master keymap from observed module data."""

from pathlib import Path

from macromapff.domain import build_master
from macromapff.io import append_merge_conflict_log
from macromapff.io import load_hop0_key_classes
from macromapff.io import load_type_to_keyid
from macromapff.io import write_master_csv


def parse_multiatom_spec(spec: str):
    parts = spec.split("::")
    if len(parts) != 2:
        raise ValueError(
            f"Invalid --multiatom-spec format: {spec}\nExpected: module_name::multiatom_csv"
        )
    module_name, csv_path = parts
    return module_name.strip(), Path(csv_path).expanduser()


class MultiatomMasterBuilder:
    def __init__(self, final_env_csv: Path, hop0_env_csv: Path) -> None:
        self.final_env_csv = final_env_csv
        self.hop0_env_csv = hop0_env_csv

    def build(self, multiatom_specs, out_prefix: Path, log_file: Path | None = None):
        type_to_keyid = load_type_to_keyid(self.final_env_csv)
        key_to_class = load_hop0_key_classes(self.hop0_env_csv)
        rows, missing_type_refs = build_master(type_to_keyid, key_to_class, multiatom_specs)

        out_prefix = Path(out_prefix).expanduser()
        out_csv = out_prefix.with_suffix(".csv")
        write_master_csv(out_csv, rows)

        if log_file is not None:
            append_merge_conflict_log(Path(log_file).expanduser(), rows)

        return out_csv, rows, missing_type_refs
