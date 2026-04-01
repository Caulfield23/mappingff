"""Command-line entrypoint for MacroMapFF workflows."""

from __future__ import annotations

import argparse
from pathlib import Path

from macromapff.pipeline import USER_DEFAULT_DB_DIR
from macromapff.pipeline import add_samples
from macromapff.pipeline import build_db
from macromapff.pipeline import parameterize


def main() -> None:
    """Parse CLI arguments and dispatch to workflow subcommands."""
    parser = argparse.ArgumentParser(
        prog="MacroMapFF",
        description="MacroMapFF simplified CLI",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_build = sub.add_parser(
        "build-db",
        help="Scan a folder of sample .lmp files and build merged databases.",
    )
    p_build.add_argument("samples", type=Path, help="Folder containing sample data")
    p_build.add_argument(
        "--db-dir",
        type=Path,
        default=USER_DEFAULT_DB_DIR,
        help="Database output folder (default: ./database)",
    )

    p_add = sub.add_parser(
        "add-samples",
        help="Append new samples and rebuild merged databases.",
    )
    p_add.add_argument("samples", type=Path, help="Folder containing new sample data")
    p_add.add_argument(
        "--db-dir",
        type=Path,
        default=USER_DEFAULT_DB_DIR,
        help="Existing database folder (default: ./database)",
    )

    p_param = sub.add_parser(
        "parameterize",
        help="Generate a parameterized LAMMPS data file for a new molecule.",
    )
    p_param.add_argument("mol", type=Path, help="Input molecule .mol file")
    p_param.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output LAMMPS data file (default: <mol_stem>_param.lmp)",
    )
    p_param.add_argument(
        "--db-dir",
        type=Path,
        default=USER_DEFAULT_DB_DIR,
        help="Database folder (default: ./database)",
    )

    args = parser.parse_args()
    db_dir = args.db_dir

    if args.command == "build-db":
        result = build_db(args.samples, db_dir)
        print(
            f"[DONE] built database from {result['samples_count']} samples, "
            f"merged total {result['merged_count']} samples -> {result['db_dir']}"
        )
        return
    if args.command == "add-samples":
        result = add_samples(args.samples, db_dir)
        print(
            f"[DONE] added {result['samples_count']} samples (overwrite on name conflict), "
            f"merged total {result['merged_count']} samples -> {result['db_dir']}"
        )
        return
    if args.command == "parameterize":
        result = parameterize(args.mol, db_dir, args.out)
        print(f"[DONE] parameterized output: {result['out']}")
        return

    raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
