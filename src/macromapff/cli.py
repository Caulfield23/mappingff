from __future__ import annotations

import argparse
from pathlib import Path

from macromapff.workflow import add_samples, build_database, parameterize_molecule


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="MacroMapFF",
        description="MacroMapFF simplified CLI",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_build = sub.add_parser(
        "build-db",
        help="Scan a folder of sample .lammps.lmp files and build merged databases.",
    )
    p_build.add_argument("samples", type=Path, help="Folder containing sample data")
    p_build.add_argument(
        "--db-dir",
        type=Path,
        default=Path("database"),
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
        default=Path("database"),
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
        default=Path("database"),
        help="Database folder (default: ./database)",
    )

    args = parser.parse_args()

    if args.command == "build-db":
        build_database(samples_root=args.samples, db_dir=args.db_dir)
        return
    if args.command == "add-samples":
        add_samples(samples_root=args.samples, db_dir=args.db_dir)
        return
    if args.command == "parameterize":
        parameterize_molecule(mol_path=args.mol, db_dir=args.db_dir, out_path=args.out)
        return

    raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
