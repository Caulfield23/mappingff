"""CLI entry point for mappingff.

This module provides the command-line interface. The actual workflow logic
is in workflow.py - this module only handles argument parsing and entry point.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from mappingff.workflow import build_db, parameterize


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="mappingff - Molecular force field parameterization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # build-db command
    build_parser = sub.add_parser(
        "build-db", help="Build parameter database from samples"
    )
    build_parser.add_argument(
        "samples_dir",
        type=Path,
        help="Directory containing sample subdirectories with .mol and .lmp files",
    )
    build_parser.add_argument(
        "-d",
        "--db",
        type=Path,
        default=Path("samples.db"),
        help="Output database file path (default: samples.db)",
    )
    build_parser.add_argument(
        "-a",
        "--append",
        action="store_true",
        help="Append to existing database instead of replacing (default: replace)",
    )

    # parameterize command
    param_parser = sub.add_parser("parameterize", help="Parameterize target molecule")
    param_parser.add_argument(
        "mol_file",
        type=Path,
        help="Path to target molecule .mol or .pdb file",
    )
    param_parser.add_argument(
        "-o",
        "--out",
        type=Path,
        help="Output LAMMPS file path (default: <mol_file>.lmp)",
    )
    param_parser.add_argument(
        "-d",
        "--db",
        type=Path,
        default=None,
        help="Path to database file (default: first .db file in mol_file's directory)",
    )
    param_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print detailed progress",
    )
    param_parser.add_argument(
        "-c",
        "--charge",
        type=float,
        default=None,
        help="Target total charge for the system (default: no adjustment)",
    )

    args = parser.parse_args()

    if args.command == "build-db":
        logging.basicConfig(
            level=logging.INFO,
            format="[%(levelname)s] %(name)s: %(message)s",
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(Path(args.db).parent / "build-db.log", mode="w"),
            ],
        )
        db_path = args.db
        build_db(args.samples_dir, db_path, append=args.append)
        print("Database build complete!")

    elif args.command == "parameterize":
        if args.out is None:
            args.out = args.mol_file.with_suffix(".lmp")

        if args.db is None:
            db_files = list(args.mol_file.parent.glob("*.db"))
            if db_files:
                args.db = db_files[0]

        if args.db is None or not args.db.exists():
            raise FileNotFoundError(f"Database not found: {args.db}")

        logging.basicConfig(
            level=logging.DEBUG if args.verbose else logging.INFO,
            format="[%(levelname)s] %(name)s: %(message)s",
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(
                    Path(args.out).parent / "parameterize.log", mode="w"
                ),
            ],
        )

        result = parameterize(
            args.mol_file,
            args.db,
            args.out,
            args.charge,
        )
        print(
            f"Parameterize complete: {result['atoms']} atoms, "
            f"{result['bonds']} bonds, {result['angles']} angles, "
            f"{result['dihedrals']} dihedrals, {result['impropers']} impropers, "
            f"{result['unique_types']} types, "
            f"hop3={result['hop3_matches']}, hop2={result['hop2_matches']}, "
            f"hop1={result['hop1_matches']}, hop0={result['hop0_matches']}, "
            f"no_match={result['no_match']}"
        )


if __name__ == "__main__":
    main()
