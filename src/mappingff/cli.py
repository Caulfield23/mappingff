"""CLI entry point for mappingff.

This module provides the command-line interface. The actual workflow logic
is in workflow.py - this module only handles argument parsing and entry point.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from mappingff.utils import USER_DEFAULT_DB_PATH, setupLogging
from mappingff.workflow import buildDb, parameterize


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="mappingff - Molecular force field parameterization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # build-db command
    build_parser = sub.add_parser("build-db", help="Build parameter database from samples")
    build_parser.add_argument(
        "samples_dir",
        type=Path,
        help="Directory containing sample subdirectories with .mol and .lmp files",
    )
    build_parser.add_argument(
        "-d", "--db",
        type=Path,
        default=USER_DEFAULT_DB_PATH,
        help="Output database file path (default: samples.db)",
    )
    build_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress",
    )

    # parameterize command
    param_parser = sub.add_parser("parameterize", help="Parameterize target molecule")
    param_parser.add_argument(
        "mol_file",
        type=Path,
        help="Path to target molecule .mol or .pdb file",
    )
    param_parser.add_argument(
        "-o", "--out",
        type=Path,
        help="Output LAMMPS file path (default: <mol_file>.lmp)",
    )
    param_parser.add_argument(
        "-d", "--db",
        type=Path,
        default=USER_DEFAULT_DB_PATH,
        help="Path to database file (default: samples.db)",
    )
    param_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress",
    )
    param_parser.add_argument(
        "-c", "--charge",
        type=float,
        default=None,
        help="Target total charge for the system (default: no adjustment)",
    )

    args = parser.parse_args()

    # Setup logging
    logLevel = logging.DEBUG if getattr(args, "verbose", False) else logging.INFO
    setupLogging(logLevel)

    if args.command == "build-db":
        dbPath = args.db
        result = buildDb(args.samples_dir, dbPath, args.verbose)
        print(f"Build complete: {result['samples_count']} samples, {result['atoms_processed']} atoms")

    elif args.command == "parameterize":
        result = parameterize(
            args.mol_file,
            args.db,
            args.out,
            args.verbose,
            args.charge,
        )
        print(f"Parameterize complete: {result['atoms']} atoms, "
              f"{result['bonds']} bonds, "
              f"{result['angles']} angles, "
              f"{result['dihedrals']} dihedrals, "
              f"{result['impropers']} impropers, "
              f"{result['unique_types']} types, "
              f"hop3={result['hop3_matches']}, "
              f"hop2={result['hop2_matches']}, "
              f"hop1={result['hop1_matches']}, "
              f"hop0={result['hop0_matches']}, "
              f"no_match={result['no_match']}")


if __name__ == "__main__":
    main()
