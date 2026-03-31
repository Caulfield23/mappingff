from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict
from typing import List

from macromapff.pipeline import HopDatabaseBuilder
from macromapff.pipeline import INTERNAL_HOP_DEPTH
from macromapff.pipeline import KeymapBuilder
from macromapff.pipeline import LammpsGenerator
from macromapff.pipeline import MultiatomExtractor
from macromapff.pipeline import MultiatomMasterBuilder
from macromapff.pipeline import build_mapping
from macromapff.pipeline import parse_multiatom_spec


class Workflow:
    """Orchestrates database building and molecule parameterization workflows."""

    def __init__(self, db_dir: Path) -> None:
        """Initialize workflow state with normalized database paths."""
        self.db_dir = db_dir.expanduser().resolve()
        self.hop_dir = self.db_dir / "hop_env"

    def discover_samples(self, samples_root: Path) -> List[Dict[str, str]]:
        """Discover sample pairs of structure and LAMMPS data files."""
        root = samples_root.expanduser().resolve()
        if not root.exists() or not root.is_dir():
            raise FileNotFoundError(f"Samples folder not found: {root}")

        lmp_files = sorted(root.rglob("*.lammps.lmp"))
        if not lmp_files:
            raise ValueError(f"No .lammps.lmp files found under: {root}")

        structure_by_stem: Dict[str, List[Path]] = {}
        for pattern in ("*.mol", "*.mol2", "*.pdb"):
            for structure in root.rglob(pattern):
                structure_by_stem.setdefault(structure.stem, []).append(structure)

        discovered: List[Dict[str, str]] = []
        used_module_ids = set()

        for lmp in lmp_files:
            stem = lmp.name.replace(".lammps.lmp", "")

            candidates = []
            for ext in (".mol", ".mol2", ".pdb"):
                candidates.extend(
                    [
                        lmp.parent / f"{stem}{ext}",
                        lmp.parent.parent / f"{stem}{ext}",
                        lmp.parent.parent.parent / f"{stem}{ext}",
                    ]
                )
            structure_path = next((p for p in candidates if p.exists()), None)
            if structure_path is None and stem in structure_by_stem:
                structure_path = sorted(structure_by_stem[stem], key=lambda p: len(str(p)))[0]
            if structure_path is None:
                raise FileNotFoundError(
                    f"Cannot find matching structure (.mol/.mol2/.pdb) for {lmp}"
                )

            module_id = stem
            idx = 2
            while module_id in used_module_ids:
                module_id = f"{stem}_{idx}"
                idx += 1
            used_module_ids.add(module_id)

            discovered.append(
                {
                    "module": module_id,
                    "mol": str(structure_path.resolve()),
                    "lmp": str(lmp.resolve()),
                }
            )

        return discovered

    def build_database(self, samples_root: Path) -> None:
        """Build all mapping databases from discovered training samples."""
        samples = self.discover_samples(samples_root)
        self._build_from_samples(samples)
        print(f"[DONE] built database from {len(samples)} samples -> {self.db_dir}")

    def add_samples(self, samples_root: Path) -> None:
        """Rebuild databases using existing and newly provided samples."""
        samples = self.discover_samples(samples_root)
        if not samples:
            raise ValueError("No samples available to build database.")
        self._build_from_samples(samples)
        print(f"[DONE] rebuilt database from {len(samples)} samples -> {self.db_dir}")

    def parameterize_molecule(self, mol_path: Path, out_path: Path | None = None) -> None:
        """Generate a parameterized LAMMPS data file for one molecule."""
        mol_path = mol_path.expanduser().resolve()

        if out_path is None:
            out_path = mol_path.with_name(f"{mol_path.stem}_param.lmp")
        out_path = out_path.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        generator = LammpsGenerator(db_dir=self.db_dir)
        generator.generate(structure=mol_path, out=out_path)

        print(f"[DONE] parameterized output: {out_path}")

    def _build_from_samples(self, samples: List[Dict[str, str]]) -> None:
        """Execute the full multi-stage build pipeline over sample metadata."""
        self.db_dir.mkdir(parents=True, exist_ok=True)
        self.hop_dir.mkdir(parents=True, exist_ok=True)

        final_specs: List[str] = []
        multi_specs: List[str] = []

        for s in samples:
            module = s["module"]
            mol = Path(s["mol"])
            lmp = Path(s["lmp"])

            env_out = self.db_dir / f"{module}_envdb"
            atom_env_csv = env_out / f"{module}_atom_env.csv"
            multi_prefix = env_out / f"{module}_multiatom_observed"

            build_mapping(
                structure_path=mol,
                out_dir=env_out,
                module=module,
                lmp_path=lmp,
                hop_depth=INTERNAL_HOP_DEPTH,
            )

            extractor = MultiatomExtractor(lmp_path=lmp, atom_env_csv=atom_env_csv)
            extractor.extract(out_prefix=multi_prefix)

            final_specs.append(f"{module}::{atom_env_csv}::{lmp}")
            multi_specs.append(f"{module}::{multi_prefix}.csv")

        final_prefix = self.db_dir / "final_env_keymap"
        external_final_log = self.db_dir.parent / "final_env_keymap.log"

        keymap_builder = KeymapBuilder(final_specs)
        final_csv, _, _, _ = keymap_builder.build(
            out_prefix=final_prefix,
            out_log=external_final_log,
        )

        hop2_out = self.hop_dir / "hop2_env_keymap.csv"
        hop1_out = self.hop_dir / "hop1_env_keymap.csv"
        hop0_out = self.hop_dir / "hop0_env_keymap.csv"
        hop_builder = HopDatabaseBuilder(final_env_csv=final_csv)
        hop_builder.build(
            hop2_out=hop2_out,
            hop1_out=hop1_out,
            hop0_out=hop0_out,
        )

        master_out_prefix = self.db_dir / "multiatom_master_keytype"
        multi_specs_parsed = [parse_multiatom_spec(s) for s in multi_specs]
        multi_builder = MultiatomMasterBuilder(
            final_env_csv=final_csv,
            hop0_env_csv=hop0_out,
        )
        multi_builder.build(
            multiatom_specs=multi_specs_parsed,
            out_prefix=master_out_prefix,
            log_file=external_final_log,
        )


def main() -> None:
    """Parse CLI arguments and dispatch to workflow subcommands."""
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
        Workflow(db_dir=args.db_dir).build_database(args.samples)
        return
    if args.command == "add-samples":
        Workflow(db_dir=args.db_dir).add_samples(args.samples)
        return
    if args.command == "parameterize":
        Workflow(db_dir=args.db_dir).parameterize_molecule(args.mol, out_path=args.out)
        return

    raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
