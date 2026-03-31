from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict
from typing import List

from macromapff.pipeline import BondedTermsSampleExtractor
from macromapff.pipeline import HopDatabaseBuilder
from macromapff.pipeline import INTERNAL_HOP_DEPTH
from macromapff.pipeline import KeymapBuilder
from macromapff.pipeline import LammpsGenerator
from macromapff.pipeline import MultiatomMasterBuilder
from macromapff.pipeline import build_sample_atommap
from macromapff.pipeline import parse_multiatom_spec


class Workflow:
    """Orchestrates database building and molecule parameterization workflows."""

    def __init__(self, db_dir: Path) -> None:
        """Initialize workflow state with normalized database paths."""
        self.db_dir = db_dir.expanduser().resolve()
        self.hop_dir = self.db_dir / "hop_env"
        self.manifest_csv = self.db_dir / "samples_manifest.csv"

    def discover_samples(self, samples_root: Path) -> List[Dict[str, str]]:
        """Discover sample pairs of structure and LAMMPS data files."""
        root = samples_root.expanduser().resolve()
        if not root.exists() or not root.is_dir():
            raise FileNotFoundError(f"Samples folder not found: {root}")

        lmp_files = sorted(root.rglob("*.lmp"))
        if not lmp_files:
            raise ValueError(f"No .lmp files found under: {root}")

        structure_by_dir: Dict[Path, List[Path]] = defaultdict(list)
        for pattern in ("*.mol", "*.mol2", "*.pdb"):
            for structure in root.rglob(pattern):
                structure_by_dir[structure.parent].append(structure)

        discovered: List[Dict[str, str]] = []
        used_module_ids = set()

        for lmp in lmp_files:
            search_dirs = [lmp.parent, lmp.parent.parent, lmp.parent.parent.parent]
            structure_path = None
            for folder in search_dirs:
                choices = sorted(structure_by_dir.get(folder, []))
                if choices:
                    structure_path = choices[0]
                    break
            if structure_path is None:
                raise FileNotFoundError(
                    f"Cannot find structure file (.mol/.mol2/.pdb) for {lmp}"
                )

            module_seed = lmp.parent.name or lmp.stem
            module_id = module_seed
            idx = 2
            while module_id in used_module_ids:
                module_id = f"{module_seed}_{idx}"
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
        records = self._build_new_sample_records(samples=samples, used_module_ids=set())
        self._write_manifest(records)
        self._merge_from_records(records)
        print(f"[DONE] built database from {len(samples)} samples -> {self.db_dir}")

    def add_samples(self, samples_root: Path) -> None:
        """Add new sample mappings and re-merge global databases."""
        existing = self._read_manifest()
        used_module_ids = {row["module"] for row in existing}

        samples = self.discover_samples(samples_root)
        if not samples:
            raise ValueError("No new samples found to add.")

        new_records = self._build_new_sample_records(
            samples=samples,
            used_module_ids=used_module_ids,
        )
        merged = existing + new_records
        self._write_manifest(merged)
        self._merge_from_records(merged)
        print(
            f"[DONE] added {len(new_records)} new samples, merged total {len(merged)} samples -> {self.db_dir}"
        )

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

    def _build_new_sample_records(
        self,
        samples: List[Dict[str, str]],
        used_module_ids: set[str],
    ) -> List[Dict[str, str]]:
        """Build per-sample artifacts for new samples and return manifest rows."""
        self.db_dir.mkdir(parents=True, exist_ok=True)
        self.hop_dir.mkdir(parents=True, exist_ok=True)

        records: List[Dict[str, str]] = []

        for s in samples:
            module_seed = s["module"]
            module = module_seed
            suffix = 2
            while module in used_module_ids:
                module = f"{module_seed}_{suffix}"
                suffix += 1
            used_module_ids.add(module)

            mol = Path(s["mol"])
            lmp = Path(s["lmp"])

            env_out = self.db_dir / f"{module}_env"
            atom_env_csv = env_out / f"{module}_AtomMap.csv"
            multi_prefix = env_out / f"{module}_BondedTerms"
            multi_csv = multi_prefix.with_suffix(".csv")

            build_sample_atommap(
                structure_path=mol,
                out_dir=env_out,
                module=module,
                lmp_path=lmp,
                hop_depth=INTERNAL_HOP_DEPTH,
            )

            extractor = BondedTermsSampleExtractor(lmp_path=lmp, atom_env_csv=atom_env_csv)
            extractor.extract(out_prefix=multi_prefix)

            records.append(
                {
                    "module": module,
                    "mol": str(mol.resolve()),
                    "lmp": str(lmp.resolve()),
                    "atom_env_csv": str(atom_env_csv.resolve()),
                    "bondedterms_csv": str(multi_csv.resolve()),
                }
            )

        return records

    def _merge_from_records(self, records: List[Dict[str, str]]) -> None:
        """Merge global databases from manifest rows."""
        if not records:
            raise ValueError("No sample records available to merge.")

        final_specs: List[str] = []
        multi_specs: List[str] = []
        for row in records:
            atom_env_csv = Path(row["atom_env_csv"]).expanduser().resolve()
            bondedterms_csv = Path(row["bondedterms_csv"]).expanduser().resolve()
            if not atom_env_csv.exists():
                raise FileNotFoundError(f"AtomMap CSV not found: {atom_env_csv}")
            if not bondedterms_csv.exists():
                raise FileNotFoundError(f"BondedTerms CSV not found: {bondedterms_csv}")
            final_specs.append(f"{row['module']}::{atom_env_csv}")
            multi_specs.append(f"{row['module']}::{bondedterms_csv}")

        final_prefix = self.db_dir / "Global_AtomMap"
        build_log = self.db_dir / "build.log"

        keymap_builder = KeymapBuilder(final_specs)
        final_csv, _, _, _ = keymap_builder.build(
            out_prefix=final_prefix,
            out_log=build_log,
        )

        hop2_out = self.hop_dir / "hop2_KeyMap.csv"
        hop1_out = self.hop_dir / "hop1_KeyMap.csv"
        hop0_out = self.hop_dir / "hop0_KeyMap.csv"
        hop_builder = HopDatabaseBuilder(final_env_csv=final_csv)
        hop_builder.build(
            hop2_out=hop2_out,
            hop1_out=hop1_out,
            hop0_out=hop0_out,
        )

        master_out_prefix = self.db_dir / "Global_BondedTerms"
        multi_specs_parsed = [parse_multiatom_spec(s) for s in multi_specs]
        multi_builder = MultiatomMasterBuilder(
            final_env_csv=final_csv,
            hop0_env_csv=hop0_out,
        )
        multi_builder.build(
            multiatom_specs=multi_specs_parsed,
            out_prefix=master_out_prefix,
            log_file=build_log,
        )

    def _read_manifest(self) -> List[Dict[str, str]]:
        """Load sample manifest rows from database directory."""
        if not self.manifest_csv.exists():
            return []

        rows: List[Dict[str, str]] = []
        with self.manifest_csv.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            required = {"module", "mol", "lmp", "atom_env_csv", "bondedterms_csv"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Invalid samples manifest, missing columns: {sorted(missing)} ({self.manifest_csv})"
                )
            for row in reader:
                if not row.get("module"):
                    continue
                rows.append(
                    {
                        "module": str(row["module"]),
                        "mol": str(row["mol"]),
                        "lmp": str(row["lmp"]),
                        "atom_env_csv": str(row["atom_env_csv"]),
                        "bondedterms_csv": str(row["bondedterms_csv"]),
                    }
                )
        return rows

    def _write_manifest(self, rows: List[Dict[str, str]]) -> None:
        """Write sample manifest rows into database directory."""
        self.manifest_csv.parent.mkdir(parents=True, exist_ok=True)
        with self.manifest_csv.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["module", "mol", "lmp", "atom_env_csv", "bondedterms_csv"],
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        "module": row["module"],
                        "mol": row["mol"],
                        "lmp": row["lmp"],
                        "atom_env_csv": row["atom_env_csv"],
                        "bondedterms_csv": row["bondedterms_csv"],
                    }
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
        help="Scan a folder of sample .lmp files and build merged databases.",
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
