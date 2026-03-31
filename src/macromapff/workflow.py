from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from macromapff.pipeline.env_build import build_mapping, INTERNAL_HOP_DEPTH
from macromapff.pipeline.multi_extract import extract_multiatom_mapping
from macromapff.pipeline.keymap_build import build_final_keymap
from macromapff.pipeline.hop_build import build_hop_databases
from macromapff.pipeline.multi_build import build_multiatom_master, parse_multiatom_spec
from macromapff.pipeline.lammps_gen import generate_lammps_data

def _discover_samples(samples_root: Path) -> List[Dict[str, str]]:
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


def _build_from_samples(samples: List[Dict[str, str]], db_dir: Path) -> None:
    db_dir = db_dir.expanduser().resolve()
    db_dir.mkdir(parents=True, exist_ok=True)

    hop_dir = db_dir / "hop_env"
    hop_dir.mkdir(parents=True, exist_ok=True)

    final_specs: List[str] = []
    multi_specs: List[str] = []

    for s in samples:
        module = s["module"]
        mol = Path(s["mol"])
        lmp = Path(s["lmp"])

        env_out = db_dir / f"{module}_envdb"
        atom_env_csv = env_out / f"{module}_atom_env.csv"
        multi_prefix = env_out / f"{module}_multiatom_observed"

        # invoke env_build
        build_mapping(
            structure_path=mol,
            out_dir=env_out,
            module=module,
            lmp_path=lmp,
            hop_depth=INTERNAL_HOP_DEPTH,
        )

        # invoke multi_extract
        multi_csv, _, _ = extract_multiatom_mapping(
            lmp_path=lmp,
            atom_env_csv=atom_env_csv,
            out_prefix=multi_prefix,
        )

        final_specs.append(f"{module}::{atom_env_csv}::{lmp}")
        multi_specs.append(f"{module}::{multi_prefix}.csv")

    final_prefix = db_dir / "final_env_keymap"
    external_final_log = db_dir.parent / "final_env_keymap.log"

    # invoke keymap_build
    final_csv, _, _, _ = build_final_keymap(
        module_specs=final_specs,
        out_prefix=final_prefix,
        out_log=external_final_log,
    )

    # invoke hop_build
    hop2_out = hop_dir / "hop2_env_keymap.csv"
    hop1_out = hop_dir / "hop1_env_keymap.csv"
    hop0_out = hop_dir / "hop0_env_keymap.csv"
    build_hop_databases(
        final_env_csv=final_csv,
        hop2_out=hop2_out,
        hop1_out=hop1_out,
        hop0_out=hop0_out,
    )

    # invoke multi_build
    master_out_prefix = db_dir / "multiatom_master_keytype"
    multi_specs_parsed = [parse_multiatom_spec(s) for s in multi_specs]
    build_multiatom_master(
        final_env_csv=final_csv,
        hop0_env_csv=hop0_out,
        multiatom_specs=multi_specs_parsed,
        out_prefix=master_out_prefix,
        log_file=external_final_log,
    )

def build_database(samples_root: Path, db_dir: Path) -> None:
    samples = _discover_samples(samples_root)
    _build_from_samples(samples=samples, db_dir=db_dir)
    print(f"[DONE] built database from {len(samples)} samples -> {db_dir}")

def add_samples(samples_root: Path, db_dir: Path) -> None:
    samples = _discover_samples(samples_root)
    if not samples:
        raise ValueError("No samples available to build database.")

    _build_from_samples(samples=samples, db_dir=db_dir)
    print(f"[DONE] rebuilt database from {len(samples)} samples -> {db_dir}")

def parameterize_molecule(mol_path: Path, db_dir: Path, out_path: Path | None = None) -> None:
    mol_path = mol_path.expanduser().resolve()
    db_dir = db_dir.expanduser().resolve()

    if out_path is None:
        out_path = mol_path.with_name(f"{mol_path.stem}_param.lmp")
    out_path = out_path.expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    generate_lammps_data(
        structure=mol_path,
        db_dir=db_dir,
        out=out_path,
    )

    print(f"[DONE] parameterized output: {out_path}")
