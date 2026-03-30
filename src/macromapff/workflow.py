from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List


def _pipeline_script(name: str) -> Path:
    return Path(__file__).resolve().parent / "pipeline" / name


def _run(cmd: List[str]) -> None:
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")


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


def _manifest_path(db_dir: Path) -> Path:
    return db_dir / "samples_manifest.json"


def _write_manifest(db_dir: Path, samples: List[Dict[str, str]]) -> None:
    path = _manifest_path(db_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"samples": samples}, indent=2), encoding="utf-8")


def _load_manifest(db_dir: Path) -> List[Dict[str, str]]:
    path = _manifest_path(db_dir)
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return list(data.get("samples", []))


def _build_from_samples(samples: List[Dict[str, str]], db_dir: Path) -> None:
    db_dir = db_dir.expanduser().resolve()
    db_dir.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    envkey_script = _pipeline_script("env_build.py")
    extract_script = _pipeline_script("multi_extract.py")
    final_script = _pipeline_script("keymap_build.py")
    hop_script = _pipeline_script("hop_build.py")
    master_script = _pipeline_script("multi_build.py")

    final_specs: List[str] = []
    multi_specs: List[str] = []

    for s in samples:
        module = s["module"]
        mol = Path(s["mol"])
        lmp = Path(s["lmp"])

        env_out = db_dir / f"{module}_envdb"
        atom_env_csv = env_out / f"{module}_atom_env.csv"
        multi_prefix = db_dir / f"{module}_multiatom_observed"

        # Internal defaults are intentionally fixed to keep user interface minimal.
        _run(
            [
                python,
                str(envkey_script),
                "--mol",
                str(mol),
                "--lmp",
                str(lmp),
                "--module",
                module,
                "--outdir",
                str(env_out),
            ]
        )

        _run(
            [
                python,
                str(extract_script),
                "--lmp",
                str(lmp),
                "--atom-env-csv",
                str(atom_env_csv),
                "--out-prefix",
                str(multi_prefix),
            ]
        )

        final_specs.append(f"{module}::{atom_env_csv}::{lmp}")
        multi_specs.append(f"{module}::{multi_prefix}.csv")

    final_prefix = db_dir / "final_env_keymap"
    cmd_final = [python, str(final_script), "--out-prefix", str(final_prefix)]
    for spec in final_specs:
        cmd_final.extend(["--module-spec", spec])
    _run(cmd_final)

    _run(
        [
            python,
            str(hop_script),
            "--final-env-csv",
            str(db_dir / "final_env_keymap.csv"),
            "--hop2-out",
            str(db_dir / "hop2_env_keymap.csv"),
            "--hop1-out",
            str(db_dir / "hop1_env_keymap.csv"),
            "--hop0-out",
            str(db_dir / "hop0_env_keymap.csv"),
        ]
    )

    cmd_master = [
        python,
        str(master_script),
        "--type-stats-csv",
        str(db_dir / "final_env_keymap_type_stats.csv"),
        "--hop0-env-csv",
        str(db_dir / "hop0_env_keymap.csv"),
        "--out-prefix",
        str(db_dir / "multiatom_master_keytype"),
    ]
    for spec in multi_specs:
        cmd_master.extend(["--multiatom-spec", spec])
    _run(cmd_master)


def build_database(samples_root: Path, db_dir: Path) -> None:
    samples = _discover_samples(samples_root)
    _build_from_samples(samples=samples, db_dir=db_dir)
    _write_manifest(db_dir=db_dir, samples=samples)
    print(f"[DONE] built database from {len(samples)} samples -> {db_dir}")


def add_samples(samples_root: Path, db_dir: Path) -> None:
    existing = _load_manifest(db_dir)
    new = _discover_samples(samples_root)

    seen_lmp = {entry["lmp"] for entry in existing}
    merged = list(existing)

    for entry in new:
        if entry["lmp"] not in seen_lmp:
            merged.append(entry)
            seen_lmp.add(entry["lmp"])

    if not merged:
        raise ValueError("No samples available to build database.")

    _build_from_samples(samples=merged, db_dir=db_dir)
    _write_manifest(db_dir=db_dir, samples=merged)
    print(f"[DONE] merged sample count: {len(merged)} -> {db_dir}")


def parameterize_molecule(mol_path: Path, db_dir: Path, out_path: Path | None = None) -> None:
    mol_path = mol_path.expanduser().resolve()
    db_dir = db_dir.expanduser().resolve()

    if out_path is None:
        out_path = mol_path.with_name(f"{mol_path.stem}_param.lmp")
    out_path = out_path.expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    generate_script = _pipeline_script("lammps_gen.py")

    _run(
        [
            python,
            str(generate_script),
            "--structure",
            str(mol_path),
            "--db-dir",
            str(db_dir),
            "--out",
            str(out_path),
        ]
    )

    print(f"[DONE] parameterized output: {out_path}")
