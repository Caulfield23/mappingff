"""Top-level workflow orchestration for build-db, add-samples and parameterize commands."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from macromapff.pipeline.atommap_sample import build_sample_atommap
from macromapff.pipeline.bondedterms_sample import extract_bondedterms_sample
from macromapff.pipeline.global_bonded import build_global_bonded
from macromapff.pipeline.keymap_hop import build_hop_databases
from macromapff.pipeline.keymap_hop import build_keymap
from macromapff.pipeline.parameterize import parameterize_lammps


# User input default.
USER_DEFAULT_DB_DIR = Path("database")

# Project layout constants.
PROJECT_HOP_DIR = "hop_env"
PROJECT_SAMPLE_ENV_SUFFIX = "_env"
PROJECT_GLOBAL_ATOMMAP_PREFIX = "Global_AtomMap"
PROJECT_GLOBAL_BONDED_PREFIX = "Global_BondedTerms"
PROJECT_BUILD_LOG = "build.log"

# Internal constants.
INTERNAL_STRUCTURE_PATTERNS = ("*.mol", "*.mol2", "*.pdb")


def discover_samples(samples_root: Path) -> list[tuple[str, Path, Path]]:
    """Discover sample triplets as (module, structure_path, lammps_path)."""
    root = samples_root.expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise FileNotFoundError(f"Samples folder not found: {root}")

    lmp_files = sorted(root.rglob("*.lmp"))
    if not lmp_files:
        raise ValueError(f"No .lmp files found under: {root}")

    structure_by_dir: dict[Path, list[Path]] = defaultdict(list)
    for pattern in INTERNAL_STRUCTURE_PATTERNS:
        for structure in root.rglob(pattern):
            structure_by_dir[structure.parent].append(structure)

    discovered: list[tuple[str, Path, Path]] = []
    used_module_ids: set[str] = set()

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

        discovered.append((module_id, structure_path.resolve(), lmp.resolve()))

    return discovered


def build_sample_envs(db_dir: Path, samples: list[tuple[str, Path, Path]]) -> None:
    """Build/overwrite sample env folders under db_dir for all discovered samples."""
    db_dir.mkdir(parents=True, exist_ok=True)
    (db_dir / PROJECT_HOP_DIR).mkdir(parents=True, exist_ok=True)

    for module, mol, lmp in samples:
        env_out = db_dir / f"{module}{PROJECT_SAMPLE_ENV_SUFFIX}"
        atom_env_csv = env_out / f"{module}_AtomMap.csv"
        multi_prefix = env_out / f"{module}_BondedTerms"

        build_sample_atommap(
            structure_path=mol,
            out_dir=env_out,
            module=module,
            lmp_path=lmp,
        )

        extract_bondedterms_sample(
            lmp_path=lmp,
            atom_env_csv=atom_env_csv,
            out_prefix=multi_prefix,
        )


def discover_sample_env_records(db_dir: Path) -> list[tuple[str, Path, Path]]:
    """Scan db_dir and return merge records as (module, atom_env_csv, bondedterms_csv)."""
    records: list[tuple[str, Path, Path]] = []
    for env_dir in sorted(db_dir.glob(f"*{PROJECT_SAMPLE_ENV_SUFFIX}")):
        if not env_dir.is_dir() or env_dir.name == PROJECT_HOP_DIR:
            continue
        module = env_dir.name[: -len(PROJECT_SAMPLE_ENV_SUFFIX)]
        atom_env_csv = env_dir / f"{module}_AtomMap.csv"
        bondedterms_csv = env_dir / f"{module}_BondedTerms.csv"
        records.append((module, atom_env_csv.resolve(), bondedterms_csv.resolve()))
    return records


def merge_database(db_dir: Path, records: list[tuple[str, Path, Path]]) -> int:
    """Merge all discovered sample env folders into global databases."""
    if not records:
        raise ValueError("No sample records available to merge.")

    final_prefix = db_dir / PROJECT_GLOBAL_ATOMMAP_PREFIX
    build_log = db_dir / PROJECT_BUILD_LOG
    hop_dir = db_dir / PROJECT_HOP_DIR
    hop_dir.mkdir(parents=True, exist_ok=True)

    final_csv, _, _, _ = build_keymap(
        module_specs=[(m, atom) for m, atom, _ in records],
        out_prefix=final_prefix,
        out_log=build_log,
    )

    hop2_out = hop_dir / "hop2_KeyMap.csv"
    hop1_out = hop_dir / "hop1_KeyMap.csv"
    hop0_out = hop_dir / "hop0_KeyMap.csv"
    build_hop_databases(
        final_env_csv=final_csv,
        hop2_out=hop2_out,
        hop1_out=hop1_out,
        hop0_out=hop0_out,
    )

    build_global_bonded(
        final_env_csv=final_csv,
        hop0_env_csv=hop0_out,
        bonded_specs=[(m, bonded) for m, _, bonded in records],
        out_prefix=db_dir / PROJECT_GLOBAL_BONDED_PREFIX,
        log_file=build_log,
    )
    return len(records)


def build_db(samples_root: Path, db_dir: Path) -> dict:
    """Build sample env folders then merge the whole database."""
    resolved_db_dir = db_dir.expanduser().resolve()
    samples = discover_samples(samples_root)
    build_sample_envs(resolved_db_dir, samples)
    merged_count = merge_database(
        resolved_db_dir,
        discover_sample_env_records(resolved_db_dir),
    )
    return {
        "samples_count": len(samples),
        "merged_count": merged_count,
        "db_dir": resolved_db_dir,
    }


def add_samples(samples_root: Path, db_dir: Path) -> dict:
    """Add/overwrite sample env folders then re-merge the whole database."""
    return build_db(samples_root, db_dir)


def parameterize(mol_path: Path, db_dir: Path, out_path: Path | None) -> dict:
    """Generate a parameterized LAMMPS data file for one molecule."""
    mol = mol_path.expanduser().resolve()
    out = (
        out_path.expanduser().resolve()
        if out_path is not None
        else mol.with_name(f"{mol.stem}_param.lmp")
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    parameterize_lammps(db_dir=db_dir, structure=mol, out=out)
    return {
        "out": out,
    }
