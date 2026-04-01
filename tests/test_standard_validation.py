from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "standard"
SEGDATA_DIR = FIXTURE_ROOT / "segdata"
TARGET_MOL = FIXTURE_ROOT / "target" / "PS-oDMS7POSS.mol"
ARTIFACT_ROOT = Path(__file__).resolve().parent / "artifacts" / "standard"


def _env_with_src() -> dict:
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(repo_root / "src") + (
        os.pathsep + existing if existing else ""
    )
    return env


def _prepare_case_dir(case_name: str) -> Path:
    case_dir = ARTIFACT_ROOT / case_name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    return case_dir


def test_standard_workflow_with_fixed_dataset() -> None:
    assert SEGDATA_DIR.exists(), f"Missing standard segdata folder: {SEGDATA_DIR}"
    assert TARGET_MOL.exists(), f"Missing standard molecule: {TARGET_MOL}"

    case_dir = _prepare_case_dir("case_standard_workflow")
    db_dir = case_dir / "database"
    out_lmp = case_dir / "PS-oDMS7POSS_param.lmp"

    build = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "build-db",
            str(SEGDATA_DIR),
            "--db-dir",
            str(db_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert build.returncode == 0, build.stdout + "\n" + build.stderr

    expected_db_files = [
        db_dir / "Global_AtomMap.csv",
        db_dir / "hop_env" / "hop2_KeyMap.csv",
        db_dir / "hop_env" / "hop1_KeyMap.csv",
        db_dir / "hop_env" / "hop0_KeyMap.csv",
        db_dir / "Global_BondedTerms.csv",
    ]
    for path in expected_db_files:
        assert path.exists(), f"Expected database artifact missing: {path}"

    atommap_header = (db_dir / "segment1_env" / "segment1_AtomMap.csv").read_text(
        encoding="utf-8", errors="ignore"
    ).splitlines()[0]
    assert "mass" in atommap_header.split(","), "AtomMap CSV must include mass column"

    param = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "parameterize",
            str(TARGET_MOL),
            "--db-dir",
            str(db_dir),
            "--out",
            str(out_lmp),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert param.returncode == 0, param.stdout + "\n" + param.stderr
    assert out_lmp.exists(), "Parameterized LAMMPS output was not generated"

    head = out_lmp.read_text(encoding="utf-8", errors="ignore").splitlines()[:5]
    assert any("LAMMPS data file" in line for line in head), (
        "Output does not look like a LAMMPS data file"
    )


def test_add_samples_keeps_existing_sample_env_and_merges() -> None:
    case_dir = _prepare_case_dir("case_add_samples_merge")
    db_dir = case_dir / "database"
    seg1_root = SEGDATA_DIR / "segment1"
    seg2_root = SEGDATA_DIR / "segment2"

    build = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "build-db",
            str(seg1_root),
            "--db-dir",
            str(db_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert build.returncode == 0, build.stdout + "\n" + build.stderr
    assert (db_dir / "segment1_env" / "segment1_AtomMap.csv").exists()

    add = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "add-samples",
            str(seg2_root),
            "--db-dir",
            str(db_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert add.returncode == 0, add.stdout + "\n" + add.stderr

    assert (db_dir / "segment1_env" / "segment1_AtomMap.csv").exists()
    assert (db_dir / "segment2_env" / "segment2_AtomMap.csv").exists()
    assert (db_dir / "Global_AtomMap.csv").exists()
    assert (db_dir / "Global_BondedTerms.csv").exists()


def test_add_samples_overwrites_on_module_name_conflict() -> None:
    case_dir = _prepare_case_dir("case_add_samples_overwrite")
    db_dir = case_dir / "database"
    seg1_root = SEGDATA_DIR / "segment1"

    build = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "build-db",
            str(seg1_root),
            "--db-dir",
            str(db_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert build.returncode == 0, build.stdout + "\n" + build.stderr

    atommap = db_dir / "segment1_env" / "segment1_AtomMap.csv"
    original = atommap.read_text(encoding="utf-8", errors="ignore")
    marker = "# MANUAL_MARKER_SHOULD_BE_OVERWRITTEN\n"
    atommap.write_text(original + marker, encoding="utf-8")

    add = subprocess.run(
        [
            sys.executable,
            "-m",
            "macromapff.cli",
            "add-samples",
            str(seg1_root),
            "--db-dir",
            str(db_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env_with_src(),
    )
    assert add.returncode == 0, add.stdout + "\n" + add.stderr

    overwritten = atommap.read_text(encoding="utf-8", errors="ignore")
    assert marker not in overwritten
    assert (db_dir / "Global_AtomMap.csv").exists()
    assert (db_dir / "Global_BondedTerms.csv").exists()
