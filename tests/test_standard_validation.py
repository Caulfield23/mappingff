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


def test_standard_workflow_with_fixed_dataset() -> None:
    assert SEGDATA_DIR.exists(), f"Missing standard segdata folder: {SEGDATA_DIR}"
    assert TARGET_MOL.exists(), f"Missing standard molecule: {TARGET_MOL}"

    if ARTIFACT_ROOT.exists():
        shutil.rmtree(ARTIFACT_ROOT)
    ARTIFACT_ROOT.mkdir(parents=True, exist_ok=True)

    db_dir = ARTIFACT_ROOT / "database"
    out_lmp = ARTIFACT_ROOT / "PS-oDMS7POSS_param.lmp"

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
        db_dir / "final_env_keymap.csv",
        db_dir / "hop2_env_keymap.csv",
        db_dir / "hop1_env_keymap.csv",
        db_dir / "hop0_env_keymap.csv",
        db_dir / "multiatom_master_keytype.csv",
        db_dir / "samples_manifest.json",
    ]
    for path in expected_db_files:
        assert path.exists(), f"Expected database artifact missing: {path}"

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
