import subprocess
import sys
from pathlib import Path
import os


def test_cli_help() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(repo_root / "src") + (
        os.pathsep + existing if existing else ""
    )

    result = subprocess.run(
        [sys.executable, "-m", "macromapff.cli", "--help"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0
    assert "MacroMapFF simplified CLI" in result.stdout
