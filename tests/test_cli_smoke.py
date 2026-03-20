import subprocess
import sys


def test_cli_help() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "macromapff.cli", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "MacroMapFF molecular parameterization pipeline CLI" in result.stdout
