"""Edge-case tests for mappingff behavior."""

from pathlib import Path

import pytest

import mappingff.workflow as workflow
from mappingff.mol import MolReader

FIXTURES = Path(__file__).parent / "standard"
SEG_DATA = FIXTURES / "samples"
TARGET_FILE = FIXTURES / "50_ps_50_pmma.mol"


def test_parameterize_skips_charge_adjustment_when_none(tmp_path, monkeypatch):
    db_path = tmp_path / "samples.db"
    out_path = tmp_path / "target.lmp"

    workflow.build_db(SEG_DATA, db_path)

    def fail_if_called():
        raise AssertionError("adjust_total_charge should not be called")

    monkeypatch.setattr(workflow, "adjust_total_charge", fail_if_called)

    workflow.parameterize(TARGET_FILE, db_path, out_path, total_charge=None)
    assert out_path.exists()


def test_molreader_rejects_unsupported_suffix(tmp_path):
    bad_path = tmp_path / "bad_input.txt"
    bad_path.write_text("not a molecule")

    with pytest.raises(ValueError, match="Unsupported file type"):
        MolReader(bad_path)
