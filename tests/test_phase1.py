"""Phase 1 acceptance tests: MolReader and MacroMapDB."""

from pathlib import Path

import pytest

from macromapff.db import MacroMapDB
from macromapff.mol import MolReader, computeHopKeys


# Base path for test fixtures (project root)
FIXTURES = Path(__file__).parent / "standard"


class TestMolReader:
    """Test MolReader can parse .mol/.pdb files."""

    @pytest.fixture
    def segment2_mol(self):
        return FIXTURES / "segdata" / "segment2" / "segment2.mol"

    def test_read_segment2_atoms(self, segment2_mol):
        """Can read atoms from segment2.mol."""
        reader = MolReader(segment2_mol)
        atoms = reader.getAtoms()
        assert len(atoms) == 55
        assert atoms[0]["symbol"] == "C"

    def test_read_segment2_bonds(self, segment2_mol):
        """Can read bonds from segment2.mol."""
        reader = MolReader(segment2_mol)
        bonds = reader.getBonds()
        assert len(bonds) == 56

    def test_read_segment2_coords(self, segment2_mol):
        """Can read 3D coordinates from segment2.mol."""
        reader = MolReader(segment2_mol)
        coords = reader.getCoords()
        assert len(coords) == 55
        # Check first atom has x, y, z
        x, y, z = coords[1]
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert isinstance(z, float)

    def test_hop0_env(self, segment2_mol):
        """Can compute hop0 environment for an atom."""
        reader = MolReader(segment2_mol)
        env = reader.computeHop0Env(1)
        assert "z" in env
        assert "neighbor_sig" in env
        assert "bond_kinds" in env

    def test_hop1_env(self, segment2_mol):
        """Can compute hop1 environment for an atom."""
        reader = MolReader(segment2_mol)
        env = reader.computeHop1Env(1)
        assert "hop1_shell" in env

    def test_hop2_env(self, segment2_mol):
        """Can compute hop2 environment for an atom."""
        reader = MolReader(segment2_mol)
        env = reader.computeHop2Env(1)
        assert "hop2_shell" in env

    def test_compute_hop_keys(self, segment2_mol):
        """Can compute hop2/hop1/hop0 keys for an atom."""
        reader = MolReader(segment2_mol)
        hop2Key, hop1Key, hop0Key = computeHopKeys(reader, 1)
        assert len(hop2Key) == 64
        assert len(hop1Key) == 64
        assert len(hop0Key) == 64
        # Keys should be different
        assert hop2Key != hop1Key
        assert hop1Key != hop0Key


class TestMacroMapDB:
    """Test MacroMapDB save/load round-trip."""

    def test_new_db_has_structure(self, tmp_path):
        """New empty db has correct structure."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        assert db._conn is not None
        assert len(db.atomTypes) == 0
        assert len(db.hop1Keymap) == 0
        assert len(db.hop0Keymap) == 0
        assert len(db.bondParams) == 0
        assert len(db.angleParams) == 0
        assert len(db.dihedralParams) == 0
        assert len(db.improperParams) == 0

    def test_insert_and_lookup_atom_type(self, tmp_path):
        """Can insert and retrieve atom type."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertAtomType("key123", {
            "element": "C",
            "hop1_key": "hop1key",
            "hop0_key": "hop0key",
            "hop2_env": {"z": 6, "neighbor_sig": "C2"},
            "mass": 12.011,
            "sigma": 3.5,
            "epsilon": 0.066,
            "source": ["seg1"],
        })
        result = db.getAtomType("key123")
        assert result is not None
        assert result["element"] == "C"
        assert result["hop1_key"] == "hop1key"
        assert result["hop2_env"]["z"] == 6

    def test_save_and_reload(self, tmp_path):
        """DB can be saved and reloaded."""
        db_path = tmp_path / "db.db"
        db = MacroMapDB(db_path)
        db.load()
        db.insertAtomType("key456", {
            "element": "O",
            "hop1_key": "hop1key2",
            "hop0_key": "hop0key2",
            "hop2_env": {"z": 8, "neighbor_sig": "C1"},
            "mass": 15.999,
            "sigma": 3.0,
            "epsilon": 0.12,
            "source": [],
        })
        db.save()

        db2 = MacroMapDB(db_path)
        db2.load()
        assert db2.getAtomType("key456")["element"] == "O"

    def test_hop1_keymap_insert(self, tmp_path):
        """Can insert hop1 key mapping."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertHop1Key("hop1key", 4)
        assert "hop1key" in db.hop1Keymap
        assert db.hop1Keymap["hop1key"]["lammps_types"] == [4]

    def test_hop0_keymap_insert(self, tmp_path):
        """Can insert hop0 key mapping."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertHop0Key("hop0key", 4)
        assert "hop0key" in db.hop0Keymap

    def test_bond_param_insert_and_lookup(self, tmp_path):
        """Can insert and lookup bond parameters."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertBondParam(("a", "b"), {"k": 340.0, "r0": 1.09})
        result = db.lookupBondParam("a", "b")
        assert result is not None
        assert result["k"] == 340.0

    def test_bond_param_canonical_order(self, tmp_path):
        """Bond param lookup is order-independent."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertBondParam(("a", "b"), {"k": 340.0, "r0": 1.09})
        assert db.lookupBondParam("a", "b") is not None
        assert db.lookupBondParam("b", "a") is not None

    def test_angle_param_canonical_order(self, tmp_path):
        """Angle param lookup is order-independent on outer atoms."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertAngleParam(("a", "b", "c"), {"k": 50.0, "theta0": 110.0})
        assert db.lookupAngleParam("a", "b", "c") is not None
        assert db.lookupAngleParam("c", "b", "a") is not None

    def test_dihedral_param_canonical_order(self, tmp_path):
        """Dihedral param lookup is order-independent on outer atoms."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertDihedralParam(("a", "b", "c", "d"), {"coeffs": [0, 0, 0.3, 0]})
        assert db.lookupDihedralParam("a", "b", "c", "d") is not None
        assert db.lookupDihedralParam("d", "c", "b", "a") is not None

    def test_improper_param_center_fixed(self, tmp_path):
        """Improper param lookup fixes center atom."""
        db = MacroMapDB(tmp_path / "db.db")
        db.load()
        db.insertImproperParam(("a", "b", "c", "d"), {"coeffs": [0, -1, 2]})
        # Center a is fixed, others sorted
        assert db.lookupImproperParam("a", "b", "c", "d") is not None
