"""Database facade for MacroMapFF parameter database.

This module provides the MacroMapDB class which is a facade for all database
operations. The database is stored as a pickle file containing a dictionary with
the following top-level keys:
    - atom_types: hop2-level atom type definitions
    - hop1_keymap: hop1-level fallback mappings
    - hop0_keymap: hop0-level fallback mappings
    - bond_params: bond parameter definitions
    - angle_params: angle parameter definitions
    - dihedral_params: dihedral parameter definitions
    - improper_params: improper parameter definitions
    - meta: metadata (version, built_at, sample_count, source_segments)
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any


class MacroMapDB:
    """Database facade providing all database operations.

    This class manages a pickle-based database for storing molecular force field
    parameters. It supports loading, saving, and all CRUD operations for
    atom types, hop keymaps, and bonded parameters.

    The database uses SHA-256 encoded environment keys as primary keys for
    efficient dictionary-based lookup during parameterization.
    """

    def __init__(self, path: Path):
        """Initialize MacroMapDB with a database file path.

        Args:
            path: Path to the pickle database file.
        """
        self._path = Path(path)
        self._data: dict[str, Any] | None = None

    def load(self) -> None:
        """Load database from pickle file.

        If the file does not exist, creates a new empty database structure.
        The database is loaded entirely into memory for fast random access.
        """
        if self._path.exists():
            with open(self._path, "rb") as f:
                self._data = pickle.load(f)
        else:
            self._data = self._newEmptyDB()

    def save(self) -> None:
        """Save database to pickle file.

        Creates parent directories if they do not exist. The entire database
        is serialized to a pickle file in binary mode.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "wb") as f:
            pickle.dump(self._data, f)

    def _newEmptyDB(self) -> dict[str, Any]:
        """Create a new empty database structure.

        Returns:
            Dictionary with all required top-level keys initialized.
        """
        return {
            "version": 1,
            "meta": {
                "built_at": "",
                "sample_count": 0,
                "source_segments": [],
            },
            "atom_types": {},
            "hop1_keymap": {},
            "hop0_keymap": {},
            "bond_params": {},
            "angle_params": {},
            "dihedral_params": {},
            "improper_params": {},
        }

    # ── atom_types operations (hop2 level) ─────────────────────────────────────

    def insertAtomType(self, hop2Key: str, info: dict) -> None:
        """Insert or update an atom type entry at hop2 level.

        Args:
            hop2Key: SHA-256 key computed from hop2 environment.
            info: Dictionary containing atom type information.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        self._data["atom_types"][hop2Key] = info

    def getAtomType(self, hop2Key: str) -> dict | None:
        """Retrieve an atom type entry by hop2 key.

        Args:
            hop2Key: SHA-256 key computed from hop2 environment.

        Returns:
            Atom type info dict if found, None otherwise.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["atom_types"].get(hop2Key)

    # ── hop1_keymap operations ─────────────────────────────────────────────────

    def insertHop1Key(self, hop1Key: str, hop0Key: str, lammpsType: int) -> None:
        """Insert or update a hop1_keymap entry.

        The hop1_keymap provides hop1-level fallback for atom typing.
        Multiple hop2 entries may map to the same hop1 key, so lammps_types
        is stored as a list to track all possible types.

        Args:
            hop1Key: SHA-256 key computed from hop1 environment.
            hop0Key: Associated hop0 key for bonded parameter lookup.
            lammpsType: LAMMPS atom type ID from a sample.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if hop1Key in self._data["hop1_keymap"]:
            entry = self._data["hop1_keymap"][hop1Key]
            entry["lammps_types"].append(lammpsType)
            entry["lammps_types"].sort()
        else:
            self._data["hop1_keymap"][hop1Key] = {
                "hop0_key": hop0Key,
                "lammps_types": [lammpsType],
            }

    # ── hop0_keymap operations ────────────────────────────────────────────────

    def insertHop0Key(self, hop0Key: str, lammpsType: int) -> None:
        """Insert or update a hop0_keymap entry.

        The hop0_keymap provides hop0-level fallback for atom typing.
        Multiple hop2/hop1 entries may map to the same hop0 key, so
        lammps_types is stored as a list.

        Args:
            hop0Key: SHA-256 key computed from hop0 environment.
            lammpsType: LAMMPS atom type ID from a sample.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if hop0Key in self._data["hop0_keymap"]:
            self._data["hop0_keymap"][hop0Key]["lammps_types"].append(lammpsType)
            self._data["hop0_keymap"][hop0Key]["lammps_types"].sort()
        else:
            self._data["hop0_keymap"][hop0Key] = {
                "lammps_types": [lammpsType],
            }

    # ── Bonded parameter operations ───────────────────────────────────────────

    def insertBondParam(self, key: tuple, params: dict) -> None:
        """Insert or merge a bond parameter entry.

        When the same key already exists (same atom type pair), parameters
        are averaged with the existing values.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB) in lexicographic order.
            params: Dict with 'k' (force constant) and 'r0' (equilibrium distance).
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if key in self._data["bond_params"]:
            existing = self._data["bond_params"][key]
            existing["k"] = (existing["k"] + params["k"]) / 2.0
            existing["r0"] = (existing["r0"] + params["r0"]) / 2.0
        else:
            self._data["bond_params"][key] = params.copy()

    def lookupBondParam(self, hop0KeyA: str, hop0KeyB: str) -> dict | None:
        """Look up bond parameter by hop0 key pair.

        Keys are stored in lexicographic order, so the input order does not matter.

        Args:
            hop0KeyA: hop0 key of first atom.
            hop0KeyB: hop0 key of second atom.

        Returns:
            Bond parameter dict if found, None otherwise.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if hop0KeyA <= hop0KeyB:
            key = (hop0KeyA, hop0KeyB)
        else:
            key = (hop0KeyB, hop0KeyA)
        return self._data["bond_params"].get(key)

    def insertAngleParam(self, key: tuple, params: dict) -> None:
        """Insert or merge an angle parameter entry.

        When the same key already exists, parameters are averaged.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB, hop0KeyC) where B is the center atom.
                 Stored in lexicographic order of outer atoms (A and C can be swapped).
            params: Dict with 'k' (force constant) and 'theta0' (equilibrium angle).
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if key in self._data["angle_params"]:
            existing = self._data["angle_params"][key]
            existing["k"] = (existing["k"] + params["k"]) / 2.0
            existing["theta0"] = (existing["theta0"] + params["theta0"]) / 2.0
        else:
            self._data["angle_params"][key] = params.copy()

    def lookupAngleParam(self, hop0KeyA: str, hop0KeyB: str, hop0KeyC: str) -> dict | None:
        """Look up angle parameter by hop0 key triple.

        The center atom (B) is fixed, outer atoms (A and C) can be swapped.

        Args:
            hop0KeyA: hop0 key of first atom.
            hop0KeyB: hop0 key of center atom.
            hop0KeyC: hop0 key of third atom.

        Returns:
            Angle parameter dict if found, None otherwise.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if hop0KeyA <= hop0KeyC:
            key = (hop0KeyA, hop0KeyB, hop0KeyC)
        else:
            key = (hop0KeyC, hop0KeyB, hop0KeyA)
        return self._data["angle_params"].get(key)

    def insertDihedralParam(self, key: tuple, params: dict) -> None:
        """Insert or merge a dihedral parameter entry.

        When the same key already exists, coefficient arrays are averaged.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD) in canonical order.
                 Outer atoms (A and D) can be swapped.
            params: Dict with 'coeffs' list (e.g., [0, 0, 0.3, 0] for OPLS dihedral).
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if key in self._data["dihedral_params"]:
            existing = self._data["dihedral_params"][key]
            coeffs = params.get("coeffs", [])
            existingCoeffs = existing.get("coeffs", [])
            if len(coeffs) == len(existingCoeffs):
                existing["coeffs"] = [
                    (e + c) / 2.0 for e, c in zip(existingCoeffs, coeffs)
                ]
        else:
            self._data["dihedral_params"][key] = params.copy()

    def lookupDihedralParam(
        self, hop0KeyA: str, hop0KeyB: str, hop0KeyC: str, hop0KeyD: str
    ) -> dict | None:
        """Look up dihedral parameter by hop0 key quadruple.

        Outer atoms (A and D) can be swapped while preserving the order of B and C.

        Args:
            hop0KeyA: hop0 key of first atom.
            hop0KeyB: hop0 key of second atom.
            hop0KeyC: hop0 key of third atom.
            hop0KeyD: hop0 key of fourth atom.

        Returns:
            Dihedral parameter dict if found, None otherwise.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if hop0KeyA <= hop0KeyD:
            key = (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD)
        else:
            key = (hop0KeyD, hop0KeyC, hop0KeyB, hop0KeyA)
        return self._data["dihedral_params"].get(key)

    def insertImproperParam(self, key: tuple, params: dict) -> None:
        """Insert or merge an improper parameter entry.

        The first atom in the key is the center atom and is fixed.
        The remaining three atoms are sorted for canonical ordering.

        Args:
            key: Tuple of (hop0KeyCenter, hop0KeyA, hop0KeyB, hop0KeyC).
            params: Dict with 'coeffs' list (e.g., [0, -1, 2] for harmonic improper).
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        if key in self._data["improper_params"]:
            existing = self._data["improper_params"][key]
            coeffs = params.get("coeffs", [])
            existingCoeffs = existing.get("coeffs", [])
            if len(coeffs) == len(existingCoeffs):
                existing["coeffs"] = [
                    (e + c) / 2.0 for e, c in zip(existingCoeffs, coeffs)
                ]
        else:
            self._data["improper_params"][key] = params.copy()

    def lookupImproperParam(
        self, hop0KeyA: str, hop0KeyB: str, hop0KeyC: str, hop0KeyD: str
    ) -> dict | None:
        """Look up improper parameter by hop0 key quadruple.

        The first atom (A) is the center atom and is fixed.
        The remaining three atoms are sorted for canonical ordering.

        Args:
            hop0KeyA: hop0 key of center atom (fixed position).
            hop0KeyB: hop0 key of second atom.
            hop0KeyC: hop0 key of third atom.
            hop0KeyD: hop0 key of fourth atom.

        Returns:
            Improper parameter dict if found, None otherwise.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        others = sorted([hop0KeyB, hop0KeyC, hop0KeyD])
        key = (hop0KeyA, others[0], others[1], others[2])
        return self._data["improper_params"].get(key)

    # ── Merge and export operations ────────────────────────────────────────────

    def mergeSample(self, sampleData: dict) -> None:
        """Merge data from a sample into the database.

        This is used during add-samples to incrementally update the database
        with new samples without rebuilding from scratch.

        Args:
            sampleData: Dictionary containing atom_types, hop1_keymap, hop0_keymap,
                       and other parameter tables to merge.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")

        for hop2Key, info in sampleData.get("atom_types", {}).items():
            if hop2Key in self._data["atom_types"]:
                existing = self._data["atom_types"][hop2Key]
                existing["source"] = list(set(existing["source"]) | set(info["source"]))
            else:
                self._data["atom_types"][hop2Key] = info.copy()

        for hop1Key, entry in sampleData.get("hop1_keymap", {}).items():
            if hop1Key in self._data["hop1_keymap"]:
                existing = self._data["hop1_keymap"][hop1Key]
                existing["lammps_types"] = list(
                    set(existing["lammps_types"]) | set(entry["lammps_types"])
                )
                existing["lammps_types"].sort()
            else:
                self._data["hop1_keymap"][hop1Key] = entry.copy()

        for hop0Key, entry in sampleData.get("hop0_keymap", {}).items():
            if hop0Key in self._data["hop0_keymap"]:
                existing = self._data["hop0_keymap"][hop0Key]
                existing["lammps_types"] = list(
                    set(existing["lammps_types"]) | set(entry["lammps_types"])
                )
                existing["lammps_types"].sort()
            else:
                self._data["hop0_keymap"][hop0Key] = entry.copy()

        self._data["meta"]["sample_count"] = len(
            set(
                seg
                for sources in self._data["atom_types"].values()
                for seg in sources.get("source", [])
            )
        )

    def export(self) -> dict:
        """Export the database as a plain dictionary.

        Returns a shallow copy of the internal data dictionary.
        The caller can modify it without affecting the internal state.

        Returns:
            Copy of the database dictionary.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data.copy()

    # ── Property accessors ─────────────────────────────────────────────────────

    @property
    def atomTypes(self) -> dict:
        """Get the atom_types table (hop2 level).

        Returns:
            Dictionary mapping hop2 keys to atom type info.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["atom_types"]

    @property
    def hop1Keymap(self) -> dict:
        """Get the hop1_keymap table.

        Returns:
            Dictionary mapping hop1 keys to {hop0_key, lammps_types}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["hop1_keymap"]

    @property
    def hop0Keymap(self) -> dict:
        """Get the hop0_keymap table.

        Returns:
            Dictionary mapping hop0 keys to {lammps_types}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["hop0_keymap"]

    @property
    def bondParams(self) -> dict:
        """Get the bond_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB) to {k, r0}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["bond_params"]

    @property
    def angleParams(self) -> dict:
        """Get the angle_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB, hop0KeyC) to {k, theta0}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["angle_params"]

    @property
    def dihedralParams(self) -> dict:
        """Get the dihedral_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD) to {coeffs}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["dihedral_params"]

    @property
    def improperParams(self) -> dict:
        """Get the improper_params table.

        Returns:
            Dictionary mapping (hop0KeyCenter, hop0KeyA, hop0KeyB, hop0KeyC) to {coeffs}.
        """
        if self._data is None:
            raise RuntimeError("Database not loaded")
        return self._data["improper_params"]
