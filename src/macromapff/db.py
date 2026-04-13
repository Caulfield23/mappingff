"""Database facade for MacroMapFF parameter database.

This module provides the MacroMapDB class which is a facade for all database
operations. The database is stored as an SQLite file with the following tables:
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

import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any


class MacroMapDB:
    """Database facade providing all database operations.

    This class manages an SQLite-based database for storing molecular force field
    parameters. It supports loading, saving, and all CRUD operations for
    atom types, hop keymaps, and bonded parameters.
    """

    def __init__(self, path: Path):
        """Initialize MacroMapDB with a database file path.

        Args:
            path: Path to the SQLite database file.
        """
        self._path = Path(path)
        self._conn: sqlite3.Connection | None = None

    def load(self) -> None:
        """Load database from SQLite file.

        If the file does not exist, creates a new database with all tables.
        """
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self._path)
        self._conn.row_factory = sqlite3.Row
        self._initSchema()

    def save(self) -> None:
        """Save database to SQLite file.

        Finalizes atom type and bonded parameter lists by computing weighted averages,
        then commits any pending transactions.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")
        self._finalizeAtomTypes()
        self._finalizeBondedParams()
        self._conn.commit()

    def _finalizeAtomTypes(self) -> None:
        """Compute averaged sigma/epsilon from lists and update the scalar columns.

        sigma is rounded to 7 decimal places, epsilon to 3 decimal places.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")
        cursor = self._conn.cursor()
        cursor.execute("SELECT hop2_key, sigma_list, epsilon_list, charge_list FROM atom_types")
        for row in cursor.fetchall():
            sigma_list = json.loads(row["sigma_list"])
            epsilon_list = json.loads(row["epsilon_list"])
            charge_list = json.loads(row["charge_list"])
            avg_sigma = round(sum(sigma_list) / len(sigma_list), 7)
            avg_epsilon = round(sum(epsilon_list) / len(epsilon_list), 3)
            avg_charge = round(sum(charge_list) / len(charge_list), 6)
            cursor.execute("""
                UPDATE atom_types
                SET sigma = ?, epsilon = ?, charge = ?
                WHERE hop2_key = ?
            """, (avg_sigma, avg_epsilon, avg_charge, row["hop2_key"]))

    def _finalizeBondedParams(self) -> None:
        """Compute averaged bonded parameters from lists and update scalar columns."""
        if self._conn is None:
            raise RuntimeError("Database not loaded")
        cursor = self._conn.cursor()

        # Bond params
        cursor.execute("SELECT rowid, coeffs_list FROM bond_params")
        for row in cursor.fetchall():
            coeffs_list = json.loads(row["coeffs_list"])
            if not coeffs_list:
                continue
            num_terms = len(coeffs_list[0])
            averaged = []
            for i in range(num_terms):
                vals = [c[i] for c in coeffs_list]
                if i < 2:
                    avg_val = round(sum(vals) / len(vals), 4)
                else:
                    avg_val = round(sum(vals) / len(vals), 3)
                averaged.append(avg_val)
            cursor.execute("UPDATE bond_params SET coeffs = ? WHERE rowid = ?",
                           (json.dumps(averaged), row["rowid"]))

        # Angle params
        cursor.execute("SELECT rowid, coeffs_list FROM angle_params")
        for row in cursor.fetchall():
            coeffs_list = json.loads(row["coeffs_list"])
            if not coeffs_list:
                continue
            num_terms = len(coeffs_list[0])
            averaged = []
            for i in range(num_terms):
                vals = [c[i] for c in coeffs_list]
                avg_val = round(sum(vals) / len(vals), 3)
                averaged.append(avg_val)
            cursor.execute("UPDATE angle_params SET coeffs = ? WHERE rowid = ?",
                           (json.dumps(averaged), row["rowid"]))

        # Dihedral params
        cursor.execute("SELECT rowid, coeffs_list FROM dihedral_params")
        for row in cursor.fetchall():
            coeffs_list = json.loads(row["coeffs_list"])
            if not coeffs_list:
                continue
            # coeffs_list is list of lists, each inner list is one sample's coeffs
            num_terms = len(coeffs_list[0])
            averaged = []
            for i in range(num_terms):
                vals = [c[i] for c in coeffs_list]
                avg_val = round(sum(vals) / len(vals), 3)
                averaged.append(avg_val)
            cursor.execute("UPDATE dihedral_params SET coeffs = ? WHERE rowid = ?",
                           (json.dumps(averaged), row["rowid"]))

        # Improper params: last two fields are categorical indices (mode, not average)
        cursor.execute("SELECT rowid, coeffs_list FROM improper_params")
        for row in cursor.fetchall():
            coeffs_list = json.loads(row["coeffs_list"])
            if not coeffs_list:
                continue

            # 1. Count occurrences of (last_two) pairs as a whole
            last_two_counts = Counter()
            for c in coeffs_list:
                key = (c[-2], c[-1])  # treat as combined key
                last_two_counts[key] += 1

            # 2. Get most common pair (if tie, most_common returns first by insertion order)
            max_count = max(last_two_counts.values())
            most_common_pairs = [k for k, v in last_two_counts.items() if v == max_count]
            mode_last_two = most_common_pairs[0]

            # 3. Filter records sharing the mode
            filtered = [c for c in coeffs_list if (c[-2], c[-1]) == mode_last_two]

            # 4. Average first term among filtered records
            first_vals = [c[0] for c in filtered]
            avg_first = round(sum(first_vals) / len(first_vals), 3)

            # 4. Result = [avg_first] + mode_last_two
            averaged = [avg_first] + list(mode_last_two)
            cursor.execute("UPDATE improper_params SET coeffs = ? WHERE rowid = ?",
                           (json.dumps(averaged), row["rowid"]))

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def _initSchema(self) -> None:
        """Create database tables if they don't exist."""
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()

        # Meta table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS meta (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        """)

        # Atom types table (hop2 level)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS atom_types (
                lammps_type INTEGER,
                element TEXT,
                mass REAL,
                sigma REAL,
                sigma_list TEXT,
                epsilon REAL,
                epsilon_list TEXT,
                charge REAL,
                charge_list TEXT,
                source TEXT,
                hop2_key TEXT PRIMARY KEY,
                hop1_key TEXT,
                hop0_key TEXT,
                hop2_env TEXT
            )
        """)

        # hop1 keymap (for external validation only)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS hop1_keymap (
                hop1_key TEXT PRIMARY KEY,
                lammps_types TEXT
            )
        """)

        # hop0 keymap
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS hop0_keymap (
                hop0_key TEXT PRIMARY KEY,
                lammps_types TEXT
            )
        """)

        # Bond params
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS bond_params (
                hop0_key_a TEXT,
                hop0_key_b TEXT,
                coeffs TEXT,
                coeffs_list TEXT,
                PRIMARY KEY (hop0_key_a, hop0_key_b)
            )
        """)

        # Angle params
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS angle_params (
                hop0_key_a TEXT,
                hop0_key_b TEXT,
                hop0_key_c TEXT,
                coeffs TEXT,
                coeffs_list TEXT,
                PRIMARY KEY (hop0_key_a, hop0_key_b, hop0_key_c)
            )
        """)

        # Dihedral params
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS dihedral_params (
                hop0_key_a TEXT,
                hop0_key_b TEXT,
                hop0_key_c TEXT,
                hop0_key_d TEXT,
                coeffs TEXT,
                coeffs_list TEXT,
                PRIMARY KEY (hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d)
            )
        """)

        # Improper params
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS improper_params (
                hop0_key_a TEXT,
                hop0_key_b TEXT,
                hop0_key_c TEXT,
                hop0_key_d TEXT,
                coeffs TEXT,
                coeffs_list TEXT,
                PRIMARY KEY (hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d)
            )
        """)

        self._conn.commit()

    # ── atom_types operations (hop2 level) ─────────────────────────────────────

    def insertAtomType(self, hop2Key: str, info: dict) -> None:
        """Insert or update an atom type entry at hop2 level.

        If hop2Key already exists, merge sigma_list, epsilon_list, and source list.
        If hop2Key is new, assign a unique lammps_type.

        Args:
            hop2Key: SHA-256 key computed from hop2 environment.
            info: Dictionary containing element, hop1_key, hop0_key, hop2_env,
                  mass, sigma, epsilon, source.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()

        # Check if hop2Key already exists
        cursor.execute("SELECT * FROM atom_types WHERE hop2_key = ?", (hop2Key,))
        existing = cursor.fetchone()

        if existing is not None:
            # hop2Key exists - merge sigma_list, epsilon_list, charge_list, and source
            existing_sigma = json.loads(existing["sigma_list"])
            existing_epsilon = json.loads(existing["epsilon_list"])
            existing_charge = json.loads(existing["charge_list"])
            existing_sources = json.loads(existing["source"])

            new_sigma = existing_sigma + [info["sigma"]]
            new_epsilon = existing_epsilon + [info["epsilon"]]
            new_charge = existing_charge + [info["charge"]]
            new_sources = list(set(existing_sources) | set(info["source"]))

            cursor.execute("""
                UPDATE atom_types
                SET epsilon_list = ?, sigma_list = ?, charge_list = ?, source = ?
                WHERE hop2_key = ?
            """, (json.dumps(new_epsilon), json.dumps(new_sigma), json.dumps(new_charge), json.dumps(new_sources), hop2Key))
        else:
            # hop2Key is new - get next available lammps_type
            cursor.execute("SELECT MAX(lammps_type) FROM atom_types")
            row = cursor.fetchone()
            max_type = row[0] if row[0] is not None else 0
            new_lammps_type = max_type + 1

            cursor.execute("""
                INSERT INTO atom_types
                (hop2_key, element, hop1_key, hop0_key, lammps_type,
                 mass, sigma, epsilon, sigma_list, epsilon_list,
                 charge, charge_list, source, hop2_env)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                hop2Key,
                info["element"],
                info["hop1_key"],
                info["hop0_key"],
                new_lammps_type,
                round(info["mass"], 3),
                round(info["sigma"], 7),
                round(info["epsilon"], 3),
                json.dumps([round(info["sigma"], 7)]),
                json.dumps([round(info["epsilon"], 3)]),
                round(info.get("charge", 0.0), 6),
                json.dumps([round(info.get("charge", 0.0), 6)]),
                json.dumps(info["source"]),
                json.dumps(info.get("hop2_env", {})),
            ))

    def getAtomType(self, hop2Key: str) -> dict | None:
        """Retrieve an atom type entry by hop2 key.

        Args:
            hop2Key: SHA-256 key computed from hop2 environment.

        Returns:
            Atom type info dict if found, None otherwise.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("SELECT * FROM atom_types WHERE hop2_key = ?", (hop2Key,))
        row = cursor.fetchone()
        if row is None:
            return None

        return {
            "element": row["element"],
            "hop1_key": row["hop1_key"],
            "hop0_key": row["hop0_key"],
            "lammps_type": row["lammps_type"],
            "mass": row["mass"],
            "sigma": row["sigma"],
            "epsilon": row["epsilon"],
            "charge": row["charge"],
            "sigma_list": json.loads(row["sigma_list"]),
            "epsilon_list": json.loads(row["epsilon_list"]),
            "charge_list": json.loads(row["charge_list"]),
            "source": json.loads(row["source"]),
            "hop2_env": json.loads(row["hop2_env"]) if row["hop2_env"] else {},
        }

    # ── hop1_keymap operations ─────────────────────────────────────────────────

    def insertHop1Key(self, hop1Key: str, lammpsType: int) -> None:
        """Insert or update a hop1_keymap entry.

        The hop1_keymap provides hop1-level fallback for atom typing.
        Multiple hop2 entries may map to the same hop1 key, so lammps_types
        is stored as a list to track all possible types.

        Args:
            hop1Key: SHA-256 key computed from hop1 environment.
            lammpsType: LAMMPS atom type ID.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("SELECT lammps_types FROM hop1_keymap WHERE hop1_key = ?", (hop1Key,))
        row = cursor.fetchone()

        if row is not None:
            lammps_types = json.loads(row["lammps_types"])
            if lammpsType not in lammps_types:
                lammps_types.append(lammpsType)
                lammps_types.sort()
            cursor.execute("""
                UPDATE hop1_keymap SET lammps_types = ? WHERE hop1_key = ?
            """, (json.dumps(lammps_types), hop1Key))
        else:
            cursor.execute("""
                INSERT INTO hop1_keymap (hop1_key, lammps_types)
                VALUES (?, ?)
            """, (hop1Key, json.dumps([lammpsType])))

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
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("SELECT lammps_types FROM hop0_keymap WHERE hop0_key = ?", (hop0Key,))
        row = cursor.fetchone()

        if row is not None:
            lammps_types = json.loads(row["lammps_types"])
            if lammpsType not in lammps_types:
                lammps_types.append(lammpsType)
                lammps_types.sort()
            cursor.execute("""
                UPDATE hop0_keymap SET lammps_types = ? WHERE hop0_key = ?
            """, (json.dumps(lammps_types), hop0Key))
        else:
            cursor.execute("""
                INSERT INTO hop0_keymap (hop0_key, lammps_types)
                VALUES (?, ?)
            """, (hop0Key, json.dumps([lammpsType])))

    # ── Bonded parameter operations ───────────────────────────────────────────

    def insertBondParam(self, key: tuple, params: dict) -> None:
        """Insert or merge a bond parameter entry.

        When the same key already exists (same atom type pair), parameters
        are averaged with the existing values.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB) in lexicographic order.
            params: Dict with 'k' (force constant) and 'r0' (equilibrium distance).
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("""
            SELECT coeffs_list FROM bond_params WHERE hop0_key_a = ? AND hop0_key_b = ?
        """, (key[0], key[1]))
        row = cursor.fetchone()

        new_coeff = [round(params["k"], 4), round(params["r0"], 4)]
        if row is not None:
            coeffs_list = json.loads(row["coeffs_list"])
            coeffs_list.append(new_coeff)
            cursor.execute("""
                UPDATE bond_params SET coeffs_list = ?
                WHERE hop0_key_a = ? AND hop0_key_b = ?
            """, (json.dumps(coeffs_list), key[0], key[1]))
        else:
            cursor.execute("""
                INSERT INTO bond_params (hop0_key_a, hop0_key_b, coeffs, coeffs_list)
                VALUES (?, ?, ?, ?)
            """, (key[0], key[1], json.dumps(new_coeff), json.dumps([new_coeff])))

    def lookupBondParam(self, hop0KeyA: str, hop0KeyB: str) -> dict | None:
        """Look up bond parameter by hop0 key pair.

        Keys are stored in lexicographic order, so the input order does not matter.

        Args:
            hop0KeyA: hop0 key of first atom.
            hop0KeyB: hop0 key of second atom.

        Returns:
            Bond parameter dict if found, None otherwise.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        if hop0KeyA <= hop0KeyB:
            cursor.execute("""
                SELECT coeffs FROM bond_params WHERE hop0_key_a = ? AND hop0_key_b = ?
            """, (hop0KeyA, hop0KeyB))
        else:
            cursor.execute("""
                SELECT coeffs FROM bond_params WHERE hop0_key_a = ? AND hop0_key_b = ?
            """, (hop0KeyB, hop0KeyA))
        row = cursor.fetchone()
        if row is None:
            return None
        coeffs = json.loads(row["coeffs"])
        return {"k": coeffs[0], "r0": coeffs[1]}

    def insertAngleParam(self, key: tuple, params: dict) -> None:
        """Insert or merge an angle parameter entry.

        When the same key already exists, parameters are averaged.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB, hop0KeyC) where B is the center atom.
                 Stored in lexicographic order of outer atoms (A and C can be swapped).
            params: Dict with 'k' (force constant) and 'theta0' (equilibrium angle).
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("""
            SELECT coeffs_list FROM angle_params
            WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ?
        """, (key[0], key[1], key[2]))
        row = cursor.fetchone()

        new_coeff = [round(params["k"], 3), round(params["theta0"], 3)]
        if row is not None:
            coeffs_list = json.loads(row["coeffs_list"])
            coeffs_list.append(new_coeff)
            cursor.execute("""
                UPDATE angle_params SET coeffs_list = ?
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ?
            """, (json.dumps(coeffs_list), key[0], key[1], key[2]))
        else:
            cursor.execute("""
                INSERT INTO angle_params (hop0_key_a, hop0_key_b, hop0_key_c, coeffs, coeffs_list)
                VALUES (?, ?, ?, ?, ?)
            """, (key[0], key[1], key[2],
                  json.dumps(new_coeff), json.dumps([new_coeff])))

    def lookupAngleParam(self, hop0KeyA: str, hop0KeyB: str, hop0KeyC: str) -> dict | None:
        """Look up angle parameter by hop0 key triple.

        The center atom (B) is fixed, outer atoms (A and C) can be swapped.

        Args:
            hop0KeyA: hop0 key of first atom.
            hop0KeyB: hop0 key of the center atom.
            hop0KeyC: hop0 key of third atom.

        Returns:
            Angle parameter dict if found, None otherwise.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        if hop0KeyA <= hop0KeyC:
            cursor.execute("""
                SELECT coeffs FROM angle_params
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ?
            """, (hop0KeyA, hop0KeyB, hop0KeyC))
        else:
            cursor.execute("""
                SELECT coeffs FROM angle_params
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ?
            """, (hop0KeyC, hop0KeyB, hop0KeyA))
        row = cursor.fetchone()
        if row is None:
            return None
        coeffs = json.loads(row["coeffs"])
        return {"k": coeffs[0], "theta0": coeffs[1]}

    def insertDihedralParam(self, key: tuple, params: dict) -> None:
        """Insert or merge a dihedral parameter entry.

        When the same key already exists, coefficient arrays are averaged.

        Args:
            key: Tuple of (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD) in canonical order.
                 Outer atoms (A and D) can be swapped.
            params: Dict with 'coeffs' list (e.g., [0, 0, 0.3, 0] for OPLS dihedral).
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("""
            SELECT coeffs_list FROM dihedral_params
            WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
        """, (key[0], key[1], key[2], key[3]))
        row = cursor.fetchone()

        new_coeffs = params.get("coeffs", [])
        rounded = [round(c, 3) for c in new_coeffs]
        if row is not None:
            coeffs_list = json.loads(row["coeffs_list"])
            coeffs_list.append(rounded)
            cursor.execute("""
                UPDATE dihedral_params SET coeffs_list = ?
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
            """, (json.dumps(coeffs_list), key[0], key[1], key[2], key[3]))
        else:
            cursor.execute("""
                INSERT INTO dihedral_params (hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d, coeffs, coeffs_list)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (key[0], key[1], key[2], key[3],
                  json.dumps(rounded),
                  json.dumps([rounded])))

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
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        if hop0KeyA <= hop0KeyD:
            cursor.execute("""
                SELECT coeffs FROM dihedral_params
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
            """, (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD))
        else:
            cursor.execute("""
                SELECT coeffs FROM dihedral_params
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
            """, (hop0KeyD, hop0KeyC, hop0KeyB, hop0KeyA))
        row = cursor.fetchone()
        if row is None:
            return None
        return {"coeffs": json.loads(row["coeffs"])}

    def insertImproperParam(self, key: tuple, params: dict) -> None:
        """Insert or merge an improper parameter entry.

        The first atom in the key is the center atom and is fixed.
        The remaining three atoms are sorted for canonical ordering.

        Args:
            key: Tuple of (hop0KeyCenter, hop0KeyA, hop0KeyB, hop0KeyC).
            params: Dict with 'coeffs' list (e.g., [0, -1, 2] for harmonic improper).
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("""
            SELECT coeffs_list FROM improper_params
            WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
        """, (key[0], key[1], key[2], key[3]))
        row = cursor.fetchone()

        new_coeffs = params.get("coeffs", [])
        rounded = [round(c, 3) for c in new_coeffs[:-2]] + [int(new_coeffs[-2]), int(new_coeffs[-1])]
        if row is not None:
            coeffs_list = json.loads(row["coeffs_list"])
            coeffs_list.append(rounded)
            cursor.execute("""
                UPDATE improper_params SET coeffs_list = ?
                WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
            """, (json.dumps(coeffs_list), key[0], key[1], key[2], key[3]))
        else:
            cursor.execute("""
                INSERT INTO improper_params (hop0_key_a, hop0_key_b, hop0_key_c, hop0_key_d, coeffs, coeffs_list)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (key[0], key[1], key[2], key[3],
                  json.dumps(rounded),
                  json.dumps([rounded])))

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
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        others = sorted([hop0KeyB, hop0KeyC, hop0KeyD])
        cursor.execute("""
            SELECT coeffs FROM improper_params
            WHERE hop0_key_a = ? AND hop0_key_b = ? AND hop0_key_c = ? AND hop0_key_d = ?
        """, (hop0KeyA, others[0], others[1], others[2]))
        row = cursor.fetchone()
        if row is None:
            return None
        return {"coeffs": json.loads(row["coeffs"])}

    # ── Merge and export operations ────────────────────────────────────────────

    def mergeSample(self, sampleData: dict) -> None:
        """Merge data from a sample into the database.

        This is used during add-samples to incrementally update the database
        with new samples without rebuilding from scratch.

        Args:
            sampleData: Dictionary containing atom_types, hop1_keymap, hop0_keymap,
                       and other parameter tables to merge.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        for hop2Key, info in sampleData.get("atom_types", {}).items():
            existing = self.getAtomType(hop2Key)
            if existing is not None:
                merged_sources = list(set(existing["source"]) | set(info["source"]))
                self.insertAtomType(hop2Key, {**info, "source": merged_sources})
            else:
                self.insertAtomType(hop2Key, info)

        for hop1Key, entry in sampleData.get("hop1_keymap", {}).items():
            cursor = self._conn.cursor()
            cursor.execute("SELECT lammps_types FROM hop1_keymap WHERE hop1_key = ?", (hop1Key,))
            row = cursor.fetchone()
            if row is not None:
                existing_types = set(json.loads(row["lammps_types"]))
                new_types = set(entry["lammps_types"])
                merged_types = sorted(existing_types | new_types)
                cursor.execute("""
                    UPDATE hop1_keymap SET lammps_types = ? WHERE hop1_key = ?
                """, (json.dumps(merged_types), hop1Key))
            else:
                cursor.execute("""
                    INSERT INTO hop1_keymap (hop1_key, hop0_key, lammps_types)
                    VALUES (?, ?, ?)
                """, (hop1Key, entry["hop0_key"], json.dumps(entry["lammps_types"])))

        for hop0Key, entry in sampleData.get("hop0_keymap", {}).items():
            cursor = self._conn.cursor()
            cursor.execute("SELECT lammps_types FROM hop0_keymap WHERE hop0_key = ?", (hop0Key,))
            row = cursor.fetchone()
            if row is not None:
                existing_types = set(json.loads(row["lammps_types"]))
                new_types = set(entry["lammps_types"])
                merged_types = sorted(existing_types | new_types)
                cursor.execute("""
                    UPDATE hop0_keymap SET lammps_types = ? WHERE hop0_key = ?
                """, (json.dumps(merged_types), hop0Key))
            else:
                cursor.execute("""
                    INSERT INTO hop0_keymap (hop0_key, lammps_types)
                    VALUES (?, ?)
                """, (hop0Key, json.dumps(entry["lammps_types"])))

        self._conn.commit()

    def getMeta(self, key: str) -> str | None:
        """Get a meta value by key.

        Args:
            key: Meta key to retrieve.

        Returns:
            Meta value as string, or None if not found.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("SELECT value FROM meta WHERE key = ?", (key,))
        row = cursor.fetchone()
        return row["value"] if row else None

    def setMeta(self, key: str, value: str) -> None:
        """Set a meta value.

        Args:
            key: Meta key to set.
            value: Meta value as string.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)
        """, (key, value))

    def export(self) -> dict:
        """Export the database as a plain dictionary.

        Returns a dictionary representation of all tables.

        Returns:
            Dictionary with all database tables.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        result: dict[str, Any] = {}

        # Export meta
        cursor = self._conn.cursor()
        cursor.execute("SELECT key, value FROM meta")
        result["meta"] = {}
        for row in cursor.fetchall():
            result["meta"][row["key"]] = row["value"]

        # Export atom_types
        result["atom_types"] = {}
        cursor.execute("SELECT * FROM atom_types")
        for row in cursor.fetchall():
            result["atom_types"][row["hop2_key"]] = {
                "element": row["element"],
                "hop0_key": row["hop0_key"],
                "lammps_type": row["lammps_type"],
                "mass": row["mass"],
                "sigma": row["sigma"],
                "epsilon": row["epsilon"],
                "source": json.loads(row["source"]),
            }

        # Export hop1_keymap
        result["hop1_keymap"] = {}
        cursor.execute("SELECT * FROM hop1_keymap")
        for row in cursor.fetchall():
            result["hop1_keymap"][row["hop1_key"]] = {
                "hop0_key": row["hop0_key"],
                "lammps_types": json.loads(row["lammps_types"]),
            }

        # Export hop0_keymap
        result["hop0_keymap"] = {}
        cursor.execute("SELECT * FROM hop0_keymap")
        for row in cursor.fetchall():
            result["hop0_keymap"][row["hop0_key"]] = {
                "lammps_types": json.loads(row["lammps_types"]),
            }

        # Export bond_params
        result["bond_params"] = {}
        cursor.execute("SELECT * FROM bond_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"])
            result["bond_params"][key] = {"k": row["k"], "r0": row["r0"]}

        # Export angle_params
        result["angle_params"] = {}
        cursor.execute("SELECT * FROM angle_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"])
            result["angle_params"][key] = {"k": row["k"], "theta0": row["theta0"]}

        # Export dihedral_params
        result["dihedral_params"] = {}
        cursor.execute("SELECT * FROM dihedral_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"], row["hop0_key_d"])
            result["dihedral_params"][key] = {"coeffs": json.loads(row["coeffs"])}

        # Export improper_params
        result["improper_params"] = {}
        cursor.execute("SELECT * FROM improper_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"], row["hop0_key_d"])
            result["improper_params"][key] = {"coeffs": json.loads(row["coeffs"])}

        return result

    # ── Property accessors ─────────────────────────────────────────────────────

    @property
    def atomTypes(self) -> dict:
        """Get the atom_types table (hop2 level).

        Returns:
            Dictionary mapping hop2 keys to atom type info.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[str, dict] = {}
        cursor.execute("SELECT * FROM atom_types")
        for row in cursor.fetchall():
            result[row["hop2_key"]] = {
                "element": row["element"],
                "hop0_key": row["hop0_key"],
                "lammps_type": row["lammps_type"],
                "mass": row["mass"],
                "sigma": row["sigma"],
                "epsilon": row["epsilon"],
                "charge": row["charge"],
                "source": json.loads(row["source"]),
            }
        return result

    @property
    def hop1Keymap(self) -> dict:
        """Get the hop1_keymap table.

        Returns:
            Dictionary mapping hop1 keys to {lammps_types}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[str, dict] = {}
        cursor.execute("SELECT * FROM hop1_keymap")
        for row in cursor.fetchall():
            result[row["hop1_key"]] = {
                "lammps_types": json.loads(row["lammps_types"]),
            }
        return result

    @property
    def hop0Keymap(self) -> dict:
        """Get the hop0_keymap table.

        Returns:
            Dictionary mapping hop0 keys to {lammps_types}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[str, dict] = {}
        cursor.execute("SELECT * FROM hop0_keymap")
        for row in cursor.fetchall():
            result[row["hop0_key"]] = {
                "lammps_types": json.loads(row["lammps_types"]),
            }
        return result

    @property
    def bondParams(self) -> dict:
        """Get the bond_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB) to {k, r0}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[tuple, dict] = {}
        cursor.execute("SELECT * FROM bond_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"])
            coeffs = json.loads(row["coeffs"])
            result[key] = {"k": coeffs[0], "r0": coeffs[1]}
        return result

    @property
    def angleParams(self) -> dict:
        """Get the angle_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB, hop0KeyC) to {k, theta0}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[tuple, dict] = {}
        cursor.execute("SELECT * FROM angle_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"])
            coeffs = json.loads(row["coeffs"])
            result[key] = {"k": coeffs[0], "theta0": coeffs[1]}
        return result

    @property
    def dihedralParams(self) -> dict:
        """Get the dihedral_params table.

        Returns:
            Dictionary mapping (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD) to {coeffs}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[tuple, dict] = {}
        cursor.execute("SELECT * FROM dihedral_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"], row["hop0_key_d"])
            result[key] = {"coeffs": json.loads(row["coeffs"])}
        return result

    @property
    def improperParams(self) -> dict:
        """Get the improper_params table.

        Returns:
            Dictionary mapping (hop0KeyCenter, hop0KeyA, hop0KeyB, hop0KeyC) to {coeffs}.
        """
        if self._conn is None:
            raise RuntimeError("Database not loaded")

        cursor = self._conn.cursor()
        result: dict[tuple, dict] = {}
        cursor.execute("SELECT * FROM improper_params")
        for row in cursor.fetchall():
            key = (row["hop0_key_a"], row["hop0_key_b"], row["hop0_key_c"], row["hop0_key_d"])
            result[key] = {"coeffs": json.loads(row["coeffs"])}
        return result
