"""Atom type resolution with four-level fallback.

This module provides the resolveAtomType function which tries to find
a LAMMPS atom type for a target atom by progressively falling back through
four levels of environment matching (all within atom_types table):

    1. hop3 exact match in atom_types (by hop3_key)
    2. hop2 fallback in atom_types (by hop2_key)
    3. hop1 fallback in atom_types (by hop1_key)
    4. hop0 fallback in atom_types (by hop0_key)

If no match is found at any level, returns (None, None, None, None) to indicate failure.
"""

from __future__ import annotations

from mappingff.db import MacroMapDB


def resolve_atom_type(
    hop3Key: str,
    hop2Key: str,
    hop1Key: str,
    hop0Key: str,
    db: MacroMapDB,
) -> tuple[int | None, str | None, str | None, str | None]:
    """Resolve atom type with four-level fallback.

    Tries each level in order using SQL index lookups:
        1. hop3_key = ? (primary key, O(1))
        2. hop2_key = ? (index, O(log n))
        3. hop1_key = ? (index, O(log n))
        4. hop0_key = ? (index, O(log n))

    Args:
        hop3Key: SHA-256 key of hop3 environment.
        hop2Key: SHA-256 key of hop2 environment.
        hop1Key: SHA-256 key of hop1 environment.
        hop0Key: SHA-256 key of hop0 environment.
        db: MacroMapDB instance with loaded database.

    Returns:
        Tuple of (lammpsType, hop0Key, hopLevel, matchedHop3Key) if found.
        hopLevel is "hop3", "hop2", "hop1", or "hop0".
        matchedHop3Key is the hop3_key from the matched database row.
        Returns (None, None, None, None) if no match.
    """
    if db._conn is None:
        raise RuntimeError("Database not loaded")

    cursor = db._conn.cursor()

    for col, level_name, key in zip(
        ["hop3_key", "hop2_key", "hop1_key", "hop0_key"],
        ["hop3", "hop2", "hop1", "hop0"],
        [hop3Key, hop2Key, hop1Key, hop0Key],
    ):
        cursor.execute(
            f"SELECT lammps_type, hop0_key, hop3_key FROM atom_types WHERE {col} = ? LIMIT 1",
            (key,),
        )
        row = cursor.fetchone()
        if row:
            return row["lammps_type"], row["hop0_key"], level_name, row["hop3_key"]

    return None, None, None, None
