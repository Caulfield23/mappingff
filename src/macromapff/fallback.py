"""Atom type resolution with three-level fallback.

This module provides the resolveAtomType function which tries to find
a LAMMPS atom type for a target atom by progressively falling back through
three levels of environment matching (all within atom_types table):

    1. hop2 exact match in atom_types (by hop2_key)
    2. hop1 match in atom_types (by hop1_key)
    3. hop0 match in atom_types (by hop0_key)

If no match is found at any level, returns (None, None) to indicate failure.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from macromapff.db import MacroMapDB


def resolveAtomType(
    hop2Key: str,
    hop1Key: str,
    hop0Key: str,
    element: str,
    db: MacroMapDB,
) -> tuple[int | None, str | None]:
    """Resolve atom type with three-level fallback.

    Tries each level in order by searching atom_types table:
        1. hop2Key in atom_types -> exact match
        2. hop1Key in atom_types (by hop1_key column)
        3. hop0Key in atom_types (by hop0_key column)

    Args:
        hop2Key: SHA-256 key of hop2 environment.
        hop1Key: SHA-256 key of hop1 environment.
        hop0Key: SHA-256 key of hop0 environment.
        element: Element symbol (e.g., 'C', 'H', 'O').
        db: MacroMapDB instance with loaded database.

    Returns:
        Tuple of (lammpsType, hop0Key) if found, (None, None) if no match.
    """
    log = logging.getLogger("fallback")

    # Level 1: hop2 exact match
    if hop2Key in db.atomTypes:
        entry = db.atomTypes[hop2Key]
        log.debug(f"  hop2 match: type={entry['lammps_type']}")
        return entry["lammps_type"], entry["hop0_key"]

    # Level 2: hop1 fallback - search by hop1_key column in atom_types
    for key, entry in db.atomTypes.items():
        if entry.get("hop1_key") == hop1Key:
            log.debug(f"  hop1 fallback: type={entry['lammps_type']}")
            return entry["lammps_type"], entry["hop0_key"]

    # Level 3: hop0 fallback - search by hop0_key column in atom_types
    for key, entry in db.atomTypes.items():
        if entry.get("hop0_key") == hop0Key:
            log.debug(f"  hop0 fallback: type={entry['lammps_type']}")
            return entry["lammps_type"], entry["hop0_key"]

    # No match found - return None to indicate failure
    log.error(f"  NO MATCH for element {element}, hop0Key={hop0Key[:16]}...")
    return None, None
