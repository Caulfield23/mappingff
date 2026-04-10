"""Atom type resolution with three-level fallback.

This module provides the resolveAtomType function which tries to find
a LAMMPS atom type for a target atom by progressively falling back through
three levels of environment matching:

    1. hop2 exact match in atom_types table
    2. hop1 match in hop1_keymap
    3. hop0 match in hop0_keymap

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

    Tries each level in order until a match is found:
        1. hop2Key in atom_types -> exact match
        2. hop1Key in hop1_keymap -> first fallback
        3. hop0Key in hop0_keymap -> second fallback

    Args:
        hop2Key: SHA-256 key of hop2 environment.
        hop1Key: SHA-256 key of hop1 environment.
        hop0Key: SHA-256 key of hop0 environment.
        element: Element symbol (e.g., 'C', 'H', 'O').
        db: MacroMapDB instance with loaded database.

    Returns:
        Tuple of (lammpsType, hop0Key) if found, (None, None) if no match.
        The hop0Key is needed for bonded parameter lookup.
        For hop2 and hop1 matches, the hop0Key comes from the mapping.
        For hop0 match, the input hop0Key is returned.
    """
    log = logging.getLogger("fallback")

    # Level 1: hop2 exact match
    if hop2Key in db.atomTypes:
        entry = db.atomTypes[hop2Key]
        log.debug(f"  hop2 match: type={entry['lammps_type']}")
        return entry["lammps_type"], entry["hop0_key"]

    # Level 2: hop1 fallback
    if hop1Key in db.hop1Keymap:
        entry = db.hop1Keymap[hop1Key]
        lammps_type = min(entry["lammps_types"])
        log.debug(f"  hop1 fallback: type={lammps_type}")
        return lammps_type, entry["hop0_key"]

    # Level 3: hop0 fallback
    if hop0Key in db.hop0Keymap:
        entry = db.hop0Keymap[hop0Key]
        lammps_type = min(entry["lammps_types"])
        log.debug(f"  hop0 fallback: type={lammps_type}")
        return lammps_type, hop0Key

    # No match found - return None to indicate failure
    log.error(f"  NO MATCH for element {element}, hop0Key={hop0Key[:16]}...")
    return None, None
