"""Utility functions: logging configuration, paths, and pickle helpers.

This module provides common utilities used across the MacroMapFF package:
    - Logging setup with consistent formatting
    - Default database path constants
    - Directory creation helpers
    - Pickle serialization helpers for database operations
"""

from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Any


# Default database directory and file path
# Users can override these via CLI arguments
USER_DEFAULT_DB_DIR = Path("./database")
USER_DEFAULT_DB_PATH = USER_DEFAULT_DB_DIR / "db.pkl"


def ensureDir(path: Path) -> None:
    """Ensure a directory exists, creating it and parents if necessary.

    This is a safe wrapper around Path.mkdir that does nothing if
    the directory already exists.

    Args:
        path: Path to the directory (must be a directory, not a file).
    """
    path.mkdir(parents=True, exist_ok=True)


def setupLogging(level: int = logging.INFO) -> None:
    """Configure the root logger with a console handler.

    Sets up logging with a consistent format that includes timestamp,
    log level, logger name, and message.

    Args:
        level: Logging level (e.g., logging.INFO, logging.DEBUG).
               Defaults to logging.INFO.
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def pickleLoad(path: Path) -> Any:
    """Load an object from a pickle file.

    Args:
        path: Path to the pickle file.

    Returns:
        The deserialized object.
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def pickleSave(obj: Any, path: Path) -> None:
    """Save an object to a pickle file.

    Creates parent directories if they do not exist before saving.

    Args:
        obj: Object to serialize.
        path: Path to the output pickle file.
    """
    ensureDir(path.parent)
    with open(path, "wb") as f:
        pickle.dump(obj, f)
