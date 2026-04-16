"""Utility functions: logging configuration and default paths.

This module provides common utilities used across the mappingff package:
    - Logging setup with consistent formatting
    - Default database path constants
"""

from __future__ import annotations

import logging
from pathlib import Path


# Default database directory and file path
# Users can override these via CLI arguments
USER_DEFAULT_DB_DIR = Path("./database")
USER_DEFAULT_DB_PATH = USER_DEFAULT_DB_DIR / "db"


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
