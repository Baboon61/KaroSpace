"""Small helpers for structured command-line progress messages."""

from __future__ import annotations

from typing import Any


def _indent(level: int) -> str:
    return "  " * max(0, int(level))


def log_step(message: Any, *, level: int = 0, flush: bool = True) -> None:
    """Print a work item in the export progress tree."""
    print(f"{_indent(level)}- {message}", flush=flush)


def log_detail(message: Any, *, level: int = 1, flush: bool = True) -> None:
    """Print a child detail/result under the current work item."""
    print(f"{_indent(level)}↳ {message}", flush=flush)


def log_warning(message: Any, *, level: int = 1, flush: bool = True) -> None:
    """Print a warning as a child detail so it stays in context."""
    log_detail(f"Warning: {message}", level=level, flush=flush)
