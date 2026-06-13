"""Logging helpers for CLI and future GUI frontends."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

DEFAULT_LOG_FORMAT = "[%(levelname)s] %(message)s"


def get_logger(name: str = "sphingolipid_toolkit") -> logging.Logger:
    """Return a package logger without changing global logging state."""

    return logging.getLogger(name)


def configure_logging(
    level: int | str = logging.INFO,
    stream: TextIO | None = None,
    logger_name: str = "sphingolipid_toolkit",
) -> logging.Logger:
    """Configure a simple stream handler for command-line runs."""

    logger = get_logger(logger_name)
    logger.setLevel(level)
    if not any(getattr(handler, "_sphingolipid_default", False) for handler in logger.handlers):
        handler = logging.StreamHandler(stream)
        handler.setFormatter(logging.Formatter(DEFAULT_LOG_FORMAT))
        handler._sphingolipid_default = True  # type: ignore[attr-defined]
        logger.addHandler(handler)
    logger.propagate = False
    return logger


def attach_file_handler(
    log_file: str | Path,
    logger: logging.Logger | None = None,
    level: int | str = logging.INFO,
) -> logging.FileHandler:
    """Attach a package log file handler, replacing any previous one for that file."""

    target = logger or get_logger()
    path = Path(log_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    resolved = path.resolve()
    for handler in list(target.handlers):
        if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename).resolve() == resolved:
            target.removeHandler(handler)
            handler.close()
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(DEFAULT_LOG_FORMAT))
    handler._sphingolipid_file = True  # type: ignore[attr-defined]
    target.addHandler(handler)
    target.setLevel(level)
    target.propagate = False
    return handler


@dataclass
class MemoryLogHandler(logging.Handler):
    """A small handler that stores formatted log messages for GUI display."""

    records: list[logging.LogRecord] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        super().__init__()
        self.setFormatter(logging.Formatter(DEFAULT_LOG_FORMAT))

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)
        self.messages.append(self.format(record))


def attach_memory_handler(logger: logging.Logger | None = None) -> MemoryLogHandler:
    """Attach and return a memory handler for a logger."""

    target = logger or get_logger()
    handler = MemoryLogHandler()
    target.addHandler(handler)
    return handler
