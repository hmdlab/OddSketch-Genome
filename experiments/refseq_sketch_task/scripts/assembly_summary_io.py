#!/usr/bin/env python3
"""Read and verify plain or gzip-compressed NCBI assembly summaries."""

from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path
from typing import BinaryIO, TextIO


def open_binary(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def content_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open_binary(path) as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def content_size(path: Path) -> int:
    total = 0
    with open_binary(path) as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            total += len(chunk)
    return total


def copy_uncompressed(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_binary(source) as input_file, destination.open("wb") as output_file:
        shutil.copyfileobj(input_file, output_file)
