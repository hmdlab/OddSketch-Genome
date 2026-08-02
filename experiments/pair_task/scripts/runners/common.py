#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_repo_root() -> Path:
    return resolve_task_root().parents[1]


def resolve_path(base: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def resolve_config_path(config_arg: str) -> Path:
    task_root = resolve_task_root()
    candidates = [
        Path(config_arg).expanduser(),
        task_root / config_arg,
        Path(__file__).resolve().parent / config_arg,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise SystemExit(f"config file not found: {config_arg}")


def load_config(config_path: Path) -> dict:
    try:
        return json.loads(config_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"failed to load config {config_path}: {error}") from error


def resolve_output_root(task_root: Path, cfg: dict, cli_outdir: str | None = None) -> Path:
    if cli_outdir:
        return Path(cli_outdir).expanduser().resolve()

    paths = cfg.get("paths", {}) if isinstance(cfg, dict) else {}
    if isinstance(paths, dict) and paths.get("outdir"):
        return resolve_path(task_root, paths["outdir"])

    return (task_root / "outputs" / "smoke").resolve()


def allocate_run_dir(base_outdir: Path, prefix: str = "run") -> Path:
    for attempt in range(1000):
        stamp = datetime.now().strftime(f"{prefix}_%Y%m%d_%H%M%S_%f")
        suffix = f"_{attempt:03d}" if attempt else ""
        candidate = base_outdir / f"{stamp}_{os.getpid()}{suffix}"
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            continue
    raise SystemExit(f"could not allocate unique run directory under {base_outdir}")
