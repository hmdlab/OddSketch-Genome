#!/usr/bin/env python3

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


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
        task_root / "config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (task_root / "config.json").resolve()


def load_config(config_path: Path) -> dict:
    try:
        return json.loads(config_path.read_text())
    except Exception:
        return {}


def resolve_output_root(task_root: Path, cfg: dict, cli_outdir: str | None = None) -> Path:
    if cli_outdir:
        return Path(cli_outdir).expanduser().resolve()

    paths = cfg.get("paths", {}) if isinstance(cfg, dict) else {}
    if isinstance(paths, dict) and paths.get("outdir"):
        return resolve_path(task_root, paths["outdir"])

    legacy = cfg.get("make_genomes", {}).get("outdir") if isinstance(cfg, dict) else None
    if legacy:
        legacy_path = resolve_path(task_root, legacy)
        return legacy_path.parent if legacy_path.name == "genomes" else legacy_path

    return (task_root / "outputs" / "default").resolve()


def display_path(raw: str | Path) -> str:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.resolve().relative_to(resolve_repo_root()))
    except Exception:
        return str(path)


def display_cmd(cmd: list[str]) -> str:
    return " ".join(display_path(part) for part in cmd)


def bindash_enabled(cfg_path: Path) -> bool:
    cfg = load_config(cfg_path)
    bindash_cfg = cfg.get("bindash", {})
    if not isinstance(bindash_cfg, dict):
        return True
    return bool(bindash_cfg.get("enabled", True))


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


def prepare_run_config(cfg_path: Path) -> tuple[Path, Path]:
    task_root = resolve_task_root()
    cfg = json.loads(cfg_path.read_text())
    base_outdir = resolve_output_root(task_root, cfg)
    base_outdir.mkdir(parents=True, exist_ok=True)

    run_dir = allocate_run_dir(base_outdir, prefix="run")
    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    cfg.setdefault("paths", {})
    cfg["paths"]["outdir"] = str(run_dir)

    used_config_path = metadata_dir / "used_config.json"
    config_text = json.dumps(cfg, indent=2) + "\n"
    used_config_path.write_text(config_text)
    (base_outdir / "latest_used_config.json").write_text(config_text)
    return run_dir, used_config_path
