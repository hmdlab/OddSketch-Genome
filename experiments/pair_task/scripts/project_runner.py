#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_repo_root() -> Path:
    return resolve_task_root().parents[1]


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


def resolve_config_path(config_arg: str) -> Path:
    task_root = resolve_task_root()
    candidates = [
        Path(config_arg),
        task_root / config_arg,
        Path(__file__).resolve().parent / config_arg,
        task_root / "config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (task_root / "config.json").resolve()


def run(cmd: list[str]) -> None:
    print("[run]", display_cmd(cmd))
    subprocess.run(cmd, check=True)


def bindash_enabled(cfg_path: Path) -> bool:
    cfg = json.loads(cfg_path.read_text())
    bindash_cfg = cfg.get("bindash", {})
    if not isinstance(bindash_cfg, dict):
        return True
    return bool(bindash_cfg.get("enabled", True))


def resolve_output_root(task_root: Path, cfg: dict) -> Path:
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    path = Path(outdir)
    return path if path.is_absolute() else (task_root / path).resolve()


def allocate_run_dir(base_outdir: Path, prefix: str = "run") -> Path:
    stamp = datetime.now().strftime(f"{prefix}_%Y%m%d_%H%M%S")
    candidate = base_outdir / stamp
    suffix = 1
    while candidate.exists():
        candidate = base_outdir / f"{stamp}_{suffix:02d}"
        suffix += 1
    return candidate


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


def generate_figures(used_config_path: Path) -> None:
    task_root = resolve_task_root()
    analysis_dir = task_root / "analysis"
    out_dir = resolve_output_root(task_root, json.loads(used_config_path.read_text()))
    results_dir = out_dir / "results"
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    odd_csv = results_dir / "comparison_results_oddsketch.csv"
    bindash_csv = results_dir / "comparison_results_bindash.csv"

    if odd_csv.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_true_vs_oddsketch.py"),
                "--csv",
                str(odd_csv),
                "--out",
                str(figures_dir / "oddsketch_true_vs_estimate.png"),
            ]
        )

    if bindash_csv.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_true_vs_bindash.py"),
                "--csv",
                str(bindash_csv),
                "--out",
                str(figures_dir / "bindash_true_vs_estimate.png"),
            ]
        )

    if odd_csv.exists() and bindash_csv.exists():
        rmse = subprocess.run(
            [
                sys.executable,
                str(analysis_dir / "compute_rmse.py"),
                "--csv",
                str(odd_csv),
                "--csv",
                str(bindash_csv),
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        (figures_dir / "rmse_summary.txt").write_text(rmse.stdout)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    run_dir, used_config_path = prepare_run_config(cfg_path)
    use_bindash = bindash_enabled(used_config_path)

    print("[run-dir]", display_path(run_dir))
    print("[used-config]", display_path(used_config_path))

    run([sys.executable, str(scripts_dir / "make_genomes.py"), "--config", str(used_config_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_true.py"), "--config", str(used_config_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_oddsketch.py"), "--config", str(used_config_path)])
    if use_bindash:
        run([sys.executable, str(scripts_dir / "cal_jaccard_bindash.py"), "--config", str(used_config_path)])
    generate_figures(used_config_path)


if __name__ == "__main__":
    main()
