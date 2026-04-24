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
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    truth_dir = out_dir / "results" / "truth"
    odd_dir = out_dir / "results" / "oddsketch"
    bindash_dir = out_dir / "results" / "bindash"

    true_pairs = truth_dir / "exact_query_db_jaccard.tsv"
    odd_pairs = odd_dir / "oddsketch_query_db_jaccard.tsv"
    bindash_pairs = bindash_dir / "bindash_query_db_jaccard.tsv"

    if true_pairs.exists() and odd_pairs.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_est_vs_true.py"),
                "--true",
                str(true_pairs),
                "--pred",
                str(odd_pairs),
                "--pred-col",
                "jaccard_oddsketch",
                "--out",
                str(figures_dir / "oddsketch_true_vs_estimate.png"),
            ]
        )

    if true_pairs.exists() and bindash_pairs.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_est_vs_true.py"),
                "--true",
                str(true_pairs),
                "--pred",
                str(bindash_pairs),
                "--pred-col",
                "jaccard_bindash",
                "--out",
                str(figures_dir / "bindash_true_vs_estimate.png"),
            ]
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    run_dir, used_config_path = prepare_run_config(cfg_path)
    cfg = json.loads(used_config_path.read_text())
    bindash_cfg = cfg.get("bindash", {})
    use_bindash = bool(bindash_cfg.get("enabled", True)) if isinstance(bindash_cfg, dict) else True

    print("[run-dir]", display_path(run_dir))
    print("[used-config]", display_path(used_config_path))

    print("\n=== Stage 1/5: Generate DB and query genomes ===")
    run([sys.executable, str(scripts_dir / "make_cluster_query_genomes.py"), "--config", str(used_config_path)])

    print("\n=== Stage 2/5: Compute exact Jaccard truth ===")
    run([sys.executable, str(scripts_dir / "true_db.py"), "--config", str(used_config_path)])

    print("\n=== Stage 3/5: Run OddSketch search ===")
    run([sys.executable, str(scripts_dir / "oddsketch_db.py"), "--config", str(used_config_path)])

    if use_bindash:
        print("\n=== Stage 4/5: Run BinDash search ===")
        run([sys.executable, str(scripts_dir / "bindash_db.py"), "--config", str(used_config_path)])
        print("\n=== Stage 5/5: Evaluate and render figures ===")
        run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(used_config_path)])
    else:
        print("\n=== Stage 4/5: Skip BinDash (disabled by config) ===")
        print("\n=== Stage 5/5: Evaluate and render figures ===")
        run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(used_config_path)])
    generate_figures(used_config_path)

    print("\n=== Run Summary ===")
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    out_path = Path(outdir) if Path(outdir).is_absolute() else (task_root / outdir).resolve()
    print("[summary] oddsketch results ->", display_path(out_path / "results" / "oddsketch" / "oddsketch_top1_neighbors.tsv"))
    if use_bindash:
        print("[summary] bindash results   ->", display_path(out_path / "results" / "bindash" / "bindash_top1_neighbors.tsv"))
    else:
        print("[summary] bindash results   -> disabled by config")
    print("[summary] figures           ->", display_path(out_path / "figures"))


if __name__ == "__main__":
    main()
