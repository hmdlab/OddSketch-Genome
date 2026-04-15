#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd, cwd=None, capture=False):
    print("[run]", " ".join(str(x) for x in cmd))
    if capture:
        p = subprocess.run(cmd, cwd=cwd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return p.stdout
    subprocess.run(cmd, cwd=cwd, check=True)
    return ""


def render_pair(exp_root: Path):
    pair_dir = exp_root / "pair_task"
    cfg = json.loads((pair_dir / "config.json").read_text())
    out_dir = (pair_dir / cfg.get("paths", {}).get("outdir", "outputs/default")).resolve()
    results_dir = out_dir / "results"
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    analysis_dir = pair_dir / "analysis"
    pair_csv = results_dir / "comparison_results_oddsketch.csv"
    bindash_csv = results_dir / "comparison_results_bindash.csv"

    if pair_csv.exists():
        run([
            sys.executable,
            "plot_true_vs_oddsketch.py",
            "--csv", str(pair_csv),
            "--out", str(figures_dir / "oddsketch_true_vs_estimate.png"),
        ], cwd=analysis_dir)

    if bindash_csv.exists():
        run([
            sys.executable,
            "plot_true_vs_bindash.py",
            "--csv", str(bindash_csv),
            "--out", str(figures_dir / "bindash_true_vs_estimate.png"),
        ], cwd=analysis_dir)

    if pair_csv.exists() and bindash_csv.exists():
        rmse = run([
            sys.executable,
            "compute_rmse.py",
            "--csv", str(pair_csv),
            "--csv", str(bindash_csv),
        ], cwd=analysis_dir, capture=True)
        (figures_dir / "rmse_summary.txt").write_text(rmse)


def resolve_search_outdir(exp_root: Path):
    cfg_path = exp_root / "search_task" / "config.json"
    if not cfg_path.exists():
        return (exp_root / "search_task" / "outputs" / "default").resolve()
    cfg = json.loads(cfg_path.read_text())
    rel = cfg.get("paths", {}).get("outdir", "outputs/default")
    return (exp_root / "search_task" / rel).resolve()


def render_search(exp_root: Path):
    search_dir = exp_root / "search_task"
    analysis_dir = search_dir / "analysis"
    out_dir = resolve_search_outdir(exp_root)
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    true_pairs = out_dir / "true_pairs.tsv"
    odd_pairs = out_dir / "oddsketch_pairs.tsv"
    bds_pairs = out_dir / "bindash_pairs.tsv"

    if true_pairs.exists() and odd_pairs.exists():
        run([
            sys.executable,
            "plot_est_vs_true.py",
            "--true", str(true_pairs),
            "--pred", str(odd_pairs),
            "--pred-col", "jaccard_oddsketch",
            "--out", str(figures_dir / "oddsketch_true_vs_estimate.png"),
        ], cwd=analysis_dir)

    if true_pairs.exists() and bds_pairs.exists():
        run([
            sys.executable,
            "plot_est_vs_true.py",
            "--true", str(true_pairs),
            "--pred", str(bds_pairs),
            "--pred-col", "jaccard_bindash",
            "--out", str(figures_dir / "bindash_true_vs_estimate.png"),
        ], cwd=analysis_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["pair", "search", "all"], default="all")
    ap.add_argument("--exp-root", default=str(Path(__file__).resolve().parent.parent))
    args = ap.parse_args()

    exp_root = Path(args.exp_root).resolve()

    if args.task in ("pair", "all"):
        render_pair(exp_root)
    if args.task in ("search", "all"):
        render_search(exp_root)


if __name__ == "__main__":
    main()
