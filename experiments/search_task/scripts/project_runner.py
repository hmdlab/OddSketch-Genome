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


def merge_tsv_files(inputs: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    with output_path.open("w") as outf:
        for in_path in inputs:
            if not in_path.exists():
                continue
            with in_path.open() as inf:
                header = next(inf, None)
                if header is None:
                    continue
                if not wrote_header:
                    outf.write(header)
                    wrote_header = True
                for line in inf:
                    if line.strip():
                        outf.write(line)
    if not wrote_header:
        raise SystemExit(f"no cluster result files to merge into {output_path}")


def manifest_has_entries(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open() as f:
        for line in f:
            if line.strip():
                return True
    return False


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


def search_scope(cfg: dict) -> str:
    return str(cfg.get("search_scope", "cluster_local"))


def run_cluster_local_workflow(scripts_dir: Path, used_config_path: Path, use_bindash: bool) -> None:
    task_root = resolve_task_root()
    cfg = json.loads(used_config_path.read_text())
    out_dir = resolve_output_root(task_root, cfg)
    manifests_root = out_dir / "data" / "manifests" / "clusters"
    cluster_dirs = sorted(path for path in manifests_root.iterdir() if path.is_dir()) if manifests_root.exists() else []
    if not cluster_dirs:
        raise SystemExit(f"cluster manifests not found under {manifests_root}")

    truth_parts = []
    truth_nn_parts = []
    odd_parts = []
    odd_nn_parts = []
    bindash_parts = []
    bindash_nn_parts = []

    print(f"\n=== Stage 2/5: Cluster-local search over {len(cluster_dirs)} clusters ===")
    for cluster_dir in cluster_dirs:
        cluster_name = cluster_dir.name
        db_list = cluster_dir / "db_genome_paths.txt"
        query_list = cluster_dir / "query_genome_paths.txt"
        if not db_list.exists() or not query_list.exists():
            raise SystemExit(f"missing cluster manifests under {cluster_dir}")
        if not manifest_has_entries(db_list) or not manifest_has_entries(query_list):
            print(f"\n--- Cluster {cluster_name}: skipped (empty db/query manifest) ---")
            continue

        truth_dir = out_dir / "results" / "truth" / "clusters" / cluster_name
        odd_dir = out_dir / "results" / "oddsketch" / "clusters" / cluster_name
        odd_intermediate_dir = out_dir / "intermediate" / "oddsketch" / "clusters" / cluster_name
        bindash_dir = out_dir / "results" / "bindash" / "clusters" / cluster_name
        bindash_intermediate_dir = out_dir / "intermediate" / "bindash" / "clusters" / cluster_name

        print(f"\n--- Cluster {cluster_name} ---")
        run(
            [
                sys.executable,
                str(scripts_dir / "true_db.py"),
                "--config",
                str(used_config_path),
                "--db-list",
                str(db_list),
                "--query-list",
                str(query_list),
                "--results-dir",
                str(truth_dir),
            ]
        )
        run(
            [
                sys.executable,
                str(scripts_dir / "oddsketch_db.py"),
                "--config",
                str(used_config_path),
                "--db-list",
                str(db_list),
                "--query-list",
                str(query_list),
                "--results-dir",
                str(odd_dir),
                "--intermediate-dir",
                str(odd_intermediate_dir),
            ]
        )
        if use_bindash:
            run(
                [
                    sys.executable,
                    str(scripts_dir / "bindash_db.py"),
                    "--config",
                    str(used_config_path),
                    "--db-list",
                    str(db_list),
                    "--query-list",
                    str(query_list),
                    "--results-dir",
                    str(bindash_dir),
                    "--intermediate-dir",
                    str(bindash_intermediate_dir),
                ]
            )

        truth_parts.append(truth_dir / "exact_query_db_jaccard.tsv")
        truth_nn_parts.append(truth_dir / "exact_top1_neighbors.tsv")
        odd_parts.append(odd_dir / "oddsketch_query_db_jaccard.tsv")
        odd_nn_parts.append(odd_dir / "oddsketch_top1_neighbors.tsv")
        if use_bindash:
            bindash_parts.append(bindash_dir / "bindash_query_db_jaccard.tsv")
            bindash_nn_parts.append(bindash_dir / "bindash_top1_neighbors.tsv")

    print("\n=== Stage 3/5: Merge cluster results ===")
    merge_tsv_files(truth_parts, out_dir / "results" / "truth" / "exact_query_db_jaccard.tsv")
    merge_tsv_files(truth_nn_parts, out_dir / "results" / "truth" / "exact_top1_neighbors.tsv")
    merge_tsv_files(odd_parts, out_dir / "results" / "oddsketch" / "oddsketch_query_db_jaccard.tsv")
    merge_tsv_files(odd_nn_parts, out_dir / "results" / "oddsketch" / "oddsketch_top1_neighbors.tsv")
    if use_bindash:
        merge_tsv_files(bindash_parts, out_dir / "results" / "bindash" / "bindash_query_db_jaccard.tsv")
        merge_tsv_files(bindash_nn_parts, out_dir / "results" / "bindash" / "bindash_top1_neighbors.tsv")


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
    scope = search_scope(cfg)

    print("[run-dir]", display_path(run_dir))
    print("[used-config]", display_path(used_config_path))
    print("[search-scope]", scope)

    print("\n=== Stage 1/5: Generate DB and query genomes ===")
    run([sys.executable, str(scripts_dir / "make_cluster_query_genomes.py"), "--config", str(used_config_path)])

    if scope == "cluster_local":
        run_cluster_local_workflow(scripts_dir, used_config_path, use_bindash)
        print("\n=== Stage 4/5: Evaluate and render figures ===")
        run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(used_config_path)])
    elif scope == "global":
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
    else:
        raise SystemExit(f"unsupported search_scope: {scope}")
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
