#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path
from time import perf_counter


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(base: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    return path if path.is_absolute() else (base / path).resolve()


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


def resolve_optional_path(task_root: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    return path if path.is_absolute() else (task_root / path).resolve()


def resolve_oddsketch_bin() -> str:
    task_root = resolve_task_root()
    repo_root = task_root.parents[1]
    candidates = [
        os.environ.get("ODDSKETCH_BIN", ""),
        str(repo_root / "src" / "oddsketch"),
        "oddsketch",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    return "oddsketch"


def run_oddsketch_sketch(list_path: Path, cfg: dict) -> list[str]:
    odd = cfg.get("oddsketch", {})
    cmd = [
        resolve_oddsketch_bin(),
        "sketch",
        f"--kmer={odd.get('kmerlen', 64)}",
        f"--sketch-size={odd.get('sketch_size', 2048)}",
        f"--j0={odd.get('j0', 0.90)}",
        f"--pos-mode={odd.get('pos_mode', 'mix')}",
        f"--threads={odd.get('threads', 1)}",
    ]
    p = subprocess.run(cmd, stdin=open(list_path, "r"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
    return [line.strip() for line in p.stdout.splitlines() if line.strip()]


def write_path_list(list_path: Path, paths: list[str]) -> None:
    list_path.parent.mkdir(parents=True, exist_ok=True)
    list_path.write_text("".join(f"{path}\n" for path in paths))


def run_oddsketch_dist_bipartite(query_list_path: Path, db_list_path: Path, cfg: dict) -> list[str]:
    odd = cfg.get("oddsketch", {})
    cmd = [
        resolve_oddsketch_bin(),
        "dist",
        f"--qlist={query_list_path}",
        f"--dblist={db_list_path}",
        f"--kmer={odd.get('kmerlen', 64)}",
        f"--sketch-size={odd.get('sketch_size', 2048)}",
        f"--j0={odd.get('j0', 0.90)}",
        f"--pos-mode={odd.get('pos_mode', 'mix')}",
        f"--threads={odd.get('threads', 1)}",
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
    return [line for line in p.stdout.splitlines() if line.strip()]


def relocate_sketches(sketch_paths: list[str], target_dir: Path) -> list[str]:
    target_dir.mkdir(parents=True, exist_ok=True)
    relocated = []
    for raw_path in sketch_paths:
        src = Path(raw_path)
        dst = target_dir / src.name
        if src.resolve() != dst.resolve():
            if dst.exists():
                dst.unlink()
            shutil.move(str(src), str(dst))
        relocated.append(str(dst.resolve()))
    return relocated


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--db-list", default=None)
    ap.add_argument("--query-list", default=None)
    ap.add_argument("--results-dir", default=None)
    ap.add_argument("--intermediate-dir", default=None)
    args = ap.parse_args()

    task_root = resolve_task_root()
    cfg = json.loads(resolve_config_path(args.config).read_text())
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))
    manifests_dir = outdir / "data" / "manifests"
    db_list = resolve_optional_path(task_root, args.db_list) or (manifests_dir / "db_genome_paths.txt")
    q_list = resolve_optional_path(task_root, args.query_list) or (manifests_dir / "query_genome_paths.txt")
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    intermediate_dir = resolve_optional_path(task_root, args.intermediate_dir) or (outdir / "intermediate" / "oddsketch")
    results_dir = resolve_optional_path(task_root, args.results_dir) or (outdir / "results" / "oddsketch")
    results_dir.mkdir(parents=True, exist_ok=True)

    t0 = perf_counter()
    db_sketches = relocate_sketches(run_oddsketch_sketch(db_list, cfg), intermediate_dir / "db_sketches")
    t1 = perf_counter()
    qry_sketches = relocate_sketches(run_oddsketch_sketch(q_list, cfg), intermediate_dir / "query_sketches")
    t2 = perf_counter()

    db_sketch_list = intermediate_dir / "db_sketch_paths.txt"
    query_sketch_list = intermediate_dir / "query_sketch_paths.txt"
    write_path_list(db_sketch_list, db_sketches)
    write_path_list(query_sketch_list, qry_sketches)

    db_map = {
        Path(path.strip()).with_suffix(Path(path.strip()).suffix + ".sketch").name: path.strip()
        for path in db_list.read_text().splitlines() if path.strip()
    }
    q_map = {
        Path(path.strip()).with_suffix(Path(path.strip()).suffix + ".sketch").name: path.strip()
        for path in q_list.read_text().splitlines() if path.strip()
    }

    nn_path = results_dir / "oddsketch_top1_neighbors.tsv"
    pairs_path = results_dir / "oddsketch_query_db_jaccard.tsv"
    lines = run_oddsketch_dist_bipartite(query_sketch_list, db_sketch_list, cfg)
    with nn_path.open("w") as outf, pairs_path.open("w") as pf:
        outf.write("query\tnn\tjaccard_oddsketch\n")
        pf.write("query\tdb\tjaccard_oddsketch\n")
        best_by_query: dict[str, tuple[float, str]] = {}
        for line in lines:
            parts = line.split("\t")
            if len(parts) != 3:
                continue
            qsketch, dsketch, val = parts
            qname = Path(qsketch).name
            dname = Path(dsketch).name
            try:
                jaccard = float(val)
            except Exception:
                continue

            query_fasta = q_map.get(qname, qname)
            db_fasta = db_map.get(dname, dname)
            pf.write(f"{query_fasta}\t{db_fasta}\t{jaccard:.10f}\n")

            current = best_by_query.get(qname)
            if current is None or jaccard > current[0]:
                best_by_query[qname] = (jaccard, dname)

        for qname in sorted(best_by_query):
            jaccard, dname = best_by_query[qname]
            outf.write(f"{q_map.get(qname, qname)}\t{db_map.get(dname, dname)}\t{jaccard:.10f}\n")

    t3 = perf_counter()
    (results_dir / "oddsketch_timing.tsv").write_text(
        f"sketch_db_sec\t{t1 - t0:.3f}\nsketch_queries_sec\t{t2 - t1:.3f}\nsearch_sec\t{t3 - t2:.3f}\n"
    )
    print(f"[oddsketch] wrote {nn_path}")


if __name__ == "__main__":
    main()
