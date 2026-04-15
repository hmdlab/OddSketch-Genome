#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import tempfile
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
    ]
    p = subprocess.run(cmd, stdin=open(list_path, "r"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
    return [line.strip() for line in p.stdout.splitlines() if line.strip()]


def run_oddsketch_dist(list_paths: list[Path], cfg: dict) -> list[str]:
    odd = cfg.get("oddsketch", {})
    cmd = [
        resolve_oddsketch_bin(),
        "dist",
        f"--kmer={odd.get('kmerlen', 64)}",
        f"--sketch-size={odd.get('sketch_size', 2048)}",
        f"--j0={odd.get('j0', 0.90)}",
        f"--pos-mode={odd.get('pos_mode', 'mix')}",
    ]
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        for path in list_paths:
            f.write(str(path) + "\n")
        temp = f.name
    try:
        p = subprocess.run(cmd, stdin=open(temp, "r"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
        return [line for line in p.stdout.splitlines() if line.strip()]
    finally:
        try:
            os.unlink(temp)
        except Exception:
            pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    task_root = resolve_task_root()
    cfg = json.loads(resolve_config_path(args.config).read_text())
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))
    db_list = outdir / "db_genomes.list"
    q_list = outdir / "queries.list"
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    t0 = perf_counter()
    db_sketches = run_oddsketch_sketch(db_list, cfg)
    (outdir / "db_genomes.sketchlist").write_text("\n".join(db_sketches) + "\n")
    t1 = perf_counter()
    qry_sketches = run_oddsketch_sketch(q_list, cfg)
    (outdir / "queries.sketchlist").write_text("\n".join(qry_sketches) + "\n")
    t2 = perf_counter()

    db_map = {
        Path(path.strip()).with_suffix(Path(path.strip()).suffix + ".sketch").name: path.strip()
        for path in db_list.read_text().splitlines() if path.strip()
    }
    q_map = {
        Path(path.strip()).with_suffix(Path(path.strip()).suffix + ".sketch").name: path.strip()
        for path in q_list.read_text().splitlines() if path.strip()
    }

    nn_path = outdir / "oddsketch_nn.tsv"
    pairs_path = outdir / "oddsketch_pairs.tsv"
    with nn_path.open("w") as outf, pairs_path.open("w") as pf:
        outf.write("query\tnn\tjaccard_oddsketch\n")
        pf.write("query\tdb\tjaccard_oddsketch\n")
        for query_sketch in qry_sketches:
            lines = run_oddsketch_dist([Path(query_sketch)] + [Path(path) for path in db_sketches], cfg)
            best = None
            qname = Path(query_sketch).name
            for line in lines:
                parts = line.split("\t")
                if len(parts) != 3:
                    continue
                f1, f2, val = parts
                try:
                    jaccard = float(val)
                except Exception:
                    continue
                if qname not in (Path(f1).name, Path(f2).name):
                    continue
                other = Path(f2).name if Path(f1).name == qname else Path(f1).name
                if other == qname:
                    continue
                query_fasta = q_map.get(qname, qname)
                db_fasta = db_map.get(other, other)
                pf.write(f"{query_fasta}\t{db_fasta}\t{jaccard:.10f}\n")
                if best is None or jaccard > best[0]:
                    best = (jaccard, other)
            if best is not None:
                outf.write(f"{q_map.get(qname, qname)}\t{db_map.get(best[1], best[1])}\t{best[0]:.10f}\n")

    t3 = perf_counter()
    (outdir / "oddsketch_times.txt").write_text(
        f"sketch_db_sec\t{t1 - t0:.3f}\nsketch_queries_sec\t{t2 - t1:.3f}\nsearch_sec\t{t3 - t2:.3f}\n"
    )
    print(f"[oddsketch] wrote {nn_path}")


if __name__ == "__main__":
    main()
