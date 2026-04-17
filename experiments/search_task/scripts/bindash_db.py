#!/usr/bin/env python3

import argparse
import json
import math
import shlex
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


def run_cmd(cmd: str, capture: bool = True) -> str:
    process = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
        text=True,
    )
    return process.stdout if capture else ""


def resolve_bindash_sketch_params(cfg: dict) -> tuple[int, int]:
    bbits = int(cfg.get("bbits", 16))
    if "sketch_size" in cfg:
        sketch_size = int(cfg["sketch_size"])
    else:
        sketch_size = 64 * int(cfg.get("sketchsize64", 256))
    if sketch_size <= 0:
        raise SystemExit("bindash.sketch_size must be positive")
    sketchsize64 = math.ceil(sketch_size / 64)
    effective_sketch_size = 64 * sketchsize64
    return effective_sketch_size, bbits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    task_root = resolve_task_root()
    cfg = json.loads(resolve_config_path(args.config).read_text())
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))
    manifests_dir = outdir / "data" / "manifests"
    db_list = manifests_dir / "db_genome_paths.txt"
    q_list = manifests_dir / "query_genome_paths.txt"
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    b = cfg.get("bindash", {})
    binpath = b.get("bindash_bin", "bindash")
    k = int(b.get("kmerlen", 64))
    sketch_size, bb = resolve_bindash_sketch_params(b)
    ss = math.ceil(sketch_size / 64)
    th = int(b.get("threads", 1))

    intermediate_dir = outdir / "intermediate" / "bindash"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    results_dir = outdir / "results" / "bindash"
    results_dir.mkdir(parents=True, exist_ok=True)

    db_prefix = intermediate_dir / "db_sketch"
    q_prefix = intermediate_dir / "query_sketch"

    t0 = perf_counter()
    run_cmd(
        f"{shlex.quote(binpath)} sketch --listfname={db_list} --nthreads={th} "
        f"--kmerlen={k} --sketchsize64={ss} --bbits={bb} --outfname={db_prefix}",
        capture=False,
    )
    t1 = perf_counter()
    run_cmd(
        f"{shlex.quote(binpath)} sketch --listfname={q_list} --nthreads={th} "
        f"--kmerlen={k} --sketchsize64={ss} --bbits={bb} --outfname={q_prefix}",
        capture=False,
    )
    t2 = perf_counter()
    out = run_cmd(f"{shlex.quote(binpath)} dist --nthreads={th} --outfname=- {q_prefix} {db_prefix}", capture=True)

    qname_to_path = {Path(path.strip()).name: path.strip() for path in q_list.read_text().splitlines() if path.strip()}
    dbname_to_path = {Path(path.strip()).name: path.strip() for path in db_list.read_text().splitlines() if path.strip()}

    best = {}
    pairs_path = results_dir / "bindash_query_db_jaccard.tsv"
    with pairs_path.open("w") as pf:
        pf.write("query\tdb\tjaccard_bindash\n")
        for line in out.splitlines():
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            q, t, _, _, jac = parts[:5]
            qn = Path(q).name
            tn = Path(t).name
            if qn == tn:
                continue
            try:
                jaccard = float(jac.split("/")[0]) / float(jac.split("/")[1]) if "/" in jac else float(jac)
            except Exception:
                continue
            pf.write(f"{qname_to_path.get(qn, qn)}\t{dbname_to_path.get(tn, tn)}\t{jaccard:.10f}\n")
            current = best.get(qn)
            if current is None or jaccard > current[0]:
                best[qn] = (jaccard, tn)

    nn_path = results_dir / "bindash_top1_neighbors.tsv"
    with nn_path.open("w") as f:
        f.write("query\tnn\tjaccard_bindash\n")
        for qn, (jaccard, tn) in best.items():
            f.write(f"{qname_to_path.get(qn, qn)}\t{dbname_to_path.get(tn, tn)}\t{jaccard:.10f}\n")

    t3 = perf_counter()
    (results_dir / "bindash_timing.tsv").write_text(
        f"sketch_db_sec\t{t1 - t0:.3f}\nsketch_queries_sec\t{t2 - t1:.3f}\nsearch_sec\t{t3 - t2:.3f}\n"
    )
    print(f"[bindash] wrote {nn_path}")


if __name__ == "__main__":
    main()
