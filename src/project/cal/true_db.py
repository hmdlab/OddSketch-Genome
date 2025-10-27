#!/usr/bin/env python3
"""
true_db.py

Compute true Jaccard for all (query, DB) pairs and label the true nearest neighbor
for each query. Uses the C++ `src/cal/true_jaccard` binary for performance.

Usage:
  cd src/project
  python cal/true_db.py --config config.json

Outputs (under paths.outdir):
  - pair_info_db_query.txt   (generated pairs for C++)
  - jaccard_true_db_query.txt (C++ output)
  - pair_map.tsv             (pair_id to file paths)
  - true_nn.tsv              (query, nn_true, jaccard_true)
"""

import argparse
import json
import subprocess
from pathlib import Path


def resolve_config_path(config_arg: str) -> Path:
    if not config_arg:
        config_arg = 'config.json'
    cands = [
        Path(config_arg),
        Path(__file__).resolve().parent.parent / config_arg,
        Path(__file__).resolve().parent / config_arg,
    ]
    for p in cands:
        if p.exists():
            return p
    return Path(config_arg)


def load_cfg(cpath: Path) -> dict:
    try:
        return json.loads(cpath.read_text())
    except Exception:
        return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()

    cpath = resolve_config_path(args.config)
    cfg = load_cfg(cpath)
    base = Path(__file__).resolve().parent.parent  # src/project
    outdir = base / cfg.get('paths', {}).get('outdir', 'data')
    outdir.mkdir(parents=True, exist_ok=True)

    db_list = outdir / 'db_genomes.list'
    q_list = outdir / 'queries.list'
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    genome_len = int(cfg.get('genome_length', 100000))

    db = [ln.strip() for ln in db_list.read_text().splitlines() if ln.strip()]
    qs = [ln.strip() for ln in q_list.read_text().splitlines() if ln.strip()]

    # Generate pair_info and pair_map (exclude self-pairs)
    pair_info = outdir / 'pair_info_db_query.txt'
    pair_map = outdir / 'pair_map.tsv'
    pid = 0
    with pair_info.open('w') as pf, pair_map.open('w') as pm:
        pf.write('pair_id\tfile1\tfile2\tmutation_count\tgenome_length\n')
        pm.write('pair_id\tquery\tdb\n')
        for q in qs:
            qname = Path(q).name
            for d in db:
                if Path(d).name == qname:
                    continue
                pid += 1
                pf.write(f"{pid}\t{q}\t{d}\t0\t{genome_len}\n")
                pm.write(f"{pid}\t{q}\t{d}\n")

    # Call C++ true_jaccard
    cpp_bin = (base.parent / 'cal' / 'true_jaccard')  # src/cal/true_jaccard
    if not cpp_bin.exists():
        raise SystemExit(f"true_jaccard binary not found: {cpp_bin} (build in src/cal)")
    out_txt = outdir / 'jaccard_true_db_query.txt'
    cmd = [
        str(cpp_bin),
        f"--config={cpath}",
        f"--pair-info={pair_info}",
        f"--out={out_txt}",
    ]
    print('[run]', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    # Parse results and pick true NN per query
    # Map pair_id -> (q, d)
    pid_to_pair = {}
    with pair_map.open() as pm:
        _ = pm.readline()
        for ln in pm:
            pid_s, q, d = ln.strip().split('\t')
            pid_to_pair[int(pid_s)] = (q, d)

    # Read true Jaccards
    # Columns: pair_id, mutation_count, genome_length, mutation_rate, jaccard_true, ...
    from collections import defaultdict
    best = defaultdict(lambda: (-1.0, None))  # q -> (best_j, best_db)
    with out_txt.open() as f:
        _ = f.readline()
        for ln in f:
            parts = ln.strip().split('\t')
            if len(parts) < 5:
                continue
            pid = int(parts[0])
            j = float(parts[4])
            q, d = pid_to_pair.get(pid, (None, None))
            if q is None:
                continue
            if j > best[q][0]:
                best[q] = (j, d)

    true_nn = outdir / 'true_nn.tsv'
    with true_nn.open('w') as f:
        f.write('query\tnn_true\tjaccard_true\n')
        for q, (j, d) in best.items():
            if d is None:
                continue
            f.write(f"{q}\t{d}\t{j:.10f}\n")
    print(f"[true] wrote {true_nn}")


if __name__ == '__main__':
    main()

