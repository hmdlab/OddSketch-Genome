#!/usr/bin/env python3
"""
true_db.py

Compute true Jaccard for all (query, DB) pairs and label the true nearest neighbor
for each query. Implements a reusable pre-processing index:

- Preprocess each FASTA once into canonical k-mer set (k=64), 2-bit encoding → 128-bit integer per kmer
  stored as sorted pairs of uint64 in .bin, with count in .idx. This avoids re-reading FASTA in later runs.
- For each (query, db) pair, count |A∩B| via two-pointer scan on sorted arrays; J = inter / (|A|+|B|-inter).

Usage:
  cd experiments/search_task
  python cal/true_db.py --config config.json [--workers N]

Outputs (under paths.outdir):
  - true_pairs.tsv    (query, db, inter, n1, n2, jaccard_true)
  - true_nn.tsv       (query, nn_true, jaccard_true)
  - Per-genome index files next to FASTA: <fasta>.k64.bin (data), <fasta>.k64.idx (count)
"""

import argparse
import json
import os
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
    base = Path(__file__).resolve().parent.parent  # experiments/search_task
    outdir = base / cfg.get('paths', {}).get('outdir', 'data')
    outdir.mkdir(parents=True, exist_ok=True)

    db_list = outdir / 'db_genomes.list'
    q_list = outdir / 'queries.list'
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    # Build a combined list for preprocessing
    combined = outdir / 'all_genomes.list'
    with combined.open('w') as cf:
        cf.write(db_list.read_text())
        cf.write(q_list.read_text())
    # Collect paths for forced rebuild
    db_paths = [ln.strip() for ln in db_list.read_text().splitlines() if ln.strip()]
    q_paths = [ln.strip() for ln in q_list.read_text().splitlines() if ln.strip()]
    k = int(cfg.get('true_jaccard', {}).get('kmerlen', 64))
    # Force: remove existing index files so that preprocess overwrites every time
    for p in set(db_paths + q_paths):
        binp = Path(p + f'.k{k}.bin')
        idxp = Path(p + f'.k{k}.idx')
        if binp.exists():
            try: binp.unlink()
            except Exception: pass
        if idxp.exists():
            try: idxp.unlink()
            except Exception: pass
    cpp = base.parents[1] / 'src' / 'cal' / 'true_index_pairs'
    if not cpp.exists():
        raise SystemExit(f"binary not found: {cpp}. Build in src/ (or src/cal)")

    # Preprocess (sequential inside C++)
    cmd1 = [str(cpp), 'preprocess', '--list', str(combined), '--k', str(k)]
    print('[run]', ' '.join(cmd1))
    subprocess.run(cmd1, check=True)

    # Pairs (sequential inside C++)
    true_pairs = outdir / 'true_pairs.tsv'
    true_nn = outdir / 'true_nn.tsv'
    cmd2 = [str(cpp), 'pairs', '--qlist', str(q_list), '--dblist', str(db_list), '--out-pairs', str(true_pairs), '--out-nn', str(true_nn), '--k', str(k)]
    print('[run]', ' '.join(cmd2))
    subprocess.run(cmd2, check=True)
    print(f"[true] wrote {true_pairs}\n[true] wrote {true_nn}")


if __name__ == '__main__':
    main()
