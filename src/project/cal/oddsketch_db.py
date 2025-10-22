#!/usr/bin/env python3
"""
oddsketch_db.py

Build OddSketch sketches for a genome DB and query set, then for each query
find the nearest neighbor in the DB by maximum Jaccard estimate (excluding self).

Usage:
  cd src/project
  python cal/oddsketch_db.py --config config.json

Inputs (paths.outdir):
  - db_genomes.list
  - queries.list

Outputs (paths.outdir):
  - db_genomes.sketchlist
  - queries.sketchlist
  - oddsketch_nn.tsv (query, nn, jaccard_oddsketch)
"""

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path
from time import perf_counter


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


def run_oddsketch_sketch(list_path: Path, cfg: dict) -> list[str]:
    odd = cfg.get('oddsketch', {})
    kmer = odd.get('kmerlen', 64)
    ssize = odd.get('sketch_size', 2048)
    j0 = odd.get('j0', 0.90)
    pos_mode = odd.get('pos_mode', 'value')
    cmd = ['../oddsketch', 'sketch', f'--kmer={kmer}', f'--sketch-size={ssize}', f'--j0={j0}', f'--pos-mode={pos_mode}']
    p = subprocess.run(cmd, stdin=open(list_path, 'r'), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
    paths = [ln.strip() for ln in p.stdout.strip().splitlines() if ln.strip()]
    return paths


def run_oddsketch_dist(list_paths: list[Path], cfg: dict) -> list[str]:
    odd = cfg.get('oddsketch', {})
    kmer = odd.get('kmerlen', 64)
    ssize = odd.get('sketch_size', 2048)
    j0 = odd.get('j0', 0.90)
    pos_mode = odd.get('pos_mode', 'value')
    cmd = ['../oddsketch', 'dist', f'--kmer={kmer}', f'--sketch-size={ssize}', f'--j0={j0}', f'--pos-mode={pos_mode}']
    # create a temp list combining query and db paths
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        for pth in list_paths:
            f.write(str(pth) + '\n')
        temp = f.name
    try:
        p = subprocess.run(cmd, stdin=open(temp, 'r'), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
        return [ln for ln in p.stdout.strip().splitlines() if ln.strip()]
    finally:
        try:
            os.unlink(temp)
        except Exception:
            pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()

    cpath = resolve_config_path(args.config)
    cfg = load_cfg(cpath)
    outdir = Path(__file__).resolve().parent.parent / cfg.get('paths', {}).get('outdir', 'data')
    db_list = outdir / 'db_genomes.list'
    q_list = outdir / 'queries.list'
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    t0 = perf_counter()
    db_sketches = run_oddsketch_sketch(db_list, cfg)
    (outdir / 'db_genomes.sketchlist').write_text('\n'.join(db_sketches) + '\n')
    t1 = perf_counter()
    qry_sketches = run_oddsketch_sketch(q_list, cfg)
    (outdir / 'queries.sketchlist').write_text('\n'.join(qry_sketches) + '\n')
    t2 = perf_counter()

    # map FASTA -> sketch for quick lookup
    # oddsketch outputs paths mirroring input with .sketch appended
    db_map = {}
    with open(db_list) as f:
        for i, ln in enumerate(f):
            pth = ln.strip()
            if not pth:
                continue
            db_map[Path(pth).with_suffix(Path(pth).suffix + '.sketch').name] = pth

    nn_path = outdir / 'oddsketch_nn.tsv'
    with nn_path.open('w') as outf:
        outf.write('query\tnn\tjaccard_oddsketch\n')
        for qs in qry_sketches:
            # dist between [query] + all db
            lines = run_oddsketch_dist([Path(qs)] + [Path(p) for p in db_sketches], cfg)
            best = None
            qname = Path(qs).name
            for ln in lines:
                parts = ln.split('\t')
                if len(parts) != 3:
                    continue
                f1, f2, val = parts[0], parts[1], parts[2]
                try:
                    j = float(val)
                except Exception:
                    continue
                # pick only pairs where one is query
                if qname in (Path(f1).name, Path(f2).name):
                    other = Path(f2).name if Path(f1).name == qname else Path(f1).name
                    # exclude self matches if any
                    if other == qname:
                        continue
                    if (best is None) or (j > best[0]):
                        best = (j, other)
            if best is not None:
                nn_fasta = db_map.get(best[1], best[1])
                outf.write(f"{qname}\t{nn_fasta}\t{best[0]:.10f}\n")

    t3 = perf_counter()
    (outdir / 'oddsketch_times.txt').write_text(
        f"sketch_db_sec\t{t1 - t0:.3f}\nsketch_queries_sec\t{t2 - t1:.3f}\nsearch_sec\t{t3 - t2:.3f}\n")
    print(f"[oddsketch] wrote {nn_path}")


if __name__ == '__main__':
    main()

