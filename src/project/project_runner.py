#!/usr/bin/env python3
"""
project_runner.py

End-to-end runner for the clustered DB search experiment.
Steps:
  1) Generate clustered genomes (DB + query sampling)
  2) Build OddSketch DB and run nearest-neighbor search
  3) Build BinDash DB and run nearest-neighbor search
  4) Summarize simple metrics

Usage:
  cd src/project
  python project_runner.py --config config.json
"""

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd: list[str]) -> None:
    print("[run]", ' '.join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()

    base = Path(__file__).resolve().parent
    cfg_path = base / args.config
    if not cfg_path.exists():
        # try project root
        cfg_path = base / 'config.json'
    # 1) genomes
    run(['python', str(base / 'make_genome' / 'make_cluster_query_genomes.py'), '--config', str(cfg_path)])
    # 2) oddsketch
    run(['python', str(base / 'cal' / 'oddsketch_db.py'), '--config', str(cfg_path)])
    # 3) bindash
    run(['python', str(base / 'cal' / 'bindash_db.py'), '--config', str(cfg_path)])

    # 4) quick summary: counts and example head
    outdir = base / 'data'
    odd = outdir / 'oddsketch_nn.tsv'
    bds = outdir / 'bindash_nn.tsv'
    print("[summary] oddsketch_nn.tsv ->", odd)
    print("[summary] bindash_nn.tsv   ->", bds)
    if odd.exists():
        print(odd.read_text().splitlines()[:5])
    if bds.exists():
        print(bds.read_text().splitlines()[:5])


if __name__ == '__main__':
    main()
