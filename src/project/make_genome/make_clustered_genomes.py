#!/usr/bin/env python3
"""
make_clustered_genomes.py

Generates n clusters of genomes with SNP mutations from cluster centers.
Outputs FASTA files under project data directory and list files for DB and queries.

Usage:
  cd src/project
  python make_genome/make_clustered_genomes.py --config config.json

Outputs (under paths.outdir):
  - genomes/cluster{cid}/g_{cid}_{idx}.fna
  - db_genomes.list (all genome FASTA paths)
  - queries.list (sampled from DB)
  - cluster_map.tsv (FASTA path, cluster_id)
"""

import argparse
import json
import os
import random
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


def rand_base(exclude: str) -> str:
    bases = ['A', 'C', 'G', 'T']
    if exclude in bases:
        bases.remove(exclude)
    return random.choice(bases)


def mutate_snp(seq: list, snps: int, rng: random.Random) -> None:
    n = len(seq)
    if n == 0 or snps <= 0:
        return
    positions = rng.sample(range(n), k=min(snps, n))
    for pos in positions:
        seq[pos] = rand_base(seq[pos])


def write_fasta(path: Path, name: str, seq: str, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    ap.add_argument('--override-queries', type=int, default=None, help='Override number of queries to sample')
    args = ap.parse_args()

    cpath = resolve_config_path(args.config)
    cfg = json.loads(cpath.read_text())

    genome_len = int(cfg.get('genome_length', 100000))
    clusters = cfg.get('clusters', {})
    n = int(clusters.get('cluster_num', 10))
    size = int(clusters.get('size', 1000))
    # 変異数は [min_snps_num, max_snps_num] の一様乱数
    min_snps = int(clusters.get('min_snps_num', 1))
    max_snps = int(clusters.get('max_snps_num', 1000))
    if max_snps < min_snps:
        max_snps = min_snps
    seed = int(clusters.get('seed', 1234))
    outdir = Path(__file__).resolve().parent.parent / cfg.get('paths', {}).get('outdir', 'data')
    qcfg = cfg.get('query', {})
    qnum = args.override_queries if args.override_queries is not None else int(qcfg.get('num_queries', 100))
    q_mut_min = int(qcfg.get('mutation_min', 1))
    q_mut_max = int(qcfg.get('mutation_max', max_snps))
    if q_mut_max < q_mut_min:
        q_mut_max = q_mut_min

    rng = random.Random(seed)

    genomes_dir = outdir / 'genomes'
    db_list = outdir / 'db_genomes.list'
    query_list = outdir / 'queries.list'
    cluster_map = outdir / 'cluster_map.tsv'

    genomes_dir.mkdir(parents=True, exist_ok=True)
    all_paths = []
    centers = []  # store center sequences per cluster id (1-based)
    with cluster_map.open('w') as cmap:
        cmap.write('path\tcluster_id\n')
        for cid in range(1, n + 1):
            # cluster center
            center = ''.join(rng.choice('ACGT') for _ in range(genome_len))
            centers.append(center)
            for idx in range(1, size + 1):
                seq_list = list(center)
                snps = rng.randint(min_snps, max_snps) if max_snps > 0 else 0
                mutate_snp(seq_list, snps, rng)
                seq = ''.join(seq_list)
                name = f"g_{cid}_{idx}"
                out_path = genomes_dir / f"cluster{cid}" / f"{name}.fna"
                write_fasta(out_path, name=name, seq=seq)
                all_paths.append(str(out_path))
                cmap.write(f"{out_path}\t{cid}\n")

    # DB list
    with db_list.open('w') as f:
        for p in all_paths:
            f.write(p + '\n')

    # Queries: mutate cluster centers (not sampled from DB)
    queries_dir = outdir / 'queries'
    queries_dir.mkdir(parents=True, exist_ok=True)
    q_paths = []
    # distribute queries roughly evenly across clusters
    per_cluster = max(1, qnum // n) if n > 0 else qnum
    made = 0
    for cid in range(1, n + 1):
        center = centers[cid - 1]
        loc = 0
        while made < qnum and loc < per_cluster:
            seq_list = list(center)
            qsnp = rng.randint(q_mut_min, q_mut_max) if q_mut_max > 0 else 0
            mutate_snp(seq_list, qsnp, rng)
            name = f"q_{cid}_{loc+1}"
            out_path = queries_dir / f"cluster{cid}" / f"{name}.fna"
            write_fasta(out_path, name=name, seq=''.join(seq_list))
            q_paths.append(str(out_path))
            made += 1
            loc += 1
        if made >= qnum:
            break
    # if still short due to truncation, fill randomly
    while made < qnum and n > 0:
        cid = rng.randint(1, n)
        center = centers[cid - 1]
        qsnp = rng.randint(q_mut_min, q_mut_max) if q_mut_max > 0 else 0
        seq_list = list(center)
        mutate_snp(seq_list, qsnp, rng)
        name = f"q_{cid}_{made+1}"
        out_path = queries_dir / f"cluster{cid}" / f"{name}.fna"
        write_fasta(out_path, name=name, seq=''.join(seq_list))
        q_paths.append(str(out_path))
        made += 1

    with query_list.open('w') as f:
        for p in q_paths:
            f.write(p + '\n')

    print(f"[make_genome] wrote {len(all_paths)} DB genomes to {genomes_dir}")
    print(f"[make_genome] wrote {len(q_paths)} query genomes to {queries_dir}")
    print(f"[make_genome] db_list={db_list}")
    print(f"[make_genome] queries={query_list} (N={len(q_paths)})")


if __name__ == '__main__':
    main()
