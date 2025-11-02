#!/usr/bin/env python3
"""
true_db.py

Compute true Jaccard for all (query, DB) pairs and label the true nearest neighbor
for each query. Implements a reusable pre-processing index:

- Preprocess each FASTA once into canonical k-mer set (k=64), 2-bit encoding → 128-bit integer per kmer
  stored as sorted pairs of uint64 in .bin, with count in .idx. This avoids re-reading FASTA in later runs.
- For each (query, db) pair, count |A∩B| via two-pointer scan on sorted arrays; J = inter / (|A|+|B|-inter).

Usage:
  cd src/project
  python cal/true_db.py --config config.json [--workers N]

Outputs (under paths.outdir):
  - true_pairs.tsv    (query, db, inter, n1, n2, jaccard_true)
  - true_nn.tsv       (query, nn_true, jaccard_true)
  - Per-genome index files next to FASTA: <fasta>.k64.bin (data), <fasta>.k64.idx (count)
"""

import argparse
import json
import os
import struct
import multiprocessing as mp
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


def encode_kmer64(seq: str) -> int:
    """Encode 64bp into 2-bit per base integer."""
    v = 0
    twobit = {'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3}
    for ch in seq:
        v = (v << 2) | twobit.get(ch, 0)
    return v


def revcomp(seq: str) -> str:
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def canonical_kmer64(seq: str) -> int:
    a = encode_kmer64(seq)
    b = encode_kmer64(revcomp(seq))
    return a if a < b else b


def build_index_for_fasta(fa_path: Path, k: int = 64) -> tuple[Path, int]:
    bin_path = Path(str(fa_path) + f'.k{k}.bin')
    idx_path = Path(str(fa_path) + f'.k{k}.idx')
    if bin_path.exists() and idx_path.exists():
        try:
            n = int(idx_path.read_text().strip())
            return bin_path, n
        except Exception:
            pass
    # build
    seqs = []
    s = []
    with fa_path.open() as f:
        for ln in f:
            if not ln:
                continue
            if ln[0] == '>':
                continue
            s.append(ln.strip())
    if not s:
        n = 0
        idx_path.write_text('0')
        bin_path.write_bytes(b'')
        return bin_path, n
    seq = ''.join(s)
    n = max(0, len(seq) - k + 1)
    kmers = set()
    for i in range(n):
        sub = seq[i:i+k]
        if len(sub) == k:
            kmers.add(canonical_kmer64(sub))
    # sort kmers and write as uint128 (hi, lo)
    arr = sorted(kmers)
    with bin_path.open('wb') as bf:
        for v in arr:
            hi = (v >> 64) & ((1<<64)-1)
            lo = v & ((1<<64)-1)
            bf.write(struct.pack('>QQ', hi, lo))
    idx_path.write_text(str(len(arr)))
    return bin_path, len(arr)


def iter_uint128_pairs(bin_path: Path):
    with bin_path.open('rb') as bf:
        data = bf.read(16)
        while data:
            hi, lo = struct.unpack('>QQ', data)
            yield (hi, lo)
            data = bf.read(16)


def intersect_count(bin_a: Path, bin_b: Path) -> tuple[int, int, int]:
    # two-pointer over sorted (hi,lo)
    # load into memory as list of tuples (acceptable for now)
    a = list(iter_uint128_pairs(bin_a))
    b = list(iter_uint128_pairs(bin_b))
    i = j = inter = 0
    na = len(a)
    nb = len(b)
    while i < na and j < nb:
        if a[i] == b[j]:
            inter += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return inter, na, nb


def compute_pairs(true_args):
    q, db, k = true_args
    qbin, qn = build_index_for_fasta(Path(q), k=k)
    results = []
    for d in db:
        if Path(d).name == Path(q).name:
            continue
        dbin, dn = build_index_for_fasta(Path(d), k=k)
        inter, na, nb = intersect_count(qbin, dbin)
        uni = na + nb - inter if (na or nb) else 1
        jac = inter / uni
        results.append((q, d, inter, na, nb, jac))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    ap.add_argument('--workers', type=int, default=max(1, os.cpu_count() or 1))
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

    db = [ln.strip() for ln in db_list.read_text().splitlines() if ln.strip()]
    qs = [ln.strip() for ln in q_list.read_text().splitlines() if ln.strip()]

    # Preprocess all genomes once (build indexes)
    k = int(cfg.get('true_jaccard', {}).get('kmerlen', cfg.get('genome_length', 64) and 64))
    for p in db + qs:
        build_index_for_fasta(Path(p), k=k)

    # Compute all pairs (parallel over queries)
    tasks = [(q, db, k) for q in qs]
    all_rows = []
    if args.workers and args.workers > 1:
        with mp.Pool(processes=args.workers) as pool:
            for res in pool.imap_unordered(compute_pairs, tasks):
                all_rows.extend(res)
    else:
        for t in tasks:
            all_rows.extend(compute_pairs(t))

    # Write true_pairs.tsv
    true_pairs = outdir / 'true_pairs.tsv'
    with true_pairs.open('w') as f:
        f.write('query\tdb\tinter\tn1\tn2\tjaccard_true\n')
        for q, d, inter, na, nb, jac in all_rows:
            f.write(f"{q}\t{d}\t{inter}\t{na}\t{nb}\t{jac:.10f}\n")
    print(f"[true] wrote {true_pairs}")

    # Pick true NN per query
    from collections import defaultdict
    best = defaultdict(lambda: (-1.0, None))
    for q, d, inter, na, nb, jac in all_rows:
        if jac > best[q][0]:
            best[q] = (jac, d)
    true_nn = outdir / 'true_nn.tsv'
    with true_nn.open('w') as f:
        f.write('query\tnn_true\tjaccard_true\n')
        for q, (jac, d) in best.items():
            f.write(f"{q}\t{d}\t{jac:.10f}\n")
    print(f"[true] wrote {true_nn}")


if __name__ == '__main__':
    main()
