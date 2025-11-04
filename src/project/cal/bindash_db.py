#!/usr/bin/env python3
"""
bindash_db.py

Build BinDash sketches for a genome DB and query set, then for each query
find the nearest neighbor in the DB by maximum Jaccard index (excluding self).

Usage:
  cd src/project
  python cal/bindash_db.py --config config.json

Inputs (paths.outdir):
  - db_genomes.list
  - queries.list

Outputs (paths.outdir):
  - bindash_db_sketch (prefix for DB sketch)
  - bindash_query_sketch (prefix for query sketch)
  - bindash_nn.tsv (query, nn, jaccard_bindash)
"""

import argparse
import json
import shlex
import subprocess
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


def run_cmd(cmd: str, capture=True) -> str:
    p = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE if capture else None,
                       stderr=subprocess.PIPE, text=True)
    return p.stdout if capture else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()

    cfg = load_cfg(resolve_config_path(args.config))
    outdir = Path(__file__).resolve().parent.parent / cfg.get('paths', {}).get('outdir', 'data')
    db_list = outdir / 'db_genomes.list'
    q_list = outdir / 'queries.list'
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    b = cfg.get('bindash', {})
    binpath = b.get('bindash_bin', 'bindash')
    k = b.get('kmerlen', 64)
    ss = b.get('sketchsize64', 256)
    bb = b.get('bbits', 16)
    th = b.get('threads', 1)

    db_prefix = outdir / 'bindash_db_sketch'
    q_prefix = outdir / 'bindash_query_sketch'

    t0 = perf_counter()
    # DB sketch
    cmd_db = (
        f"{shlex.quote(binpath)} sketch --listfname={db_list} "
        f"--nthreads={th} --kmerlen={k} --sketchsize64={ss} --bbits={bb} --outfname={db_prefix}"
    )
    run_cmd(cmd_db, capture=False)
    t1 = perf_counter()
    # Query sketch
    cmd_q = (
        f"{shlex.quote(binpath)} sketch --listfname={q_list} "
        f"--nthreads={th} --kmerlen={k} --sketchsize64={ss} --bbits={bb} --outfname={q_prefix}"
    )
    run_cmd(cmd_q, capture=False)
    t2 = perf_counter()
    # Dist (query vs db)
    cmd_dist = f"{shlex.quote(binpath)} dist --nthreads={th} --outfname=- {q_prefix} {db_prefix}"
    out = run_cmd(cmd_dist, capture=True)

    # Parse results: lines: query target mutdist pval jaccard
    # Keep max jaccard per query, excluding self
    best = {}
    qname_to_path = {Path(p.strip()).name: p.strip() for p in q_list.read_text().splitlines() if p.strip()}
    dbname_to_path = {Path(p.strip()).name: p.strip() for p in db_list.read_text().splitlines() if p.strip()}
    pairs_path = outdir / 'bindash_pairs.tsv'
    with pairs_path.open('w') as pf:
        pf.write('query\tdb\tjaccard_bindash\n')
    for ln in out.strip().splitlines():
        parts = ln.strip().split('\t')
        if len(parts) < 5:
            continue
        q, t, mutdist, pval, jac = parts[:5]
        # Map names back to FASTA if possible
        qn = Path(q).name
        tn = Path(t).name
        # Exclude self-matches
        if qn == tn:
            continue
        try:
            if '/' in jac:
                a, b = jac.split('/')
                j = float(a) / float(b)
            else:
                j = float(jac)
        except Exception:
            continue
        # record pair
        with pairs_path.open('a') as pf:
            pf.write(f"{qname_to_path.get(qn, qn)}\t{dbname_to_path.get(tn, tn)}\t{j:.10f}\n")
        cur = best.get(qn)
        if (cur is None) or (j > cur[0]):
            best[qn] = (j, tn)

    nn_path = outdir / 'bindash_nn.tsv'
    with nn_path.open('w') as f:
        f.write('query\tnn\tjaccard_bindash\n')
        for qn, (j, tn) in best.items():
            f.write(f"{qname_to_path.get(qn, qn)}\t{dbname_to_path.get(tn, tn)}\t{j:.10f}\n")

    t3 = perf_counter()
    (outdir / 'bindash_times.txt').write_text(
        f"sketch_db_sec\t{t1 - t0:.3f}\nsketch_queries_sec\t{t2 - t1:.3f}\nsearch_sec\t{t3 - t2:.3f}\n")
    print(f"[bindash] wrote {nn_path}")


if __name__ == '__main__':
    main()
