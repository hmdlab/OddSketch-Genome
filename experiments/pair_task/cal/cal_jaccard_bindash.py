#!/usr/bin/env python3
"""
cal_diverse_bindash.py
外部クローンした BinDash を使って、生成済みFASTAペアに対する推定Jaccardを取得し、
`jaccard_bindash_results.txt` と比較用CSVを出力します。

前提:
- BinDash を `external/bindash` などに clone & build 済み
- 入力FASTAは `experiments/pair_task/make_genomes/make_diverse_genomes.py` により生成済み

設定は必ず `experiments/pair_task/pipeline_config.json` の `bindash` セクションから読み込みます。

使い方（例）:
  cd experiments/pair_task
  python cal/cal_diverse_bindash.py

備考:
- BinDash の具体的CLIはリポジトリ/バージョンに依存します。`pair_cmd`/`dist_cmd` を
  ご環境の実コマンドに合わせて設定してください。
  - 例1: "{bin} sketch --list {genome_list} --out {sketch_prefix} --threads {threads}"
  - 例2: "{bin} dist   --list {sketch_list}  --threads {threads}"
"""

import csv
import json
import os
import shlex
import subprocess
from pathlib import Path
import argparse


def read_pair_info(pair_info_path: Path):
    pairs = []
    with pair_info_path.open() as f:
        header = f.readline()
        for line in f:
            pid, f1, f2, mut, glen = line.strip().split('\t')
            pairs.append({
                'pair_id': int(pid), 'file1': f1, 'file2': f2,
                'mutation_count': int(mut), 'genome_length': int(glen)
            })
    return pairs


def run_cmd(cmd: str, capture=True):
    print(f"[run] {cmd}")
    if capture:
        p = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.stderr.strip():
            print(p.stderr.strip())
        return p.stdout
    else:
        subprocess.run(cmd, shell=True, check=True)
        return ""


def parse_estimate_output(text: str, mode: str) -> float:
    """BinDashの出力1行から推定Jaccardを抽出する。
    mode: 'tsv3_or_fraction' を想定（汎用）
    - 3列TSV: file1  file2  value
    - 5列TSV（最後が a/b 形式）: query target dist p jaccard
    """
    parts = text.strip().split('\t')
    if not parts:
        raise ValueError("empty output line")
    if mode == 'tsv3_or_fraction':
        if len(parts) == 3:
            return float(parts[2])
        if len(parts) >= 5:
            j = parts[4]
            if '/' in j:
                a, b = j.split('/')
                return float(a) / float(b)
            return float(j)
    raise ValueError(f"Unsupported parse mode: {mode}")


def ensure_outdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def _resolve_config_path(config_arg: str) -> Path:
    candidates = []
    if config_arg:
        candidates.append(Path(config_arg))
        candidates.append(Path(__file__).resolve().parent / config_arg)
    candidates.append(Path(__file__).resolve().parent / 'pipeline_config.json')
    candidates.append(Path('pipeline_config.json'))
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='pipeline_config.json', help='Path to pipeline config JSON')
    args = ap.parse_args()

    # 既定
    cfg = {
        'bindash_bin': 'bindash',
        'threads': 8,
        'outdir': 'data/test_genomes',
        'mode': 'sketch_dist',  # 'sketch_dist' (推奨) or 'pairwise'
        'pair_cmd': '{bin} exact --nthreads={threads} {f1} {f2}',
        'parse': 'tsv3_or_fraction',
        'kmerlen': 64,
        'sketchsize64': 32,
        'bbits': 16,
    }
    # config を読み込み（--config 明示指定に対応）
    cpath = _resolve_config_path(args.config)
    try:
        loaded = json.loads(cpath.read_text())
    except Exception:
        loaded = None
    if isinstance(loaded, dict):
        if 'bindash' in loaded and isinstance(loaded['bindash'], dict):
            cfg.update(loaded['bindash'])
        else:
            cfg.update(loaded)

    base = Path(__file__).resolve().parent.parent
    outdir = (base / cfg['outdir']).resolve()
    ensure_outdir(outdir)

    pair_info = base / 'data' / 'test_genomes' / 'pair_info.txt'
    if not pair_info.exists():
        raise SystemExit(f"pair_info not found: {pair_info}")
    pairs = read_pair_info(pair_info)

    # 推定を実行
    results = []
    if cfg['mode'] == 'pairwise':
        for p in pairs:
            cmd = cfg['pair_cmd'].format(bin=cfg['bindash_bin'], f1=p['file1'], f2=p['file2'], threads=cfg['threads'])
            out = run_cmd(cmd, capture=True)
            # 最終行 or 最初の非空行を採用
            line = None
            for ln in out.strip().splitlines()[::-1]:
                if ln.strip():
                    line = ln
                    break
            if not line:
                print(f"warn: no output for pair {p['pair_id']}")
                continue
            try:
                j = parse_estimate_output(line, cfg['parse'])
            except Exception as e:
                print(f"warn: parse failed for pair {p['pair_id']}: {e}")
                continue
            results.append({
                'pair_id': p['pair_id'],
                'mutation_count': p['mutation_count'],
                'genome_length': p['genome_length'],
                'jaccard_bindash': j,
                'file1': p['file1'],
                'file2': p['file2']
            })
    elif cfg['mode'] == 'sketch_dist':
        # すべてのFASTAを1つのリストにまとめて sketch -> dist
        listfile = outdir / 'bindash_sketch_list.txt'
        with listfile.open('w') as f:
            for p in pairs:
                f.write(p['file1'] + '\n')
                f.write(p['file2'] + '\n')
        sketch_path = outdir / 'bindash_sketch'
        sketch_cmd = (
            f"{cfg['bindash_bin']} sketch --listfname={listfile} "
            f"--nthreads={cfg['threads']} --kmerlen={cfg['kmerlen']} "
            f"--sketchsize64={cfg['sketchsize64']} --bbits={cfg['bbits']} "
            f"--outfname={sketch_path}"
        )
        run_cmd(sketch_cmd, capture=False)

        dist_cmd = f"{cfg['bindash_bin']} dist --nthreads={cfg['threads']} --outfname=- {sketch_path}"
        dist_out = run_cmd(dist_cmd, capture=True)

        # ペア辞書（順序無視）
        pair_index = {(p['file1'], p['file2']): p for p in pairs}
        pair_index_rev = {(p['file2'], p['file1']): p for p in pairs}

        for ln in dist_out.strip().splitlines():
            parts = ln.strip().split('\t')
            if len(parts) < 5:
                continue
            q, t, mutdist, pval, jac = parts[:5]
            if (q, t) in pair_index:
                base = pair_index[(q, t)]
            elif (t, q) in pair_index_rev:
                base = pair_index_rev[(t, q)]
            else:
                continue
            if '/' in jac:
                a, b = jac.split('/')
                j = float(a) / float(b)
            else:
                try:
                    j = float(jac)
                except ValueError:
                    continue
            results.append({
                'pair_id': base['pair_id'],
                'mutation_count': base['mutation_count'],
                'genome_length': base['genome_length'],
                'jaccard_bindash': j,
                'file1': base['file1'],
                'file2': base['file2'],
            })
    else:
        raise SystemExit("Unsupported mode. Use 'pairwise' or 'sketch_dist'.")

    # 結果の保存（テキスト）
    out_txt = outdir / 'jaccard_bindash_results.txt'
    with out_txt.open('w') as f:
        f.write('pair_id\tmutation_count\tgenome_length\tjaccard_bindash\tfile1\tfile2\n')
        for r in results:
            f.write(f"{r['pair_id']}\t{r['mutation_count']}\t{r['genome_length']}\t{r['jaccard_bindash']:.10f}\t{r['file1']}\t{r['file2']}\n")
    print(f"wrote {out_txt}")

    # Trueとマージして比較用CSV
    true_path = outdir / 'jaccard_true_results.txt'
    if true_path.exists():
        true = {}
        with true_path.open() as f:
            rd = csv.reader(f, delimiter='\t')
            next(rd)
            for row in rd:
                if not row: continue
                pid = int(row[0])
                true[pid] = {'mutation_count': int(row[1]), 'jaccard_true': float(row[4])}
        out_csv = outdir / 'comparison_results_bindash.csv'
        with out_csv.open('w') as f:
            w = csv.writer(f)
            w.writerow(['pair_id','mutation_count','jaccard_true','jaccard_bindash'])
            for r in results:
                pid = r['pair_id']
                if pid in true:
                    w.writerow([pid, true[pid]['mutation_count'], true[pid]['jaccard_true'], r['jaccard_bindash']])
        print(f"wrote {out_csv}")
    else:
        print(f"note: true results not found ({true_path}); skipped CSV merge.")


if __name__ == '__main__':
    main()
