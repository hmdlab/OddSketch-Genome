#!/usr/bin/env python3
import argparse, subprocess, sys, os, csv

def cmd_download(args):
    # assembly_summary.txt からコメント・ヘッダ除去し、ftp_path 列(20列目)を抽出
    with open(args.summary) as f, open(args.out, 'w') as fw:
        reader = csv.reader((l for l in f if not l.startswith('#')), delimiter='\t')
        next(reader)  # ヘッダをスキップ
        for i,row in enumerate(reader):
            if args.filter and row[4] != args.filter:
                continue
            ftp = row[19]
            fw.write(ftp + "\n")
            if args.limit and i+1 >= args.limit:
                break
    # 実ダウンロードは一括 or 後述の sketch でパイプしても OK

def cmd_sketch(args):
    # oddsketch 実行。stdin で .fna.gz ファイル名リストを受け取り、
    # 出力はバイナリ *.sketch または TSV など好みで。
    # ここでは簡易に「oddsketch sketch」を想定。
    proc = subprocess.Popen(
        ['oddsketch', 'sketch',
         '--threads', str(args.threads),
         '--out', args.out],
        stdin=open(args.list),
        stdout=sys.stdout
    )
    proc.wait()

def cmd_dist(args):
    # oddsketch dist 実行。stdin でスケッチファイル名リストを受け取り、
    # TSV (query, target, distance) を標準出力へ。
    proc = subprocess.Popen(
        ['oddsketch', 'dist',
         '--threads', str(args.threads)],
        stdin=open(args.list),
        stdout=sys.stdout
    )
    proc.wait()

def main():
    p = argparse.ArgumentParser(prog='oddpipe')
    sub = p.add_subparsers(dest='cmd')
    # download
    d = sub.add_parser('download')
    d.add_argument('--summary', required=True)
    d.add_argument('--filter', default=None,
                   help='assembly_summary の refseq_category フィルタ')
    d.add_argument('--limit', type=int, default=None,
                   help='取得する行数上限')
    d.add_argument('--out', required=True,
                   help='ftp_path 一覧を出力するファイル')
    d.set_defaults(func=cmd_download)
    # sketch
    s = sub.add_parser('sketch')
    s.add_argument('--list', required=True,
                   help='.fna/.fna.gz のファイル名リスト（1行1パス）')
    s.add_argument('--threads', type=int, default=8)
    s.add_argument('--out', required=True,
                   help='スケッチ出力プレフィックス or ファイル')
    s.set_defaults(func=cmd_sketch)
    # dist
    t = sub.add_parser('dist')
    t.add_argument('--list', required=True,
                   help='.sketch ファイル名リスト（1行1パス）')
    t.add_argument('--threads', type=int, default=8)
    t.set_defaults(func=cmd_dist)

    args = p.parse_args()
    if not args.cmd:
        p.print_help(); sys.exit(1)
    args.func(args)

if __name__=='__main__':
    main()
