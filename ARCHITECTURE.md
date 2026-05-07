# アーキテクチャ概要

## コンポーネント
- `oddsketch`（C++/std::thread）: ゲノム列からスケッチを生成し、スケッチ間距離を計算するコア。
- `oddpipe.py`（Python CLI）: パイプラインの司令塔。入力リストの生成と、`oddsketch` のサブコマンド呼び出しを行う。
- 解析補助（Python）: `extract_jaccard_values.py`、`src/test/*` にある比較/可視化スクリプト。

## データフロー
1) 参照取得リスト作成（download）
- 入力: `assembly_summary.txt`
- 出力: FTP/ファイルパス一覧（例: `refgenomes.list`）

2) スケッチ生成（sketch）
- 入力: `.fna/.fna.gz` のパスを列挙したリストファイル（標準入力）
- 実行: `oddsketch sketch --threads N`
- 出力: 各入力パスに対応する `<input>.sketch`

3) 距離計算（dist）
- 入力:
  - all-to-all: `.sketch` のパスを列挙したリスト（標準入力）
  - bipartite: `--qlist` と `--dblist`
  - pairlist: 2 列 TSV の `--pairlist`
- 実行:
  - `oddsketch dist --threads N`
  - `oddsketch dist --threads N --qlist queries.txt --dblist db.txt`
  - `oddsketch dist --threads N --pairlist pairs.tsv`
- 出力: 標準出力に TSV（`sketch1<TAB>sketch2<TAB>jaccard_estimate`）

```
assembly_summary.txt ──> refgenomes.list ──> .sketch 群 ──> dist.tsv
          (oddpipe download)      (oddpipe sketch)       (oddpipe dist)
```

## ディレクトリ規約
- コア: `src`（`oddsketch.cpp`, `Makefile`, `oddpipe.py`）
- ゴールデン/検証データ: `src/test/data/test_genomes/`
- 図表/スクリプト: `src/test/`
- 大容量データはコミットしない。パスは `.list`/`.sketch` で参照。

## 並列化と性能
- `oddsketch` は `--threads` でスレッド数を制御する。
- `sketch` は入力ファイル単位で並列化する。
- `dist` は all-to-all では外側 index 単位、bipartite では query 単位、pairlist では left sketch 単位で並列化する。
- 再現性のため、比較時はスレッド数・CPU/メモリ・I/O 環境を PR 説明に明記。

## 入出力の要点
- リスト形式: 1 行 1 パス（絶対/相対どちらも可）。all-to-all / sketch では標準入力で渡す。
- pairlist 形式: 1 行 2 パス、タブ区切り。
- 出力 TSV: タブ区切り。必要に応じて先頭数行をサンプルとして共有。
- 例:
  - `cat refgenomes.list | ./src/oddsketch sketch --threads 8`
  - `cat refgenomes.sketch.list | ./src/oddsketch dist --threads 8 > refgenomes.dist.tsv`
  - `./src/oddsketch dist --threads 8 --qlist queries.sketch.list --dblist db.sketch.list > query_db.dist.tsv`
  - `./src/oddsketch dist --threads 8 --pairlist pairs.tsv > pairwise.dist.tsv`

## 拡張ポイント
- 前処理/取得: `download` 段階でフィルタやミラー取得スクリプトを差し替え可能。
- 特徴量/距離: `oddsketch` に新しいスケッチ方式や距離指標を追加。
- 解析: 既存の比較スクリプトにメトリクス（RMSE, R² など）を継ぎ足し可能。
