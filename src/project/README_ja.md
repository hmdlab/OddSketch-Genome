# プロジェクト: DB 検索ベンチマーク（OddSketch vs BinDash）

このプロジェクトは、クラスタ中心から SNP を導入して合成ゲノムを生成し、ゲノム DB を構築したうえで、DB 内近傍探索（最近傍）を OddSketch と BinDash で比較します。

## ディレクトリ構成
- `config.json`: 実験設定（ゲノム長・クラスタ数・ツール設定など）。
- `make_genome/make_clustered_genomes.py`: クラスタ中心から DB 用ゲノムとクエリ用ゲノムを生成し、リストを出力。
- `cal/oddsketch_db.py`: OddSketch による DB/クエリのスケッチ化と検索。
- `cal/bindash_db.py`: BinDash による DB/クエリのスケッチ化と検索。
- `data/`: 生成物（FASTA、リスト、スケッチ、結果）。Git には無視されます。

## クイックスタート
1) DB とクエリの生成
- 例: 10 クラスタ × 各 1000 DB ゲノム、長さ 1e5 bp。DB はクラスタ中心を含み、各クラスタで合計 `cluster_size` 件（センター1 + 派生 `cluster_size-1`）。DBの SNP 数は U[clusters.mutation_min,clusters.mutation_max]、クエリの SNP 数は U[query.query_mutation_min,query.query_mutation_max]
- `cd src/project`
- `python make_genome/make_clustered_genomes.py --config config.json`

2) 真値と評価（任意）
- クエリはクラスタ中心から生成しているため、概念的な最近傍はクラスタ中心です。概念ラベルでの評価と、厳密Jaccardによるラベルでの評価を切り替え可能です。
- 厳密ラベルの作成（任意）: `python cal/true_db.py --config config.json` → `data/true_nn.tsv`
- 精度評価（top-1）:（既定は概念ラベル）
  - 概念ラベル（既定）: `python cal/evaluate_nn.py --config config.json`
    - 明示する場合は `--labels conceptual`
  - 厳密ラベル: `python cal/evaluate_nn.py --config config.json --labels true`
    - `data/true_nn.tsv` が無い場合は概念ラベルにフォールバックします。
  - 出力: `data/nn_eval.tsv`（端末に精度を表示）

3a) OddSketch 検索
- `python cal/oddsketch_db.py --config config.json`
- 出力: `data/oddsketch_nn.tsv`、所要時間 `data/oddsketch_times.txt`

3b) BinDash 検索
- `python cal/bindash_db.py --config config.json`
- 出力: `data/bindash_nn.tsv`、所要時間 `data/bindash_times.txt`

4) 一括実行
- `python project_runner.py --config config.json`（1 → 3a → 3b を実行。真値/評価は任意のため既定では含みません）

## 設定のポイント（config.json）
- `genome_length`: 各ゲノム長（例: 100000）
- `clusters`: `num_clusters`（クラスタ数）, `cluster_size`（各クラスタの個体数）, `mutation_min` / `mutation_max`（DB 側のSNP数の範囲）, `seed`
- `paths.outdir`: 生成物の出力先（既定 `data`）
- `oddsketch`: `kmerlen`, `sketch_size`, `j0`, `pos_mode`
- `bindash`: `bindash_bin`, `kmerlen`, `sketchsize64`, `bbits`, `threads`
- `query`: `num_queries`, `query_mutation_min`, `query_mutation_max`（クエリはクラスタ中心から独立に変異生成）

## メモ
- クエリはクラスタ中心から独立に変異生成し、最近傍選択時は自己一致を除外します。真値（厳密Jaccard）は任意です。
- BinDash を PATH に置けない場合は `bindash.bindash_bin` に絶対パスを指定。高い類似度を扱う場合は `sketchsize64 ≥ 188` を推奨。
- 規模が大きいため（例: 1 万 × 1e5 bp ≈ 1e9 塩基）、十分なディスク容量を確保し、初回は小規模で試すことを推奨します。
