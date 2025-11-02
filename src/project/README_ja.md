# プロジェクト: DB 検索ベンチマーク（OddSketch vs BinDash）

このプロジェクトは、クラスタ中心から SNP を導入して合成ゲノムを生成し、ゲノム DB を構築したうえで、DB 内近傍探索（最近傍）を OddSketch と BinDash で比較します。

## ディレクトリ構成
- `config.json`: 実験設定（ゲノム長・クラスタ数・ツール設定など）。
- `make_genome/make_cluster_query_genomes.py`: クラスタ中心から DB 用ゲノムとクエリ用ゲノムを生成し、リストを出力。
- `cal/oddsketch_db.py`: OddSketch による DB/クエリのスケッチ化と検索。
- `cal/bindash_db.py`: BinDash による DB/クエリのスケッチ化と検索。
- `data/`: 生成物（FASTA、リスト、スケッチ、結果）。Git には無視されます。
  - true_pairs.tsv / oddsketch_pairs.tsv / bindash_pairs.tsv: クエリ×DBのペアごとのJaccard（真値/推定）
  - nn_eval.tsv: ツールごとのtop-1精度

## クイックスタート
1) DB とクエリの生成
- 例: 10 クラスタ × 各 1000 DB ゲノム、長さ 1e5 bp。DB はクラスタ中心を含まず、各クラスタで派生 `cluster_size` 件。DB の SNP 数は U[clusters.mutation_min,clusters.mutation_max]、クエリの SNP 数は U[query.query_mutation_min,query.query_mutation_max]
- `cd src/project`
- `python make_genome/make_cluster_query_genomes.py --config config.json`

2) 真値と評価（必須）
- 厳密Jaccard（クエリ×DB全組）: `python cal/true_db.py --config config.json`
  - 出力: `data/true_pairs.tsv`, `data/true_nn.tsv`
- 精度評価（top-1）: `python cal/evaluate_nn.py --config config.json`
  - 出力: `data/nn_eval.tsv`（端末に精度を表示）
- 真値と推定の図示:
  - OddSketch: `python analysis/plot_est_vs_true.py --true data/true_pairs.tsv --pred data/oddsketch_pairs.tsv --pred-col jaccard_oddsketch --out data/oddsketch_true_vs_estimate.png`
  - BinDash:   `python analysis/plot_est_vs_true.py --true data/true_pairs.tsv --pred data/bindash_pairs.tsv   --pred-col jaccard_bindash   --out data/bindash_true_vs_estimate.png`

3a) OddSketch 検索
- `python cal/oddsketch_db.py --config config.json`
- 出力: `data/oddsketch_nn.tsv`、所要時間 `data/oddsketch_times.txt`

3b) BinDash 検索
- `python cal/bindash_db.py --config config.json`
- 出力: `data/bindash_nn.tsv`、所要時間 `data/bindash_times.txt`

4) 一括実行
- `python project_runner.py --config config.json`（1 → 2 → 3a → 3b → 評価 まで実行）

## 設定のポイント（config.json）
- `genome_length`: 各ゲノム長（例: 100000）
- `clusters`: `num_clusters`（クラスタ数）, `cluster_size`（各クラスタの個体数）, `mutation_min` / `mutation_max`（DB 側のSNP数の範囲）, `seed`
- `paths.outdir`: 生成物の出力先（既定 `data`）
- `oddsketch`: `kmerlen`, `sketch_size`, `j0`, `pos_mode`
- `bindash`: `bindash_bin`, `kmerlen`, `sketchsize64`, `bbits`, `threads`
- `query`: `num_queries`, `query_mutation_min`, `query_mutation_max`（クエリはクラスタ中心から独立に変異生成）

## メモ
- クエリはクラスタ中心から独立に変異生成し、最近傍選択時は自己一致を除外します。評価には必ず厳密Jaccardに基づく真値を用います。
- BinDash を PATH に置けない場合は `bindash.bindash_bin` に絶対パスを指定。高い類似度を扱う場合は `sketchsize64 ≥ 188` を推奨。
- 規模が大きいため（例: 1 万 × 1e5 bp ≈ 1e9 塩基）、十分なディスク容量を確保し、初回は小規模で試すことを推奨します。
