# プロジェクト: DB 検索ベンチマーク（OddSketch vs BinDash）

このプロジェクトは、クラスタ中心から SNP を導入して合成ゲノムを生成し、ゲノム DB を構築したうえで、DB 内近傍探索（最近傍）を OddSketch と BinDash で比較します。

## ディレクトリ構成
- `config.json`: 実験設定（ゲノム長・クラスタ数・ツール設定など）。
- `make_genome/make_clustered_genomes.py`: クラスタ化ゲノムとリストを生成。
- `cal/oddsketch_db.py`: OddSketch による DB/クエリのスケッチ化と検索。
- `cal/bindash_db.py`: BinDash による DB/クエリのスケッチ化と検索。
- `data/`: 生成物（FASTA、リスト、スケッチ、結果）。Git には無視されます。

## クイックスタート
1) ゲノム生成（既定: 10 クラスタ × 各 1000 ゲノム、長さ 1e5 bp。各ゲノムのSNP数は一様分布 U[min_snps_num, max_snps_num]）
- `cd src/project`
- `python make_genome/make_clustered_genomes.py --config config.json`

2) OddSketch 検索
- `python cal/oddsketch_db.py --config config.json`
- 出力: `data/oddsketch_nn.tsv`、所要時間 `data/oddsketch_times.txt`

3) BinDash 検索
- `python cal/bindash_db.py --config config.json`
- 出力: `data/bindash_nn.tsv`、所要時間 `data/bindash_times.txt`

4) 一括実行
- `python project_runner.py --config config.json`

## 設定のポイント（config.json）
- `genome_length`: 各ゲノム長（例: 100000）
- `clusters`: `n`（クラスタ数）, `size`（各クラスタの個体数）, `min_snps_num` / `max_snps_num`（各ゲノムのSNP数の一様分布範囲）, `seed`
- `paths.outdir`: 生成物の出力先（既定 `data`）
- `oddsketch`: `kmerlen`, `sketch_size`, `j0`, `pos_mode`
- `bindash`: `bindash_bin`, `kmerlen`, `sketchsize64`, `bbits`, `threads`
- `query.num_queries`: クエリとして DB からサンプリングする件数

## メモ
- クエリは DB からサンプリングし、最近傍選択時は自己一致を除外します。
- BinDash を PATH に置けない場合は `bindash.bindash_bin` に絶対パスを指定。高い類似度を扱う場合は `sketchsize64 ≥ 188` を推奨。
- 規模が大きいため（例: 1 万 × 1e5 bp ≈ 1e9 塩基）、十分なディスク容量を確保し、初回は小規模で試すことを推奨します。
