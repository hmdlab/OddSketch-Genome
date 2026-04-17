# Search Task

この task は、クラスタ化した合成ゲノムを生成してゲノム DB を構築し、OddSketch と BinDash の最近傍検索を比較します。

## 構成
- `config.json`: task 設定
- `scripts/`: 生成、真値計算、検索、評価、一括実行
- `analysis/`: 図示と figure 生成
- `outputs/default/`: 既定の生成物

## 基本手順
```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

一括実行:

```bash
uv run python scripts/project_runner.py --config config.json
```

繰り返し実行:

```bash
uv run python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234
```

`project_runner.py` は、設定された出力ルート配下に run ごとのディレクトリを作り、その run で使った設定を `<run>/metadata/used_config.json` に保存し、図生成まで行います。
`repeat_runner.py` は batch ディレクトリを作成し、各 run の設定を `runs/run_XXX/metadata/used_config.json` に保存します。
既定設定では、最新の設定を `outputs/default/latest_used_config.json` にも保存します。

## config.json の説明
`config.json` では、合成 DB 生成、query 生成、検索パラメータを設定します。

- `genome_length`
  - 合成ゲノム 1 本あたりの長さです。
- `clusters.num_clusters`
  - DB ゲノム生成に使うクラスタ数です。
- `clusters.cluster_size`
  - 1 クラスタあたりの DB ゲノム数です。
- `clusters.mutation_min`, `clusters.mutation_max`
  - 各クラスタ内のゲノム生成で使う変異数範囲です。
- `clusters.seed`
  - 再現性のための乱数 seed です。
- `query.num_queries`
  - 生成する query ゲノム数です。
- `query.query_mutation_min`, `query.query_mutation_max`
  - query ゲノム生成に使う変異数範囲です。
- `paths.outdir`
  - 生成物の出力ルートです。
  - 既定値: `outputs/default`
- `oddsketch.kmerlen`
  - OddSketch で使う `k` です。
- `oddsketch.sketch_size`
  - OddSketch 検索で使う sketch size です。
- `oddsketch.j0`
  - OddSketch の類似度しきい値パラメータです。
- `oddsketch.pos_mode`
  - OddSketch に渡す positional sampling mode です。
- `bindash.bindash_bin`
  - BinDash 実行ファイルの名前またはパスです。
- `bindash.enabled`
  - BinDash の検索と評価を実行するかどうかです。
  - `false` にすると OddSketch のみで実験します。
- `bindash.kmerlen`, `bindash.sketch_size`, `bindash.bbits`
  - BinDash の主要な sketch パラメータです。
  - `sketch_size` は指定したいビン数です。
  - 内部では BinDash の `--sketchsize64` に変換するため、実効値は 64 ビン単位に切り上がります。
- `bindash.threads`
  - BinDash のスレッド数です。

よくある変更:

- 出力先を変えたい:
  - `paths.outdir` を変更
- 軽いスモークテストにしたい:
  - `genome_length`, `clusters.num_clusters`, `clusters.cluster_size`, `query.num_queries` を小さくする
- 検索設定を比較したい:
  - `oddsketch.sketch_size`, `oddsketch.j0`, `bindash.sketch_size` などを変更

## 出力
既定の出力ルートは `outputs/default/` です。

- `data/db_genomes/`, `data/queries/`
- `data/manifests/db_genome_paths.txt`, `data/manifests/query_genome_paths.txt`
- `data/manifests/cluster_map.tsv`, `data/manifests/genome_mutations.tsv`
- `intermediate/oddsketch/`, `intermediate/bindash/`
- `results/truth/exact_query_db_jaccard.tsv`
- `results/truth/exact_top1_neighbors.tsv`
- `results/oddsketch/oddsketch_query_db_jaccard.tsv`
- `results/oddsketch/oddsketch_top1_neighbors.tsv`
- `results/bindash/bindash_query_db_jaccard.tsv`
- `results/bindash/bindash_top1_neighbors.tsv`
- `results/evaluation/top1_accuracy_comparison.tsv`
- `figures/`
