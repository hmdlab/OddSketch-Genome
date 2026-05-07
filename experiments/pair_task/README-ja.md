# Pair Task

この task では、合成ゲノムペアを生成し、厳密 Jaccard と OddSketch / BinDash の推定値を比較します。

## 構成
- `config.json`: task 設定
- `scripts/`: ゲノム生成、Jaccard 計算、task runner
- `analysis/`: 図示、RMSE 集計、figure 生成
- `outputs/default/`: 既定の生成物

## 基本手順
```bash
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
```

個別実行:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
```

`project_runner.py` は、設定された出力ルート配下に run ごとのディレクトリを作り、その run で使った設定を `<run>/metadata/used_config.json` に保存し、図生成まで行います。
既定設定では、最新の解決済み設定を `outputs/default/latest_used_config.json` にも保存します。

この task の OddSketch は現在、次の一括処理で動きます。
- `sketch`: `pair_info.txt` に含まれるユニークな genome 全体を 1 回の `oddsketch sketch` で sketch 化
- `dist`: 生成された sketch pair を `oddsketch dist --pairlist=...` で一括評価

RMSE 集計:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

## config.json の説明
`config.json` では、合成データ生成と推定パラメータの両方を設定します。

- `paths.outdir`
  - 生成物の出力ルートです。
  - 既定値: `outputs/default`
- `make_genomes.genome_length`
  - 合成ゲノム 1 本あたりの長さです。
- `make_genomes.num_pairs`
  - 生成するゲノムペア数です。
- `make_genomes.mutation_min`, `make_genomes.mutation_max`
  - 各ペアに入れる変異数の最小値と最大値です。
  - 各ペアの変異数はこの範囲でサンプリングされます。
- `make_genomes.seed_base`
  - 再現性のための乱数 seed の基準値です。
- `true_jaccard.kmerlen`
  - 真値 Jaccard 計算で使う `k` です。
- `oddsketch.kmerlen`
  - OddSketch で使う `k` です。
- `oddsketch.sketch_size`
  - OddSketch の sketch size です。
- `oddsketch.j0`
  - OddSketch の類似度しきい値パラメータです。
- `oddsketch.pos_mode`
  - OddSketch に渡す positional sampling mode です。
- `oddsketch.canonical`
  - canonical k-mer を使うかどうかです。
  - pair_task の OddSketch script では、各 genome を 1 回だけ sketch し、必要な pair だけを `dist --pairlist` で評価します。
- `bindash.bindash_bin`
  - BinDash 実行ファイルの名前またはパスです。
- `bindash.enabled`
  - この task で BinDash の処理を実行するかどうかです。
  - `false` にすると OddSketch のみで実験します。
- `bindash.threads`
  - BinDash のスレッド数です。
- `bindash.mode`
  - BinDash の実行モードです。現在の既定値は `sketch_dist` です。
- `bindash.kmerlen`, `bindash.sketch_size`, `bindash.bbits`
  - BinDash の主要な sketch パラメータです。
  - `sketch_size` はスケッチ全体の目標 bit 数として解釈されます。
  - 内部では `64 * bbits` で割って BinDash の `--sketchsize64` に変換するため、実効メモリは `64 * bbits` bit 単位に切り上がります。
- `bindash.pair_cmd`
  - 1 ペア評価時に使う BinDash コマンドテンプレートです。

よくある変更:

- 出力先を変えたい:
  - `paths.outdir` を変更
- 軽いスモークテストにしたい:
  - `make_genomes.genome_length` と `make_genomes.num_pairs` を小さくする
- sketch 設定を比較したい:
  - `oddsketch.sketch_size`, `oddsketch.j0`, `bindash.sketch_size` などを変更

## 出力
既定の出力ルートは `outputs/default/` です。

- `genomes/`
- `pair_info.txt`, `genome_paths.txt`
- `results/jaccard_true_results.txt`
- `results/jaccard_oddsketch_results.txt`
- `results/jaccard_bindash_results.txt`
- `results/comparison_results_oddsketch.csv`
- `results/comparison_results_bindash.csv`
- `figures/`
