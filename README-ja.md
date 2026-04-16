# genome-oddsketch

OddSketch 本体は `src/` に置き、単独でビルド・利用できるようにしています。比較実験のワークフローは `experiments/` 配下に分離しています。

## 構成
- `src/`, `include/`: OddSketch 本体実装と CLI
- `experiments/pair_task`: 合成ゲノムペアに対する Jaccard 比較
- `experiments/search_task`: クラスタ化合成ゲノムに対する検索比較
- `experiments/tools/src/`, `experiments/tools/bin/`: ベンチマーク用の実験補助ツール

## 必要環境
- C++17 コンパイラ
- Python 3.8+
- 推奨セットアップ: `uv sync`
- BinDash は比較実験を行う場合のみ必要

## OddSketch のビルド
```bash
cd src
make
```

生成されるもの:
- `src/oddsketch`
- `experiments/tools/bin/true_jaccard`
- `experiments/tools/bin/true_index_pairs`

`experiments/tools/bin/true_jaccard` と `experiments/tools/bin/true_index_pairs` は実験ワークフロー用の補助バイナリであり、OddSketch 本体 CLI そのものではありません。
ソースは `experiments/tools/src/` にあります。

## pair_task
既定では `experiments/pair_task/outputs/default/` 以下に入出力します。別の場所を使いたい場合は `experiments/pair_task/config.json` の `paths.outdir` を変更するか、生成ステップで `--outdir` を指定してください。

```bash
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

個別実行:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
uv run python analysis/make_figures.py
```

BinDash を使わない場合:

```bash
uv run python scripts/project_runner.py --config config.json --skip-bindash
uv run python analysis/make_figures.py
```

RMSE 集計:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

生成物:
- FASTA ペア: `experiments/pair_task/outputs/default/genomes/`
- ペア情報: `experiments/pair_task/outputs/default/pair_info.txt`
- 結果テーブル: `experiments/pair_task/outputs/default/results/`
- 図: `experiments/pair_task/outputs/default/figures/`

## search_task
既定の出力先は `experiments/search_task/outputs/default/` です。必要なら `experiments/search_task/config.json` の `paths.outdir` を変更してください。

```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

詳細は `experiments/README-ja.md`, `experiments/pair_task/README-ja.md`, `experiments/search_task/README_ja.md` を参照してください。
