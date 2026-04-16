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
uv run python analysis/make_figures.py
```

一括実行:

```bash
uv run python scripts/project_runner.py --config config.json
uv run python scripts/project_runner.py --config config.json --skip-bindash
uv run python analysis/make_figures.py
```

繰り返し実行:

```bash
uv run python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234
uv run python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234 --skip-bindash
```

## 出力
既定の出力ルートは `outputs/default/` です。

- `genomes/`, `queries/`
- `db_genomes.list`, `queries.list`
- `true_pairs.tsv`, `true_nn.tsv`
- `oddsketch_pairs.tsv`, `oddsketch_nn.tsv`
- `bindash_pairs.tsv`, `bindash_nn.tsv`
- `nn_eval.tsv`
- `figures/`
