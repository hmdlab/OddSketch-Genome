# 実験ディレクトリ

このディレクトリには、`src/` の OddSketch 本体を使った比較実験ワークフローをまとめています。

## 構成
- `pair_task/`: 合成ゲノムペアに対する Jaccard 比較
- `search_task/`: クラスタ化合成ゲノムに対する最近傍検索比較
- `tools/`: 実験ワークフロー専用の C++ 補助ツールと外部ツール準備スクリプト
- `env/Dockerfile`: 実験用コンテナ環境

各 task はそれぞれ次を持ちます。
- 設定: `<task>/config.json`
- コード: `<task>/scripts/`, `<task>/analysis/`
- 生成物: `<task>/outputs/default/`

## クイック実行
各 task のディレクトリで実行します。

```bash
uv sync
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

または:

```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

## メモ
- 特定の OddSketch バイナリを使う場合は `ODDSKETCH_BIN` を設定してください。
- BinDash は外部ツールであり、このリポジトリには同梱していません。
- 各 task の既定出力先は `config.json` の `paths.outdir` で定義しています。
- 厳密 Jaccard 計算用の補助バイナリは `experiments/tools/src/` から `make -C src` で `experiments/tools/bin/` に生成されます。
