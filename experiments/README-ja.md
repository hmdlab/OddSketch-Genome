# 実験ディレクトリ

このディレクトリには、`src/` の OddSketch 本体を使った比較実験ワークフローをまとめています。

## 構成
- `pair_task/`: 合成ゲノムペアに対する Jaccard 比較
- `search_task/`: クラスタ化合成ゲノムに対する最近傍検索比較
- `scripts/run_benchmark.sh`: 片方または両方の task を実行
- `scripts/make_figures.py`: 各 task の出力から図を再生成
- `env/Dockerfile`: 実験用コンテナ環境

各 task はそれぞれ次を持ちます。
- 設定: `<task>/config.json`
- コード: `<task>/scripts/`, `<task>/analysis/`
- 生成物: `<task>/outputs/default/`

## クイック実行
リポジトリルートで実行します。

```bash
uv sync
bash experiments/scripts/run_benchmark.sh --mode all
```

任意オプション:

```bash
bash experiments/scripts/run_benchmark.sh --mode pair
bash experiments/scripts/run_benchmark.sh --mode search
bash experiments/scripts/run_benchmark.sh --mode all --skip-bindash
```

## メモ
- 特定の OddSketch バイナリを使う場合は `ODDSKETCH_BIN` を設定してください。
- BinDash は外部ツールであり、このリポジトリには同梱していません。
- 各 task の既定出力先は `config.json` の `paths.outdir` で定義しています。
- 厳密 Jaccard 計算用の補助バイナリは `make -C src` で `tools/bin/` に生成されます。
