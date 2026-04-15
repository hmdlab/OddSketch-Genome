# 実験再現パッケージ

このディレクトリは、**OddSketch vs BinDash** のベンチマーク再現用パッケージです。

- 本体実装はリポジトリルート（`src/`, `include/`）にあります。
- 実験用スクリプトと設定はこのディレクトリに分離しています。
- 生成物は `experiments/outputs/` に出力します。

## 重要: BinDash の扱い

このリポジトリには、BinDash のソースコードやバイナリを同梱しません。

- BinDash は著者による外部研究ツールです。
- 再現時は `experiments/scripts/build_tools.sh` でインストール/ビルドします。
- 必要に応じて `--ref <tag-or-commit>` で参照バージョンを固定できます。

## ディレクトリ構成

- `env/Dockerfile`: コンテナ実行環境
- `scripts/download_data.sh`: （任意）データ取得用の雛形
- `scripts/build_tools.sh`: 外部ツール（BinDash）の導入/ビルド
- `scripts/run_benchmark.sh`: ベンチマーク実行
- `scripts/make_figures.py`: 比較図の生成
- `configs/`: ベンチマーク設定
- `outputs/`: 生成物（gitignore対象）

## クイック実行（ホスト環境）

リポジトリルートで実行:

```bash
uv sync
bash experiments/scripts/build_tools.sh --method auto
bash experiments/scripts/run_benchmark.sh --mode all
```

任意オプション:

```bash
bash experiments/scripts/run_benchmark.sh --mode pair
bash experiments/scripts/run_benchmark.sh --mode search
bash experiments/scripts/run_benchmark.sh --mode all --skip-bindash
```

## Docker での再現（推奨）

リポジトリルートで実行:

```bash
docker build -f experiments/env/Dockerfile -t oddsketch-exp .
docker run --rm -it oddsketch-exp
```

必要に応じてホスト側の出力先をマウント:

```bash
docker run --rm -it \
  -v "$(pwd)/experiments/outputs:/workspace/experiments/outputs" \
  oddsketch-exp
```

## メモ

- ペア比較ベンチマークは `experiments/pair_task/data/test_genomes/` を作業データとして使用します。
- 検索ベンチマークの出力先は `configs/search_task.config.json` で設定し、既定は `experiments/outputs/search_task` です。
- 特定の OddSketch バイナリを使う場合は `ODDSKETCH_BIN` を設定して実行してください。

## 引用

BinDash 比較を使う場合は、BinDash プロジェクトの原著論文を引用してください。
