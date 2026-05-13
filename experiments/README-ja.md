# 実験ディレクトリ

このディレクトリには、`src/` の OddSketch 本体を使った比較実験ワークフローをまとめています。

## 構成
- `pair_task/`: 合成ゲノムペアに対する Jaccard 比較
- `search_task/`: クラスタ化合成ゲノムに対する最近傍検索比較
- `refseq_sketch_task/`: 実 RefSeq ゲノムの OddSketch DB 構築時間・メモリ・サイズ計測
- `tools/`: 実験ワークフロー専用の C++ 補助ツールと外部ツール準備スクリプト
- リポジトリ root の `Dockerfile` / `docker-compose.yml`: `pair_task` と OddSketch CLI 用のコンテナ環境

各 task はそれぞれ次を持ちます。
- 設定: `<task>/config.json`
- コード: `<task>/scripts/`, `<task>/analysis/`
- 生成物: `<task>/outputs/default/`

## クイック実行
各 task のディレクトリで実行します。

```bash
uv sync
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

または:

```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

実 RefSeq ゲノムのスケッチ構築ベンチは root ディレクトリから `qsub` します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
コンテナ実行はリポジトリ root で行います:

```bash
docker compose build
docker compose run --rm pair-task
```

リポジトリ root の `docker-data/` 配下のファイルに対して OddSketch CLI を直接使う例:

```bash
printf '%s\n' /data/genome_001.fna /data/genome_002.fna | docker compose run --rm -T oddsketch sketch --threads=8
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | docker compose run --rm -T oddsketch dist --all-to-all --threads=8
docker compose run --rm oddsketch dist --bipartite --qlist /data/queries.sketch.list --dblist /data/db.sketch.list --threads=8
docker compose run --rm oddsketch dist --pairlist /data/sketch_pairs.tsv --threads=8
```

## メモ
- 特定の OddSketch バイナリを使う場合は `ODDSKETCH_BIN` を設定してください。
- BinDash は外部ツールであり、このリポジトリには同梱していません。
- 各 task の既定出力先は `config.json` の `paths.outdir` で定義しています。
- `pair_task/scripts/batch_project_runner.py` は config ごとに `paths.outdir` 配下へ run ディレクトリを作り、その run で使った設定を `<run>/metadata/used_config.json` に保存します。
- 厳密 Jaccard 計算用の補助バイナリは `experiments/tools/src/` から `make -C src` で `experiments/tools/bin/` に生成されます。
