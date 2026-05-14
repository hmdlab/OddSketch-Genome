# Docker usage

このファイルは Docker を使った OddSketch CLI と pair_task 実験の実行方法をまとめます。
通常の README では概要だけを示し、Docker の詳しい使い方はここに分けています。

## 役割

このリポジトリの Docker 構成は、1 つの image と複数の compose service を使います。

| 名前 | 役割 |
| --- | --- |
| `genome-oddsketch:latest` | OddSketch 本体、Python 依存、実験補助ツール、BinDash を含む Docker image |
| `docker run --rm genome-oddsketch` | image の既定動作。`oddsketch --help` を表示 |
| `oddsketch` | OddSketch CLI を直接使う compose service |
| `pair-task` | `experiments/pair_task/config.json` で pair_task を 1 run 実行 |
| `pair-task-sketchsize` | `experiments/pair_task/configs_sketchsize/` の config 群を実行 |
| `pair-task-bbits` | `experiments/pair_task/configs_bbits/` の config 群を実行 |

## Build

最初に image を build します。

```bash
docker compose build
```

BinDash の tag を変えたい場合:

```bash
BINDASH_TAG=v2.6 docker compose build
```

## ローカル出力 mount

Docker container の中は、ホストとは別の filesystem です。
何も mount せずに container 内へ結果を書くと、`--rm` で container を消したときに結果も消えます。

このリポジトリでは、ホスト側のディレクトリを container 内の以下の場所にマウントしています。


| ローカル | container 内 | 用途 |
| --- | --- | --- |
| `./experiments/pair_task` | `/workspace/experiments/pair_task` | pair_task の config、script、出力 |
| `./docker-data` | `/data` | OddSketch CLI に渡す FASTA、sketch、list |


初回だけ、必要ならローカルディレクトリを作成します。

```bash
mkdir -p experiments/pair_task/outputs docker-data
```

`docker-data/` は OddSketch CLI に自分の FASTA や list を渡すための置き場です。
例えばホスト側の:

```text
./docker-data/genome_001.fna
```

は container 内では:

```text
/data/genome_001.fna
```

として見えます。


## OddSketch CLI

help:

```bash
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
```

`docker-data/` に FASTA を置くと、container 内では `/data/...` として見えます。

sketch:

```bash
printf '/data/genome_001.fna\tgenome_001\t1\n/data/genome_002.fna\tgenome_002\t1\n' > docker-data/genomes.tsv
docker compose run --rm oddsketch sketch --listfname /data/genomes.tsv --threads=8
```

dist all-to-all:

```bash
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | \
docker compose run --rm -T oddsketch dist --all-to-all --threads=8
```

dist bipartite:

```bash
docker compose run --rm oddsketch dist --bipartite \
  --qlist /data/queries.sketch.list \
  --dblist /data/db.sketch.list \
  --threads=8
```

dist pairlist:

```bash
docker compose run --rm oddsketch dist \
  --pairlist /data/sketch_pairs.tsv \
  --threads=8
```

## pair_task

default config を 1 run 実行:

```bash
docker compose run --rm pair-task
```

出力先:

```text
experiments/pair_task/outputs/default/
```

## pair_task 実験再現: sketchsize

`experiments/pair_task/configs_sketchsize/` の config 群を実行します。

```bash
docker compose run --rm pair-task-sketchsize
```

並列数を指定する場合:

```bash
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-sketchsize
```

出力先:

```text
experiments/pair_task/outputs/sketchsize/
```

## pair_task実験再現：b-bits

`experiments/pair_task/configs_bbits/` の config 群を実行します。

```bash
docker compose run --rm pair-task-bbits
```

並列数を指定する場合:

```bash
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-bbits
```

出力先:

```text
experiments/pair_task/outputs/bbits/
```

## docker run で実験を直接指定する場合

compose service を使わず、image に直接コマンドを渡すこともできます。

```bash
docker run --rm \
  -v "$PWD/experiments/pair_task/outputs:/workspace/experiments/pair_task/outputs" \
  -v "$PWD/docker-data:/data" \
  genome-oddsketch \
  uv run python experiments/pair_task/scripts/batch_project_runner.py \
    --config experiments/pair_task/config.json
```

ただし通常は `docker compose run --rm pair-task` の方が短く、mount も compose 側で揃うので扱いやすいです。
