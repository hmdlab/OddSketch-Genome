# Docker Usage

This file describes the Docker workflow for the OddSketch CLI and the `pair_task` benchmark. The main README keeps only the short overview; Docker-specific details live here.

## Components

This repository uses one image and several Compose services.

| Name | Purpose |
| --- | --- |
| `genome-oddsketch:latest` | Docker image containing OddSketch, Python dependencies, benchmark helper tools, and BinDash |
| `docker run --rm genome-oddsketch` | Default image command; shows `oddsketch --help` |
| `oddsketch` | Compose service for direct OddSketch CLI use |
| `pair-task` | Runs one `pair_task` job with `experiments/pair_task/config.json` |
| `pair-task-sketchsize` | Runs configs under `experiments/pair_task/configs/sketchsize/` |
| `pair-task-bbits` | Runs configs under `experiments/pair_task/configs/bbits/` |

## Build

Build the image:

```bash
docker compose build
```

To use a different BinDash tag:

```bash
BINDASH_TAG=v2.6 docker compose build
```

By default, the image builds BinDash from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

with `BINDASH_TAG=v2.6`. Override `BINDASH_REPO` or `BINDASH_TAG` when a different BinDash source or revision is required.

## Local Output Mounts

Container filesystems are separate from the host filesystem. Results written only inside a `--rm` container disappear when the container exits.

The Compose services mount these host directories:

| Host path | Container path | Purpose |
| --- | --- | --- |
| `./experiments/pair_task` | `/workspace/experiments/pair_task` | `pair_task` configs, scripts, and outputs |
| `./docker-data` | `/data` | FASTA files, sketches, and path lists for the OddSketch CLI |

Create local directories when needed:

```bash
mkdir -p experiments/pair_task/outputs docker-data
```

A host file such as:

```text
./docker-data/genome_001.fna
```

is visible in the container as:

```text
/data/genome_001.fna
```

## OddSketch CLI

Show help:

```bash
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
```

Build sketches:

```bash
printf '/data/genome_001.fna\n/data/genome_002.fna\n' > docker-data/genomes.list
docker compose run --rm oddsketch sketch --input-paths /data/genomes.list --threads=8
```

Run all-to-all distances:

```bash
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | \
docker compose run --rm -T oddsketch dist --all-to-all --threads=8
```

Run bipartite distances:

```bash
docker compose run --rm oddsketch dist --bipartite \
  --qlist /data/queries.sketch.list \
  --dblist /data/db.sketch.list \
  --threads=8
```

Run pair-list distances:

```bash
docker compose run --rm oddsketch dist \
  --pairlist /data/sketch_pairs.tsv \
  --threads=8
```

## pair_task

Run the default config once:

```bash
docker compose run --rm pair-task
```

Output is written under:

```text
experiments/pair_task/outputs/default/
```

## Reproduce Sketch-Size Runs

Run configs under `experiments/pair_task/configs/sketchsize/`:

```bash
docker compose run --rm pair-task-sketchsize
```

Set parallelism when needed:

```bash
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-sketchsize
```

Output is written under:

```text
experiments/pair_task/outputs/sketchsize/
```

## Reproduce b-bits Runs

Run configs under `experiments/pair_task/configs/bbits/`:

```bash
docker compose run --rm pair-task-bbits
```

Set parallelism when needed:

```bash
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-bbits
```

Output is written under:

```text
experiments/pair_task/outputs/bbits/
```

## Direct docker run

You can run the image directly instead of using Compose:

```bash
docker run --rm \
  -v "$PWD/experiments/pair_task/outputs:/workspace/experiments/pair_task/outputs" \
  -v "$PWD/docker-data:/data" \
  genome-oddsketch \
  uv run python experiments/pair_task/scripts/batch_project_runner.py \
    --config experiments/pair_task/config.json
```

The Compose services are usually shorter and keep the mounts consistent.
