# genome-oddsketch

OddSketch itself lives under `src/` and can be built and used independently. Benchmark and comparison workflows live under `experiments/`.

## Overview
This repository provides the OddSketch CLI and reproducible workflows for evaluating genome Jaccard similarity estimation.

OddSketch can:
- `sketch`: build sketches from genome FASTA files listed on stdin.
- `dist`: estimate Jaccard similarity between sketches using one of three explicit modes:
  - `--all-to-all` / `--alltoall`: compare all sketches listed on stdin.
  - `--bipartite --qlist ... --dblist ...`: compare every query sketch against every database sketch.
  - `--pairlist ...`: compare only sketch pairs listed in a two-column TSV.

Reproducible experiment workflows are under `experiments/`:
- `pair_task`: synthetic genome-pair benchmarks for exact Jaccard, OddSketch, and BinDash.
- `search_task`: clustered synthetic-genome search benchmarks.
- `refseq_sketch_task`: real RefSeq genome sketch-build benchmarks.

## Layout
- `src/`, `include/`: OddSketch implementation and CLI
- `experiments/pair_task`: pairwise Jaccard benchmark on synthetic genome pairs
- `experiments/search_task`: database search benchmark on clustered synthetic genomes
- `experiments/refseq_sketch_task`: real RefSeq OddSketch database build benchmark
- `experiments/tools/src/`, `experiments/tools/bin/`: experimental helper tools used by benchmarks

## Requirements
Use either the local `uv` workflow or Docker.

Local workflow:
- C++17 compiler
- Python 3.10+
- `uv sync`
- BinDash only if you want comparison runs

Docker workflow:
- Docker with Compose
- See [`README-docker.md`](README-docker.md)

## BinDash Setup
On Linux/HPC, this repository expects BinDash to live at `experiments/tools/bin/bindash`.
The following script clones BinDash, checks out `v2.6`, builds it, and installs the binary there.

```bash
bash scripts/bootstrap.sh
```

You can override the tag if needed.

```bash
BINDASH_TAG=v2.6 bash scripts/bootstrap.sh
```

## Build OddSketch
```bash
cd src
make CXX=g++ LDFLAGS=-lstdc++fs
```

This builds:
- `src/oddsketch`
- `experiments/tools/bin/true_jaccard`
- `experiments/tools/bin/true_index_pairs`

`experiments/tools/bin/true_jaccard` and `experiments/tools/bin/true_index_pairs` are experimental helper binaries used by the benchmark workflows, not part of the core OddSketch CLI surface.
Their source files live in `experiments/tools/src/`.

Core CLI examples:

```bash
src/oddsketch sketch --listfname data/oddsketch_cli_sample/lists/sample_genomes.tsv --threads=8
src/oddsketch dist --all-to-all --threads=8 < data/oddsketch_cli_sample/lists/sample_sketches.list
src/oddsketch dist --bipartite \
  --qlist data/oddsketch_cli_sample/lists/sample_queries.sketch.list \
  --dblist data/oddsketch_cli_sample/lists/sample_db.sketch.list \
  --threads=8
src/oddsketch dist \
  --pairlist data/oddsketch_cli_sample/lists/sample_sketch_pairs.tsv \
  --threads=8
```

`--listfname` expects a tab-separated file:

```text
Path-to-sequence-file<TAB>genome-name<TAB>number-of-consecutive-sequences
```

The examples above use the small FASTA samples under `data/oddsketch_cli_sample/fastas/`.
After `sketch`, generated `*.sketch` files are written next to those FASTA files.

`dist` supports three explicit modes:
- `--all-to-all` (or `--alltoall`): compare all sketches listed on stdin
- `--bipartite --qlist ... --dblist ...`: compare every query sketch against every database sketch
- `--pairlist ...`: compare only the two-column TSV pairs in the pairlist

`--pairlist` expects a tab-separated file with one sketch pair per line.

## Docker
Docker usage is documented separately in [`README-docker.md`](README-docker.md).
The repository root contains a `Dockerfile` and `docker-compose.yml` that build one image with Python dependencies, `src/oddsketch`, helper tools, and BinDash.

Main roles:
- `docker run --rm genome-oddsketch`: show `oddsketch --help`
- `docker compose run --rm oddsketch --help`: use the OddSketch CLI service
- `docker compose run --rm pair-task`: run the default pairwise benchmark
- `docker compose run --rm pair-task-sketchsize`: reproduce the sketch-size sweep
- `docker compose run --rm pair-task-bbits`: reproduce the b-bits sweep

```bash
docker compose build
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
docker compose run --rm pair-task
docker compose run --rm pair-task-sketchsize
docker compose run --rm pair-task-bbits
```

## Pairwise Benchmark
Default outputs are written under `experiments/pair_task/outputs/default/`. You can override the root output directory in `experiments/pair_task/config.json` via `paths.outdir`, or pass `--outdir` to the genome generation step.

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

Step-by-step:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
```

RMSE summary:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/<run>/results/comparison_results_oddsketch.csv \
  --csv outputs/default/<run>/results/comparison_results_bindash.csv
```

Generated files:
- FASTA pairs: `experiments/pair_task/outputs/default/<run>/genomes/`
- Pair metadata: `experiments/pair_task/outputs/default/<run>/pair_info.txt`
- Result tables: `experiments/pair_task/outputs/default/<run>/results/`
- Figures: `experiments/pair_task/outputs/default/<run>/figures/`

## Search Benchmark
Default outputs are written under `experiments/search_task/outputs/default/`. Override with `paths.outdir` in `experiments/search_task/config.json`.

```bash
bash scripts/bootstrap.sh
make -C src CXX=g++ LDFLAGS=-lstdc++fs
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

## RefSeq Sketch Benchmark
This task sketches real genomes with OddSketch and stores database size, build time, peak memory, RefSeq version/fetch metadata, and `assembly_summary_refseq.txt` under `/data`.

To download every assembly listed in `experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` first. By default this keeps only `.fna.gz` files:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

Submit it from the repository root:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

With an explicit config:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh experiments/refseq_sketch_task/config.json
```

See `experiments/README.md`, `experiments/pair_task/README.md`, `experiments/search_task/README.md`, and `experiments/refseq_sketch_task/README.md` for task-specific details.
