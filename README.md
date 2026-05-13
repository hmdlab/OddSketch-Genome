# genome-oddsketch

OddSketch itself lives under `src/` and can be built and used independently. Benchmark and comparison workflows live under `experiments/`.

## Layout
- `src/`, `include/`: OddSketch implementation and CLI
- `experiments/pair_task`: pairwise Jaccard benchmark on synthetic genome pairs
- `experiments/search_task`: database search benchmark on clustered synthetic genomes
- `experiments/refseq_sketch_task`: real RefSeq OddSketch database build benchmark
- `experiments/tools/src/`, `experiments/tools/bin/`: experimental helper tools used by benchmarks

## Requirements
- C++17 compiler
- Python 3.10+
- `uv sync` recommended
- BinDash only if you want comparison runs

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
printf '%s\n' genome_001.fna genome_002.fna | src/oddsketch sketch --threads=8
printf '%s\n' genome_001.fna.sketch genome_002.fna.sketch | src/oddsketch dist --all-to-all --threads=8
src/oddsketch dist --bipartite --qlist queries.sketch.list --dblist db.sketch.list --threads=8
src/oddsketch dist --pairlist sketch_pairs.tsv --threads=8
```

`dist` supports three explicit modes:
- `--all-to-all` (or `--alltoall`): compare all sketches listed on stdin
- `--bipartite --qlist ... --dblist ...`: compare every query sketch against every database sketch
- `--pairlist ...`: compare only the two-column TSV pairs in the pairlist

`--pairlist` expects a tab-separated file with one sketch pair per line.

## Docker
The repository root now contains a `Dockerfile` and `docker-compose.yml` that build one image with:
- Python dependencies installed via `uv sync`
- `src/oddsketch` built
- BinDash installed into `experiments/tools/bin/bindash`

Build once:

```bash
docker compose build
```

Show the OddSketch CLI help:

```bash
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
```

Run the default pairwise benchmark:

```bash
docker compose run --rm pair-task
```

Run the reproducibility sweeps currently used for `pair_task`:

```bash
docker compose run --rm pair-task-sketchsize
docker compose run --rm pair-task-bbits
```

To run those sweeps in parallel, set `PAIR_TASK_JOBS`:

```bash
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-sketchsize
PAIR_TASK_JOBS=4 docker compose run --rm pair-task-bbits
```

Local output mounts:
- Docker containers have their own filesystem. Files written only inside a container disappear when the container is removed.
- This `docker-compose.yml` mounts the local directory `./experiments/pair_task/outputs` into the container at `/workspace/experiments/pair_task/outputs`.
- Because of that mount, pair-task results created in Docker are visible on your host machine under `experiments/pair_task/outputs/`.
- The local directory `./docker-data` is mounted at `/data` for OddSketch CLI input files.

Create the local directories before running Docker if they do not exist:

```bash
mkdir -p experiments/pair_task/outputs docker-data
```

Output locations:
- default run: `experiments/pair_task/outputs/default/`
- sketch-size sweep: `experiments/pair_task/outputs/sketchsize/`
- b-bits sweep: `experiments/pair_task/outputs/bbits/`

Use the OddSketch CLI directly against files under `./docker-data`:

```bash
printf '%s\n' /data/genome_001.fna /data/genome_002.fna | docker compose run --rm -T oddsketch sketch --threads=8
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | docker compose run --rm -T oddsketch dist --all-to-all --threads=8
docker compose run --rm oddsketch dist --bipartite --qlist /data/queries.sketch.list --dblist /data/db.sketch.list --threads=8
docker compose run --rm oddsketch dist --pairlist /data/sketch_pairs.tsv --threads=8
```

Notes:
- `pair-task` writes to `./experiments/pair_task/outputs`
- the Docker image default command is `oddsketch --help`
- `pair-task` runs `experiments/pair_task/scripts/batch_project_runner.py --config experiments/pair_task/config.json`
- `pair-task-sketchsize` runs all configs under `experiments/pair_task/configs_sketchsize`
- `pair-task-bbits` runs all configs under `experiments/pair_task/configs_bbits`
- `oddsketch` mounts `./docker-data` at `/data`
- override the BinDash version at build time with `BINDASH_TAG=... docker compose build`

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
