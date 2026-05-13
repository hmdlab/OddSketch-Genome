# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure
- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `search_task/`: nearest-neighbor search benchmark on clustered synthetic genomes
- `refseq_sketch_task/`: real RefSeq OddSketch database build timing, memory, and size benchmark
- `tools/`: C++ helper tools and external-tool setup scripts used only by experiment workflows
- repository-root `Dockerfile` / `docker-compose.yml`: container environment for `pair_task` and the OddSketch CLI

Each task keeps its own defaults:
- config: `<task>/config.json`
- code: `<task>/scripts/`, `<task>/analysis/`
- generated data: `<task>/outputs/default/`

## Quick Run
Run each task from its own directory:

```bash
uv sync
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

Or:

```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

The real RefSeq sketch-build benchmark is submitted from the repository root:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
Containerized run from the repository root:

```bash
docker compose build
docker compose run --rm pair-task
```

OddSketch CLI examples against files under repository-root `docker-data/`:

```bash
printf '%s\n' /data/genome_001.fna /data/genome_002.fna | docker compose run --rm -T oddsketch sketch --threads=8
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | docker compose run --rm -T oddsketch dist --all-to-all --threads=8
docker compose run --rm oddsketch dist --bipartite --qlist /data/queries.sketch.list --dblist /data/db.sketch.list --threads=8
docker compose run --rm oddsketch dist --pairlist /data/sketch_pairs.tsv --threads=8
```

## Notes
- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository.
- Task-local `config.json` files define the default output roots via `paths.outdir`.
- `pair_task/scripts/batch_project_runner.py` creates a unique run directory under `paths.outdir` for each config and saves the resolved config to `<run>/metadata/used_config.json`.
- Exact-Jaccard helper binaries are built from `experiments/tools/src/` into `experiments/tools/bin/` by `make -C src`.
