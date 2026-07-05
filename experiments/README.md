# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure
- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `refseq_sketch_task/`: real RefSeq OddSketch database build timing, memory, and size benchmark
- `tools/`: C++ helper tools and external-tool setup scripts used only by experiment workflows
- repository-root `Dockerfile` / `docker-compose.yml`: container environment for `pair_task` and the OddSketch CLI

Each task keeps its own defaults:
- config: `<task>/config.json`
- code: `<task>/scripts/`, `<task>/analysis/`
- generated data: `<task>/outputs/default/`

## Task Size
- `pair_task/` is the lighter synthetic benchmark. It can be used as a local smoke test after reducing `make_genomes.genome_length` and `make_genomes.num_pairs` in `pair_task/config.json`.
- `refseq_sketch_task/` is the heavy real-data benchmark. It downloads and sketches hundreds of thousands of RefSeq bacterial genomes, so it is intended for an HPC or server environment with substantial storage and runtime.

BinDash is not needed for ordinary OddSketch-Genome CLI use or OddSketch-only experiment runs. Install BinDash only when reproducing workflows that include the BinDash baseline:
- `pair_task/scripts/cal_jaccard_bindash.py`
- `pair_task` runs with `bindash.enabled=true`
- `refseq_sketch_task/scripts/refseq_bindash_sketch_runner.py`
- `refseq_sketch_task/jobs/qsub_refseq_bindash_sketch.sh`

## Quick Run
Install Python dependencies and build helper binaries from the repository root:

```bash
uv sync
make -C src CXX=g++ LDFLAGS=-lstdc++fs
```

Install BinDash only if you will run one of the BinDash baseline workflows listed above:

```bash
bash scripts/bootstrap.sh
```

Run the pairwise synthetic benchmark from its task directory:

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

Run the RefSeq sketch-build benchmark from the repository root:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json
```

Or submit it to Grid Engine:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

Containerized run from the repository root:

```bash
docker compose build
docker compose run --rm pair-task
```

OddSketch CLI examples against files under repository-root `docker-data/`:

```bash
printf '/data/genome_001.fna\n/data/genome_002.fna\n' > docker-data/genomes.list
docker compose run --rm oddsketch sketch --input-paths /data/genomes.list --threads=8
printf '%s\n' /data/genome_001.fna.sketch /data/genome_002.fna.sketch | docker compose run --rm -T oddsketch dist --all-to-all --threads=8
docker compose run --rm oddsketch dist --bipartite --qlist /data/queries.sketch.list --dblist /data/db.sketch.list --threads=8
docker compose run --rm oddsketch dist --pairlist /data/sketch_pairs.tsv --threads=8
```

## Notes
- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository. The default helper script builds it from `https://github.com/zhaoxiaofei/bindash.git` at tag `v2.6`.
- Task-local `config.json` files define the default output roots via `paths.outdir`.
- `pair_task/scripts/batch_project_runner.py` creates a unique run directory under `paths.outdir` for each config and saves the resolved config to `<run>/metadata/used_config.json`.
- Exact-Jaccard helper binaries are built from `experiments/tools/src/` into `experiments/tools/bin/` by `make -C src`.
