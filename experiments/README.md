# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure
- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `refseq_sketch_task/`: real RefSeq OddSketch database build timing, memory, and size benchmark
- `tools/`: C++ helper tools and external-tool setup scripts used only by experiment workflows

Each task keeps its own defaults:
- config: `<task>/config.json`
- code: `<task>/scripts/`, `<task>/analysis/`
- generated data: `<task>/outputs/default/`

## Task Size
- `pair_task/` is the lighter synthetic benchmark. It can be used as a local smoke test after reducing `make_genomes.genome_length` and `make_genomes.num_pairs` in `pair_task/config.json`.
- `refseq_sketch_task/` is the heavy real-data benchmark. It downloads and sketches hundreds of thousands of RefSeq bacterial genomes, so it is intended for an HPC or server environment with substantial storage and runtime.

The experiment workflows include BinDash baseline comparisons. Install BinDash before running the default `pair_task` configs or the RefSeq BinDash sketch benchmark.

## Quick Run
Install Python dependencies and build helper binaries from the repository root:

```bash
uv sync
make -C src CXX=g++ LDFLAGS=-lstdc++fs
```

Install BinDash for the baseline comparison steps:

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

## Notes
- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository. The default helper script builds it from `https://github.com/zhaoxiaofei/bindash.git` at tag `v2.6`.
- Task-local `config.json` files define the default output roots via `paths.outdir`.
- `pair_task/scripts/batch_project_runner.py` creates a unique run directory under `paths.outdir` for each config and saves the resolved config to `<run>/metadata/used_config.json`.
- Exact-Jaccard helper binaries are built from `experiments/tools/src/` into `experiments/tools/bin/` by `make -C src`.
