# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure

- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `refseq_sketch_task/`: real RefSeq OddSketch database build timing, memory, and size benchmark
- `tools/`: C++ helper tools and external-tool setup scripts used only by experiment workflows

Each task keeps its own code, configuration, and ignored output directories.
`pair_task/` separates its paired and supplemental settings under
`pair_task/configs/`; `refseq_sketch_task/` uses a task-level `config.json`.

## Benchmark Workflows

- `pair_task/` is the synthetic benchmark. Its paper-scale settings are in
  `configs/paired.json` and `configs/supplementary.json`.
- `refseq_sketch_task/` is the real-data benchmark. It downloads and sketches
  hundreds of thousands of RefSeq bacterial genomes, so it requires substantial
  storage and runtime and is intended for a server or HPC environment.

Both workflows include BinDash baseline comparisons.

## Prerequisites

Run the setup commands from the repository root:

```bash
uv sync
make -C src CXX=g++ LDFLAGS=-lstdc++fs
scripts/bootstrap.sh
```

BinDash is required by the paired and supplemental pair-task experiments and
the RefSeq BinDash sketch benchmark.

## Running the Benchmarks

Run the paper-scale paired comparison from the repository root:

```bash
experiments/pair_task/jobs/paired_sketchsize.sh
```

Run all supplemental pair-task experiments with:

```bash
experiments/pair_task/jobs/supplementary.sh
```

Run the RefSeq OddSketch and BinDash sketch-build benchmarks with:

```bash
experiments/refseq_sketch_task/jobs/refseq_oddsketch_sketch.sh
experiments/refseq_sketch_task/jobs/refseq_bindash_sketch.sh
```

See each task's README for individual supplemental experiments, RefSeq dataset
acquisition, and runner options.

## Notes

- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository. The default helper script builds it from `https://github.com/zhaoxiaofei/bindash.git` at tag `v2.6`.
- Task configs define their default output roots via `paths.outdir`.
- Pair-task runs create a fresh timestamped directory and record the resolved
  settings alongside the observations and summaries.
- Exact-Jaccard helper binaries are built from `experiments/tools/src/` into `experiments/tools/bin/` by `make -C src`.
