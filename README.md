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
printf '%s\n' genome_001.fna.sketch genome_002.fna.sketch | src/oddsketch dist --threads=8
src/oddsketch dist --threads=8 --pairlist sketch_pairs.tsv
```

`--pairlist` expects a tab-separated file with one sketch pair per line.

## Pairwise Benchmark
Default outputs are written under `experiments/pair_task/outputs/default/`. You can override the root output directory in `experiments/pair_task/config.json` via `paths.outdir`, or pass `--outdir` to the genome generation step.

```bash
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
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
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

Generated files:
- FASTA pairs: `experiments/pair_task/outputs/default/genomes/`
- Pair metadata: `experiments/pair_task/outputs/default/pair_info.txt`
- Result tables: `experiments/pair_task/outputs/default/results/`
- Figures: `experiments/pair_task/outputs/default/figures/`

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

To download every assembly listed in `experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` first:

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
