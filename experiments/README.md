# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure
- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `search_task/`: nearest-neighbor search benchmark on clustered synthetic genomes
- `scripts/run_benchmark.sh`: runs one or both tasks
- `scripts/make_figures.py`: regenerates figures from each task output directory
- `env/Dockerfile`: container environment

Each task keeps its own defaults:
- config: `<task>/config.json`
- code: `<task>/scripts/`, `<task>/analysis/`
- generated data: `<task>/outputs/default/`

## Quick Run
From repository root:

```bash
uv sync
bash experiments/scripts/run_benchmark.sh --mode all
```

Optional:

```bash
bash experiments/scripts/run_benchmark.sh --mode pair
bash experiments/scripts/run_benchmark.sh --mode search
bash experiments/scripts/run_benchmark.sh --mode all --skip-bindash
```

## Notes
- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository.
- Task-local `config.json` files define the default output roots via `paths.outdir`.
- Exact-Jaccard helper binaries are built into `tools/bin/` by `make -C src`.
