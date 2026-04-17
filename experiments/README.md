# Experiments

This directory contains benchmark workflows built around the standalone OddSketch implementation in `src/`.

## Structure
- `pair_task/`: pairwise Jaccard benchmark on synthetic genome pairs
- `search_task/`: nearest-neighbor search benchmark on clustered synthetic genomes
- `tools/`: C++ helper tools and external-tool setup scripts used only by experiment workflows
- `env/Dockerfile`: container environment

Each task keeps its own defaults:
- config: `<task>/config.json`
- code: `<task>/scripts/`, `<task>/analysis/`
- generated data: `<task>/outputs/default/`

## Quick Run
Run each task from its own directory:

```bash
uv sync
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
```

Or:

```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

## Notes
- To use a specific OddSketch binary, set `ODDSKETCH_BIN`.
- BinDash is external and is not vendored in this repository.
- Task-local `config.json` files define the default output roots via `paths.outdir`.
- `project_runner.py` creates a unique run directory under `paths.outdir` and saves the resolved config to `<run>/metadata/used_config.json`.
- Exact-Jaccard helper binaries are built from `experiments/tools/src/` into `experiments/tools/bin/` by `make -C src`.
