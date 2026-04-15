# Experiments Repro Package

This directory is the reproducibility package for benchmarking **OddSketch vs BinDash**.

Japanese version: `README-ja.md`

- Core implementation lives in repository root (`src/`, `include/`).
- Experimental scripts/configs are isolated here.
- Generated benchmark outputs should go to `experiments/outputs/`.

## Important: BinDash Handling

This repository does **not** vendor BinDash source code or binaries.

- BinDash is an external research tool by its original authors.
- For reproducibility we install/build BinDash in `experiments/scripts/build_tools.sh`.
- If needed, pin a specific reference with `--ref <tag-or-commit>`.

## Directory Layout

- `env/Dockerfile`: containerized benchmark environment
- `scripts/download_data.sh`: placeholder for optional dataset download
- `scripts/build_tools.sh`: install/build external tools (BinDash)
- `scripts/run_benchmark.sh`: run the benchmark pipeline
- `scripts/make_figures.py`: generate comparison plots
- `configs/`: benchmark config files
- `outputs/`: generated artifacts (gitignored)

## Quick Run (Host)

From repository root:

```bash
uv sync
bash experiments/scripts/build_tools.sh --method auto
bash experiments/scripts/run_benchmark.sh --mode all
```

Optional flags:

```bash
bash experiments/scripts/run_benchmark.sh --mode pair
bash experiments/scripts/run_benchmark.sh --mode search
bash experiments/scripts/run_benchmark.sh --mode all --skip-bindash
```

## Docker Repro (Recommended)

From repository root:

```bash
docker build -f experiments/env/Dockerfile -t oddsketch-exp .
docker run --rm -it oddsketch-exp
```

Mount host outputs if needed:

```bash
docker run --rm -it \
  -v "$(pwd)/experiments/outputs:/workspace/experiments/outputs" \
  oddsketch-exp
```

## Notes

- Pairwise benchmark scripts use `experiments/pair_task/data/test_genomes/` as working data.
- Search benchmark outputs are configured via `configs/search_task.config.json` and default to `experiments/outputs/search_task`.
- To use a specific OddSketch binary, set `ODDSKETCH_BIN` before running.

## Citation

If you use BinDash in comparisons, cite the original BinDash papers from the BinDash project.
