# Pair Task

This task generates synthetic genome pairs and compares exact Jaccard, OddSketch, and BinDash.

## Layout
- `config.json`: task settings
- `scripts/`: genome generation, Jaccard calculation, and task runner
- `analysis/`: plotting, RMSE utilities, and figure generation
- `outputs/default/`: default generated data

## Quick Start
```bash
cd experiments/pair_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

Skip BinDash:

```bash
uv run python scripts/project_runner.py --config config.json --skip-bindash
uv run python analysis/make_figures.py
```

Step-by-step:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
uv run python analysis/make_figures.py
```

RMSE summary:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

## Outputs
Default root: `outputs/default/`

- `genomes/`
- `pair_info.txt`, `genome_paths.txt`
- `results/jaccard_true_results.txt`
- `results/jaccard_oddsketch_results.txt`
- `results/jaccard_bindash_results.txt`
- `results/comparison_results_oddsketch.csv`
- `results/comparison_results_bindash.csv`
- `figures/`
