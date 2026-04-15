# Pair Task

This task generates synthetic genome pairs and compares exact Jaccard, OddSketch, and BinDash.

## Layout
- `config.json`: task settings
- `scripts/`: genome generation and Jaccard calculation
- `analysis/`: plotting and RMSE utilities
- `outputs/default/`: default generated data

## Quick Start
```bash
cd experiments/pair_task
python scripts/make_genomes.py --config config.json
python scripts/cal_jaccard_true.py --config config.json
python scripts/cal_jaccard_oddsketch.py --config config.json
python scripts/cal_jaccard_bindash.py --config config.json
python analysis/plot_true_vs_oddsketch.py
python analysis/plot_true_vs_bindash.py
```

RMSE summary:

```bash
python analysis/compute_rmse.py \
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
