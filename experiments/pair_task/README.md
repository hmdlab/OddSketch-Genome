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
```

Step-by-step:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
```

`project_runner.py` creates a unique run directory under the configured output root, saves the resolved config to `<run>/metadata/used_config.json`, and generates figures for that run.
When you use the default config, it also updates `outputs/default/latest_used_config.json` so you can easily inspect the latest resolved config.

RMSE summary:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

## Config
`config.json` controls both the synthetic data generation and the estimation settings.

- `paths.outdir`
  - Root directory for all generated files.
  - Default: `outputs/default`
- `make_genomes.genome_length`
  - Length of each synthetic genome in base pairs.
- `make_genomes.num_pairs`
  - Number of genome pairs to generate.
- `make_genomes.mutation_min`, `make_genomes.mutation_max`
  - Minimum and maximum mutation counts applied to each pair.
  - Each pair is sampled within this range.
- `make_genomes.seed_base`
  - Base random seed for reproducible genome generation.
- `true_jaccard.kmerlen`
  - `k` used for exact Jaccard calculation.
- `oddsketch.kmerlen`
  - `k` used by OddSketch.
- `oddsketch.sketch_size`
  - Sketch size for OddSketch estimation.
- `oddsketch.j0`
  - OddSketch similarity threshold parameter.
- `oddsketch.pos_mode`
  - Positional sampling mode passed to OddSketch.
- `oddsketch.canonical`
  - Whether canonical k-mers are used.
- `bindash.bindash_bin`
  - BinDash executable name or path.
- `bindash.enabled`
  - Whether BinDash steps are executed in this task.
  - Set `false` to run OddSketch-only experiments.
- `bindash.threads`
  - Number of threads used by BinDash.
- `bindash.mode`
  - BinDash execution mode. The current default is `sketch_dist`.
- `bindash.kmerlen`, `bindash.sketch_size`, `bindash.bbits`
  - Main BinDash sketch parameters.
  - `sketch_size` is interpreted as the target total sketch size in bits.
  - Internally the script converts this to BinDash's `--sketchsize64` by dividing by `64 * bbits`, so the effective memory is rounded up to the nearest `64 * bbits` bits.
- `bindash.pair_cmd`
  - Command template used when evaluating one genome pair with BinDash.

Common edits:

- Change output location:
  - set `paths.outdir`
- Make a smaller smoke test:
  - reduce `make_genomes.genome_length` and `make_genomes.num_pairs`
- Compare different sketch settings:
  - change `oddsketch.sketch_size`, `oddsketch.j0`, or `bindash.sketch_size`

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
