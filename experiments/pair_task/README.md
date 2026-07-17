# Pair Task

This task generates synthetic genome pairs and compares exact Jaccard, OddSketch, and BinDash.

This is the lighter benchmark workflow in this repository. The default
`config.json` generates 10 pairs of 1 Mbp synthetic genomes. The paper's
sketch-size experiment is substantially heavier: it evaluates eight sketch
sizes with 1,000 genome pairs per configuration.

BinDash is an external dependency and is not vendored in this repository. The default helper script installs it from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

with `BINDASH_TAG=v2.6`.

BinDash is required only for `cal_jaccard_bindash.py` or task runs with `bindash.enabled=true`. Set `bindash.enabled=false` for OddSketch-only runs.

For the benchmark baseline reported here, tag `v2.6` corresponds to commit `ce2d16816beade65db992b8cd6eced00b54ca9ef`, and the executable reports `version 2.2.0 commit ce2d168-clean`.

## Reproducing Paper Figures
The sketch-size workflow requires BinDash and sufficient compute time and
storage for eight configurations of 1,000 genome pairs each. From this
directory, run:

```bash
uv run python scripts/batch_project_runner.py --config-dir configs/sketchsize
```

After all configurations finish successfully, `batch_project_runner.py`
automatically creates `outputs/sketchsize/RMSEvsSKETCHSIZE.tsv` and regenerates
the sketch-size summary and RMSE-by-true-Jaccard figures. Only runs created by
the current batch are included, so earlier runs under the same output directory
do not affect the figures.

## Quick Start
Run the default 10-pair configuration as a local smoke test:

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

For a faster smoke test, reduce `make_genomes.genome_length` and
`make_genomes.num_pairs` in `config.json`.

`batch_project_runner.py` creates a unique run directory under the configured output root for each config, saves the resolved config to `<run>/metadata/used_config.json`, and generates figures for that run.
When you use the default config, it also updates `outputs/default/latest_used_config.json` so you can easily inspect the latest resolved config.

OddSketch in this task now runs in batch mode:
- `sketch`: one `oddsketch sketch` invocation over all unique genomes in `pair_info.txt`
- `dist`: one `oddsketch dist --pairlist=...` invocation over the generated sketch pairs

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
- `oddsketch.threads`
  - Number of threads used by OddSketch for both `sketch` and `dist`.
- `oddsketch.sketch_size`
  - Sketch size for OddSketch estimation.
- `oddsketch.j0`
  - OddSketch similarity threshold parameter.
- `oddsketch.pos_mode`
  - Positional sampling mode passed to OddSketch.
- `oddsketch.canonical`
  - Whether canonical k-mers are used.
  - The pair-task OddSketch script sketches each unique genome once and then evaluates only the requested pairs via `dist --pairlist`.
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
Default root: `outputs/default/<run>/`

- `genomes/`
- `pair_info.txt`, `genome_paths.txt`
- `results/jaccard_true_results.txt`
- `results/jaccard_oddsketch_results.txt`
- `results/jaccard_bindash_results.txt`
- `results/comparison_results_oddsketch.csv`
- `results/comparison_results_bindash.csv`
- `figures/`

## Manual Execution
Run the other bundled configuration groups with:

```bash
uv run python scripts/batch_project_runner.py --config-dir configs/default
uv run python scripts/batch_project_runner.py --config-dir configs/bbits
```

Recompute RMSE summaries from the result CSV files of a completed run:

```bash
uv run python analysis/per_run/compute_rmse.py \
  --csv outputs/default/<run>/results/comparison_results_oddsketch.csv \
  --csv outputs/default/<run>/results/comparison_results_bindash.csv
```

Inspect sketch storage for one completed sketch-size run:

```bash
uv run python analysis/per_run/report_sketch_memory.py \
  --run-dir outputs/sketchsize/<run>
```

The Grid Engine job script retained from the paper experiments can be submitted
with:

```bash
qsub jobs/qsub_project_runner.sh
```

Review the queue, resource, environment, and path settings in the job script
before using it on another cluster.

## Layout
- `config.json`: task settings
- `configs/`: configuration groups for the bundled experiments
- `scripts/`: genome generation, Jaccard calculation, and task runner
- `jobs/`: Grid Engine job script used for the paper experiments
- `analysis/per_run/`: per-run plotting, RMSE, and sketch-memory utilities
- `analysis/aggregate/`: summary plots across multiple runs
- `outputs/`: generated data, result tables, and figures
