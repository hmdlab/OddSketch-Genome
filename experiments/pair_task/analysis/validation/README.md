# Seeded validation workflows

These workflows generate independent synthetic datasets and independent OddSketch hash streams. Each qsub job runs the requested experiment, aggregates replicate-level statistics, and writes PNG, TSV, and Markdown outputs.

## Submit full experiments

Submit from the repository root (submission from `experiments/pair_task` is also supported):

```bash
qsub experiments/pair_task/jobs/qsub_validation_repeats.sh
qsub experiments/pair_task/jobs/qsub_k_sensitivity.sh
qsub experiments/pair_task/jobs/qsub_oph_baseline.sh
```

The defaults are defined in `experiments/pair_task/configs/validation/config.json`.

## Small cluster smoke tests

Arguments after the qsub script are passed to the Python runner:

```bash
qsub experiments/pair_task/jobs/qsub_validation_repeats.sh \
  --replicates 2 --num-pairs 20 --genome-length 20000 \
  --mutation-min 1 --mutation-max 100 --bootstrap 100
```

The same overrides work for the k-sensitivity and OPH jobs.

## Outputs

Each experiment creates a unique run below:

```text
outputs/validation/<experiment>/run_<timestamp>_<pid>/
```

`latest_run.txt` in each experiment directory records the latest run path. Every run contains:

- `metadata/used_config.json`: resolved design and runtime overrides
- `datasets/`: FASTA, exact-Jaccard, and raw seeded evaluation files
- `observations.tsv`: combined pair-level results
- `summary/report.md`: publication-oriented tables and method notes
- `summary/*.tsv`: machine-readable aggregate tables
- `summary/*.png`: final figures

The repeats workflow reports RMSE confidence intervals, clipping rates (`D >= n/2`), empty-bucket rates, and empirical 1%/5% saturation-safe boundaries. Confidence intervals resample independent replicates, not pooled individual rows.

The OPH workflow is memory matched at the payload level. For an `n`-bit OddSketch it retains `L=n/64` full 64-bit densified OPH minima, for exactly `n` OPH payload bits, and estimates Jaccard as `1 - S/L`. File-container and header overhead are excluded for both methods.

## Paired BinDash repeats

`qsub_bindash_validation_repeats.sh` reuses the saved 20-replicate OddSketch
datasets and evaluates BinDash on exactly the same 1,000 pairs in every
replicate and at all eight payload sizes. It writes checkpoint files per
replicate and sketch size, then produces paired RMSE confidence intervals,
RMSE-difference tables, clip-rate tables, PNG figures, and a Markdown report.

Submit from the repository root:

```bash
qsub experiments/pair_task/jobs/qsub_bindash_validation_repeats.sh
```

By default, the job reads the path in
`outputs/validation/repeats/latest_run.txt`. Pin a particular source run with:

```bash
qsub -v ODD_REPEATS_RUN=/absolute/path/to/run \
  experiments/pair_task/jobs/qsub_bindash_validation_repeats.sh
```

To resume an interrupted output run, also set `BINDASH_REPEATS_RUN_DIR` to its
absolute path. Outputs are under
`outputs/validation/bindash_repeats/<run>/summary/`.
