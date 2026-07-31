# Seeded validation workflows

These workflows generate independent synthetic datasets and independent OddSketch hash streams for the k-mer sensitivity and memory-matched OPH experiments. Each qsub job runs the requested experiment, aggregates replicate-level statistics, and writes PNG, TSV, and Markdown outputs.

## Submit full experiments

Submit all supplemental experiments from the repository root:

```bash
qsub experiments/pair_task/jobs/supplementary.sh
```

Omitting the experiment argument runs all supplemental experiments.
The defaults are defined in `experiments/pair_task/configs/validation/config.json`.
The main 20-replicate OddSketch/BinDash comparison, including RMSE confidence
intervals and OddSketch clipping rates, is run by
`experiments/pair_task/jobs/paired_sketchsize.sh`.

Submit one supplemental experiment with:

```bash
qsub experiments/pair_task/jobs/supplementary.sh bindash_recommended
qsub experiments/pair_task/jobs/supplementary.sh k_sensitivity
qsub experiments/pair_task/jobs/supplementary.sh oph
```

## Small cluster smoke tests

Arguments after the experiment name are passed to its Python runner:

```bash
qsub experiments/pair_task/jobs/supplementary.sh k_sensitivity \
  --replicates 2 --num-pairs 20 --genome-length 20000 \
  --mutation-min 1 --mutation-max 100 --bootstrap 100
```

The same overrides work for the other supplemental experiments.

## Outputs

Each submission creates a unique parent run below:

```text
outputs/validation/run_<timestamp>_<job-id>/
├── bindash_recommended/
├── k_sensitivity/
└── oph/
```

An individual submission contains only the selected experiment directory.
Runs are not resumed or updated after creation. The main and BinDash workflows
use a temporary `.work/` directory during computation and remove it after
successful completion.

Each k-mer sensitivity and OPH experiment contains:

- `metadata/used_config.json`: resolved design and runtime overrides
- `datasets/`: FASTA, exact-Jaccard, and raw seeded evaluation files
- `observations.tsv`: combined pair-level results
- `summary/report.md`: publication-oriented tables and method notes
- `summary/*.tsv`: machine-readable aggregate tables
- `summary/*.png`: final figures

The OPH workflow is memory matched at the payload level. For an `n`-bit OddSketch it retains `L=n/64` full 64-bit densified OPH minima, for exactly `n` OPH payload bits, and estimates Jaccard as `1 - S/L`. File-container and header overhead are excluded for both methods.
