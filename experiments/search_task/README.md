# Search Task

This task generates clustered synthetic genomes, builds a genome DB, and compares nearest-neighbor search with OddSketch and BinDash.

## Layout
- `config.json`: task settings
- `scripts/`: generation, exact labels, search, evaluation, runners
- `analysis/`: plotting
- `outputs/default/`: default generated data

## Quick Start
```bash
cd experiments/search_task
python scripts/make_cluster_query_genomes.py --config config.json
python scripts/true_db.py --config config.json
python scripts/oddsketch_db.py --config config.json
python scripts/bindash_db.py --config config.json
python scripts/evaluate_nn.py --config config.json
```

End-to-end:

```bash
python scripts/project_runner.py --config config.json
```

Repeated runs:

```bash
python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234
```

## Outputs
Default root: `outputs/default/`

- `genomes/`, `queries/`
- `db_genomes.list`, `queries.list`
- `true_pairs.tsv`, `true_nn.tsv`
- `oddsketch_pairs.tsv`, `oddsketch_nn.tsv`
- `bindash_pairs.tsv`, `bindash_nn.tsv`
- `nn_eval.tsv`
- `figures/`

Example plot:

```bash
python analysis/plot_est_vs_true.py \
  --true outputs/default/true_pairs.tsv \
  --pred outputs/default/oddsketch_pairs.tsv \
  --pred-col jaccard_oddsketch \
  --out outputs/default/figures/oddsketch_true_vs_estimate.png
```
