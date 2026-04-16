# Search Task

This task generates clustered synthetic genomes, builds a genome DB, and compares nearest-neighbor search with OddSketch and BinDash.

## Layout
- `config.json`: task settings
- `scripts/`: generation, exact labels, search, evaluation, runners
- `analysis/`: plotting and figure generation
- `outputs/default/`: default generated data

## Quick Start
```bash
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
uv run python analysis/make_figures.py
```

End-to-end:

```bash
uv run python scripts/project_runner.py --config config.json
uv run python scripts/project_runner.py --config config.json --skip-bindash
uv run python analysis/make_figures.py
```

Repeated runs:

```bash
uv run python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234
uv run python scripts/repeat_runner.py --config config.json --runs 10 --seed-base 1234 --skip-bindash
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
uv run python analysis/make_figures.py
```
