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

## Config
`config.json` controls the synthetic database generation, query generation, and search parameters.

- `genome_length`
  - Length of each synthetic genome in base pairs.
- `clusters.num_clusters`
  - Number of clusters used to generate the database genomes.
- `clusters.cluster_size`
  - Number of database genomes per cluster.
- `clusters.mutation_min`, `clusters.mutation_max`
  - Mutation-count range used when generating genomes inside each cluster.
- `clusters.seed`
  - Random seed for reproducible cluster generation.
- `query.num_queries`
  - Number of query genomes to generate.
- `query.query_mutation_min`, `query.query_mutation_max`
  - Mutation-count range used to derive query genomes from their source clusters.
- `paths.outdir`
  - Root directory for all generated files.
  - Default: `outputs/default`
- `oddsketch.kmerlen`
  - `k` used by OddSketch.
- `oddsketch.sketch_size`
  - Sketch size for OddSketch search.
- `oddsketch.j0`
  - OddSketch similarity threshold parameter.
- `oddsketch.pos_mode`
  - Positional sampling mode passed to OddSketch.
- `bindash.bindash_bin`
  - BinDash executable name or path.
- `bindash.kmerlen`, `bindash.sketch_size`, `bindash.bbits`
  - Main BinDash sketch parameters.
  - `sketch_size` is the requested number of bins.
  - Internally the script converts this to BinDash's `--sketchsize64`, so the effective size is rounded up to a multiple of 64 bins.
- `bindash.threads`
  - Number of threads used by BinDash.

Common edits:

- Change output location:
  - set `paths.outdir`
- Make a smaller smoke test:
  - reduce `genome_length`, `clusters.num_clusters`, `clusters.cluster_size`, and `query.num_queries`
- Compare search settings:
  - change `oddsketch.sketch_size`, `oddsketch.j0`, or `bindash.sketch_size`

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
