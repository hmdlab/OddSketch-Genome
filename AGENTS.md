# Repository Guidelines

## Project Structure & Module Organization
- `src/`: core implementation.
  - `oddsketch.cpp` + `Makefile`: builds C++17/OpenMP binary `oddsketch`.
  - `oddpipe.py`: pipeline CLI (download → sketch → dist).
  - `extract_jaccard_values.py`, `cal/`, `test/`: analysis, plotting, and sample data/scripts.
- Test data: `src/test/data/` (e.g., `test_genomes/`). Generated FASTA/sketch/figures/CSV are gitignored.

## Build, Test, and Development Commands
- Build core: `cd src && make` (produces `oddsketch` next to sources). Clean: `make clean`.
- Pipeline (example):
  - Download: `cd src && ./oddpipe.py download --summary assembly_summary.txt --filter "reference genome" --limit 120 --out refgenomes.list`
  - Sketch: `cd src && ./oddpipe.py sketch --list refgenomes.list --threads 8 --out refgenomes.sketch`
  - Distances: `cd src && ./oddpipe.py dist --list refgenomes.sketch --threads 8 > refgenomes.dist.tsv`
- Python deps: `pip install -r requirements.txt`.
- Optional plots: `cd src/test/analysis_images && python true_vs_oddsketch_compare.py`.

## Coding Style & Naming Conventions
- C++: C++17, compile with warnings (`-Wall -Wextra`). Prefer RAII; functions/variables `snake_case`, types `CamelCase`; indent 2–4 spaces; headers live in `src/`.
- Python: 3.8+, `snake_case`, 4‑space indent, rich `argparse` `--help`. Executables include `#!/usr/bin/env python3`.
- Files: large data stays out of repo; reference via text lists (`.list`, `.sketch`) under `src/test/data/` when possible.

## Testing Guidelines
- Approach: repository uses analysis scripts, not a unit test framework.
- Golden data: validate Jaccard outputs against `src/test/data/test_genomes/` results.
- Reproduce: regenerate `.sketch` lists, then run `src/test/text_bindash_oddsketch_compare.py`; report RMSE/R² and attach the first lines of key TSVs.
- Performance: state thread count and machine info for comparisons.

## Commit & Pull Request Guidelines
- Commits: imperative, concise summary; group related changes. Examples: `sketch and test complete`, `distpart complete`. Reference issues in body (e.g., `#123`).
- PRs: include purpose/background, exact reproduction commands, sample output (first ~20 TSV lines), and figures under `src/images/` if helpful. Specify OS/compiler and threads used.

## Security & Configuration Tips
- Requirements: OpenMP‑capable `g++`, Python 3, sufficient disk for genome downloads.
- Do not commit large datasets; rely on `.list`/`.sketch` references. Test artifacts are already gitignored.

