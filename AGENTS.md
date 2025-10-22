# Repository Guidelines

## Project Structure & Module Organization
- `src/`: core implementation and scripts.
  - `oddsketch.cpp` + `Makefile` → builds C++17/OpenMP binary `oddsketch`.
  - `oddpipe.py`: pipeline CLI (download → sketch → dist).
  - `extract_jaccard_values.py`, `cal/`, `test/`: analysis, plots, and sample data/scripts.
- Test data lives under `src/test/data/` (e.g., `test_genomes/`). Generated FASTA/sketch/figures/CSV are gitignored and should not be committed.

## Build, Test, and Development Commands
- Build core: `cd src && make` (produces `oddsketch` next to sources). Clean: `make clean`.
- Install Python deps: `pip install -r requirements.txt`.
- Pipeline examples:
  - Download: `cd src && ./oddpipe.py download --summary assembly_summary.txt --filter "reference genome" --limit 120 --out refgenomes.list`
  - Sketch: `cd src && ./oddpipe.py sketch --list refgenomes.list --threads 8 --out refgenomes.sketch`
  - Distances: `cd src && ./oddpipe.py dist --list refgenomes.sketch --threads 8 > refgenomes.dist.tsv`
- Optional plots: `cd src/test/analysis_images && python true_vs_oddsketch_compare.py`.

## Coding Style & Naming Conventions
- C++: C++17, compile with warnings (`-Wall -Wextra`). Prefer RAII. Functions/variables `snake_case`; types `CamelCase`. Indent 2–4 spaces. Headers in `src/`.
- Python: 3.8+, `snake_case`, 4‑space indent, rich `argparse` `--help`. Executables include `#!/usr/bin/env python3`.
- Files: keep large data out of repo; reference via `.list`/`.sketch` under `src/test/data/` when possible.

## Testing Guidelines
- This repo uses analysis scripts, not a unit test framework.
- Validate Jaccard outputs using `src/test/data/test_genomes/` as golden references.
- Reproduce: use the synthetic pipeline under `src/test` (generate → oddsketch → true → plots); report RMSE/R² and paste the first ~20 lines of key TSVs.
- Performance notes: state thread count (`--threads`) and machine info for comparisons.

## Commit & Pull Request Guidelines
- Commits: imperative, concise summaries; group related changes. Examples: `sketch and test complete`, `distpart complete`. Reference issues in body (e.g., `#123`).
- PRs: include purpose/background, exact reproduction commands, sample output (first ~20 TSV lines), and any figures under `src/images/` if helpful. Specify OS, compiler, and threads used.

## Security & Configuration Tips
- Requirements: OpenMP‑capable `g++`, Python 3, and sufficient disk for genome downloads.
- Use `src/test/pipeline_config.json` for all settings (including BinDash `bindash` section). Do not use a separate `bindash_config.json`.
- Do not commit large datasets; rely on `.list`/`.sketch` references. Test artifacts are already gitignored.
