#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import tempfile
from pathlib import Path

from common import load_config, resolve_config_path, resolve_output_root, resolve_task_root


def resolve_oddsketch_bin() -> str:
    task_root = resolve_task_root()
    repo_root = task_root.parents[1]
    candidates = [
        os.environ.get("ODDSKETCH_BIN", ""),
        str(repo_root / "src" / "oddsketch"),
        "oddsketch",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    return "oddsketch"


def build_oddsketch_base_cmd(subcommand: str, cfg: dict) -> list[str]:
    oddcfg = cfg.get("oddsketch", {}) if isinstance(cfg, dict) else {}
    cmd = [
        resolve_oddsketch_bin(),
        subcommand,
        f"--kmer={oddcfg.get('kmerlen', 64)}",
        f"--sketch-size={oddcfg.get('sketch_size', 8192)}",
        f"--j0={oddcfg.get('j0', 0.75)}",
        f"--pos-mode={oddcfg.get('pos_mode', 'mix')}",
        f"--threads={oddcfg.get('threads', 1)}",
    ]
    canonical = oddcfg.get("canonical")
    if canonical is not None:
        cmd.append(f"--canonical={1 if canonical else 0}")
    return cmd


def run_oddsketch_sketch(genome_files: list[str], cfg: dict) -> list[str]:
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".list") as f:
        for genome_file in genome_files:
            f.write(f"{genome_file}\n")
        temp_input_paths = f.name

    try:
        cmd = build_oddsketch_base_cmd("sketch", cfg)
        cmd.append(f"--input-paths={temp_input_paths}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    except subprocess.CalledProcessError:
        return []
    finally:
        os.unlink(temp_input_paths)


def run_oddsketch_dist_pairlist(pairlist_path: Path, cfg: dict) -> list[dict]:
    cmd = build_oddsketch_base_cmd("dist", cfg)
    cmd.append(f"--pairlist={pairlist_path}")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        return []

    distances = []
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 3:
            file1, file2, jaccard_dist = parts
            distances.append({
                "file1": file1,
                "file2": file2,
                "jaccard_estimate": float(jaccard_dist),
            })
    return distances


def read_pair_info(pair_info_file: Path) -> list[dict]:
    pairs = []
    with pair_info_file.open("r") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            pair_id, file1, file2, mutation_count, genome_length = parts
            pairs.append({
                "pair_id": int(pair_id),
                "file1": file1,
                "file2": file2,
                "mutation_count": int(mutation_count),
                "genome_length": int(genome_length),
            })
    return pairs


def unique_genome_files(pairs: list[dict]) -> list[str]:
    seen = set()
    genomes = []
    for pair in pairs:
        for key in ("file1", "file2"):
            path = pair[key]
            if path not in seen:
                seen.add(path)
                genomes.append(path)
    return genomes


def write_pairlist_file(pairlist_path: Path, pairs: list[dict], sketch_map: dict[str, str]) -> None:
    with pairlist_path.open("w") as f:
        for pair in pairs:
            f.write(f"{sketch_map[pair['file1']]}\t{sketch_map[pair['file2']]}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json", help="Path to task config JSON")
    args = ap.parse_args()

    task_root = resolve_task_root()
    config_path = resolve_config_path(args.config)
    cfg = load_config(config_path)
    output_root = resolve_output_root(task_root, cfg)
    pair_info_file = output_root / "pair_info.txt"
    results_dir = output_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_file = results_dir / "jaccard_oddsketch_results.txt"
    true_path = results_dir / "jaccard_true_results.txt"
    out_csv = results_dir / "comparison_results_oddsketch.csv"

    if not pair_info_file.exists():
        raise SystemExit(f"pair_info not found: {pair_info_file}")

    print(f"Config file: {config_path}")
    print(f"Input: {pair_info_file}")
    print(f"Result output: {output_file}")

    pairs = read_pair_info(pair_info_file)
    genome_files = unique_genome_files(pairs)
    sketch_files = run_oddsketch_sketch(genome_files, cfg)
    if len(sketch_files) != len(genome_files):
        raise SystemExit(
            f"oddsketch sketch failed or returned an unexpected number of files: "
            f"expected {len(genome_files)}, got {len(sketch_files)}"
        )

    sketch_map = dict(zip(genome_files, sketch_files))
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv") as f:
        pairlist_path = Path(f.name)
    try:
        write_pairlist_file(pairlist_path, pairs, sketch_map)
        distances = run_oddsketch_dist_pairlist(pairlist_path, cfg)
    finally:
        os.unlink(pairlist_path)

    if len(distances) != len(pairs):
        raise SystemExit(
            f"oddsketch dist failed or returned an unexpected number of rows: "
            f"expected {len(pairs)}, got {len(distances)}"
        )

    results = []
    for pair, distance in zip(pairs, distances):
        results.append({
            "pair_id": pair["pair_id"],
            "mutation_count": pair["mutation_count"],
            "genome_length": pair["genome_length"],
            "jaccard_estimate": distance["jaccard_estimate"],
            "sketch_file1": distance["file1"],
            "sketch_file2": distance["file2"],
        })

    with output_file.open("w") as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_estimate\tsketch_file1\tsketch_file2\n")
        for result in results:
            f.write(
                f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                f"{result['jaccard_estimate']:.10f}\t{result['sketch_file1']}\t{result['sketch_file2']}\n"
            )

    if true_path.exists():
        true = {}
        with true_path.open() as tf:
            rd = csv.reader(tf, delimiter="\t")
            next(rd, None)
            for row in rd:
                if row:
                    true[int(row[0])] = {
                        "mutation_count": int(row[1]),
                        "jaccard_true": float(row[4]),
                    }
        with out_csv.open("w") as cf:
            writer = csv.writer(cf)
            writer.writerow(["pair_id", "mutation_count", "jaccard_true", "jaccard_oddsketch"])
            for result in sorted(results, key=lambda item: item["pair_id"]):
                truth = true.get(result["pair_id"])
                if truth:
                    writer.writerow([
                        result["pair_id"],
                        truth["mutation_count"],
                        truth["jaccard_true"],
                        result["jaccard_estimate"],
                    ])


if __name__ == "__main__":
    main()
