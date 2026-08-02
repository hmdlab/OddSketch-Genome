#!/usr/bin/env python3
"""Run seeded OddSketch validation experiments and generate final figures/tables."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from common import resolve_repo_root, resolve_task_root


OBSERVATION_FIELDS = (
    "experiment",
    "replicate",
    "genome_seed",
    "hash_seed",
    "pair_id",
    "mutation_count",
    "genome_length",
    "kmerlen",
    "sketch_size",
    "j0",
    "jaccard_true",
    "jaccard_oddsketch",
    "jaccard_oph",
    "odd_raw_estimate",
    "hamming_distance",
    "num_buckets",
    "oph_mismatches",
    "empty_buckets_left",
    "empty_buckets_right",
    "oph_num_buckets",
    "oph_storage_bits",
    "oph_empty_buckets_left",
    "oph_empty_buckets_right",
    "clipped",
)


def resolve_config(raw: str) -> Path:
    task_root = resolve_task_root()
    candidates = (
        Path(raw),
        task_root / raw,
        task_root / "configs" / "supplementary.json",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise SystemExit(f"validation config not found: {raw}")


def resolve_binary(env_name: str, default: Path) -> Path:
    raw = os.environ.get(env_name)
    path = Path(raw).expanduser() if raw else default
    if not path.is_absolute():
        path = (resolve_repo_root() / path).resolve()
    if not path.exists() or not os.access(path, os.X_OK):
        raise SystemExit(f"required executable not found: {path}")
    return path


def allocate_run_dir(base: Path, experiment: str) -> Path:
    root = base
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    parent = root / f"{stamp}_{os.getpid()}"
    run_dir = parent / experiment
    run_dir.mkdir(parents=True)
    (root / "latest_run.txt").write_text(f"{parent}\n")
    return run_dir


def run(cmd: list[str], *, stdout_path: Path | None = None) -> None:
    print("[run]", " ".join(str(part) for part in cmd), flush=True)
    if stdout_path is None:
        subprocess.run(cmd, check=True)
        return
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w") as output:
        subprocess.run(cmd, check=True, stdout=output)


def write_dataset_config(
    base_cfg: dict,
    dataset_dir: Path,
    genome_seed: int,
    genome_length: int,
    num_pairs: int,
    mutation_min: int,
    mutation_max: int,
) -> Path:
    cfg = json.loads(json.dumps(base_cfg))
    cfg.setdefault("paths", {})["outdir"] = str(dataset_dir)
    make_cfg = cfg.setdefault("make_genomes", {})
    make_cfg["seed_base"] = genome_seed
    make_cfg["genome_length"] = genome_length
    make_cfg["num_pairs"] = num_pairs
    make_cfg["mutation_min"] = mutation_min
    make_cfg["mutation_max"] = mutation_max
    metadata_dir = dataset_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    path = metadata_dir / "used_config.json"
    path.write_text(json.dumps(cfg, indent=2) + "\n")
    return path


def read_pair_info(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fasta_pairlist(path: Path, pairs: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for pair in pairs:
            handle.write(f"{pair['file1']}\t{pair['file2']}\n")


def keyed_rows(path: Path, key: str = "pair_id") -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        return {int(row[key]): row for row in csv.DictReader(handle, delimiter="\t")}


def eval_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def experiment_matrix(config: dict, experiment: str) -> tuple[int, list[int], list[int], int]:
    section = config[experiment]
    replicates = int(section["replicates"])
    seed_base = int(section["hash_seed_base"])
    if experiment == "k_sensitivity":
        k_values = [int(value) for value in section["k_values"]]
        sketch_sizes = [int(section["sketch_size"])]
    else:
        k_values = [int(section["kmerlen"])]
        sketch_sizes = [int(value) for value in section["sketch_sizes"]]
    return replicates, k_values, sketch_sizes, seed_base


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True, choices=("k_sensitivity", "oph"))
    parser.add_argument("--config", default="configs/supplementary.json")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Create this directory for the selected experiment",
    )
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--num-pairs", type=int, default=None)
    parser.add_argument("--genome-length", type=int, default=None)
    parser.add_argument("--mutation-min", type=int, default=None)
    parser.add_argument("--mutation-max", type=int, default=None)
    parser.add_argument("--bootstrap", type=int, default=None)
    args = parser.parse_args()

    task_root = resolve_task_root()
    repo_root = resolve_repo_root()
    config_path = resolve_config(args.config)
    cfg = json.loads(config_path.read_text())
    validation_cfg = cfg["validation"]
    analysis_cfg = cfg["analysis"]
    replicates, k_values, sketch_sizes, hash_seed_base = experiment_matrix(
        cfg, args.experiment
    )
    if args.replicates is not None:
        replicates = args.replicates
    if replicates <= 0:
        raise SystemExit("replicates must be positive")

    make_cfg = cfg["make_genomes"]
    num_pairs = args.num_pairs if args.num_pairs is not None else int(make_cfg["num_pairs"])
    genome_length = (
        args.genome_length if args.genome_length is not None else int(make_cfg["genome_length"])
    )
    mutation_min = args.mutation_min if args.mutation_min is not None else int(make_cfg["mutation_min"])
    mutation_max = args.mutation_max if args.mutation_max is not None else int(make_cfg["mutation_max"])
    if mutation_min < 0 or mutation_max < mutation_min or mutation_max > genome_length:
        raise SystemExit("require 0 <= mutation_min <= mutation_max <= genome_length")
    genome_seed_base = int(make_cfg["seed_base"])
    threads = int(validation_cfg.get("threads", 1))
    j0 = float(validation_cfg["j0"])
    pos_mode = str(validation_cfg.get("pos_mode", "mix"))
    canonical = bool(validation_cfg.get("canonical", True))
    experiment_cfg = cfg[args.experiment]
    memory_matched_oph = bool(experiment_cfg.get("memory_matched", False))
    bootstrap = args.bootstrap if args.bootstrap is not None else int(analysis_cfg["bootstrap"])

    output_base_raw = Path(cfg["paths"]["outdir"])
    output_base = output_base_raw if output_base_raw.is_absolute() else (task_root / output_base_raw)
    if args.output_dir:
        output_path = Path(args.output_dir).expanduser()
        run_dir = (
            output_path.resolve()
            if output_path.is_absolute()
            else (repo_root / output_path).resolve()
        )
        if run_dir.exists():
            raise SystemExit(f"output directory already exists: {run_dir}")
        run_dir.mkdir(parents=True)
    else:
        run_dir = allocate_run_dir(output_base.resolve(), args.experiment)
    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    resolved = json.loads(json.dumps(cfg))
    resolved["runtime"] = {
        "experiment": args.experiment,
        "replicates": replicates,
        "num_pairs": num_pairs,
        "genome_length": genome_length,
        "mutation_min": mutation_min,
        "mutation_max": mutation_max,
        "bootstrap": bootstrap,
        "k_values": k_values,
        "sketch_sizes": sketch_sizes,
        "memory_matched_oph": memory_matched_oph,
    }
    (metadata_dir / "used_config.json").write_text(json.dumps(resolved, indent=2) + "\n")

    oddsketch_bin = resolve_binary("ODDSKETCH_BIN", repo_root / "src" / "oddsketch")
    true_bin = resolve_binary(
        "TRUE_JACCARD_BIN", repo_root / "experiments" / "tools" / "bin" / "true_jaccard"
    )
    make_genomes_script = task_root / "scripts" / "runners" / "make_genomes.py"
    observations_path = run_dir / "observations.tsv"

    print(f"[validation] run_dir={run_dir}")
    print(f"[validation] experiment={args.experiment} replicates={replicates}")
    print(f"[validation] k_values={k_values} sketch_sizes={sketch_sizes}")

    with observations_path.open("w", newline="") as observations_handle:
        writer = csv.DictWriter(observations_handle, fieldnames=OBSERVATION_FIELDS, delimiter="\t")
        writer.writeheader()

        for replicate in range(1, replicates + 1):
            genome_seed = genome_seed_base + replicate * 100000
            hash_seed = hash_seed_base + replicate
            dataset_dir = run_dir / "datasets" / f"replicate_{replicate:03d}"
            dataset_cfg = write_dataset_config(
                cfg,
                dataset_dir,
                genome_seed,
                genome_length,
                num_pairs,
                mutation_min,
                mutation_max,
            )
            print(f"\n=== Replicate {replicate}/{replicates} ===", flush=True)
            run([sys.executable, str(make_genomes_script), "--config", str(dataset_cfg)])

            pair_info_path = dataset_dir / "pair_info.txt"
            pairs = read_pair_info(pair_info_path)
            if len(pairs) != num_pairs:
                raise SystemExit(
                    f"pair count mismatch in {pair_info_path}: expected {num_pairs}, got {len(pairs)}"
                )
            pairlist_path = dataset_dir / "fasta_pairs.tsv"
            write_fasta_pairlist(pairlist_path, pairs)

            for kmerlen in k_values:
                truth_path = dataset_dir / "truth" / f"true_jaccard_k{kmerlen}.tsv"
                truth_path.parent.mkdir(parents=True, exist_ok=True)
                run([
                    str(true_bin),
                    f"--kmer={kmerlen}",
                    f"--pair-info={pair_info_path}",
                    f"--out={truth_path}",
                ])
                truths = keyed_rows(truth_path)

                for sketch_size in sketch_sizes:
                    eval_path = (
                        dataset_dir
                        / "raw"
                        / f"eval_k{kmerlen}_n{sketch_size}_seed{hash_seed}.tsv"
                    )
                    command = [
                        str(oddsketch_bin),
                        "eval",
                        f"--pairlist={pairlist_path}",
                        f"--kmer={kmerlen}",
                        f"--sketch-size={sketch_size}",
                        f"--j0={j0}",
                        f"--pos-mode={pos_mode}",
                        f"--canonical={1 if canonical else 0}",
                        f"--hash-seed={hash_seed}",
                        f"--threads={threads}",
                    ]
                    if memory_matched_oph:
                        if sketch_size % 64 != 0:
                            raise SystemExit(
                                f"memory-matched 64-bit OPH requires sketch_size divisible by 64: "
                                f"{sketch_size}"
                            )
                        oph_buckets = sketch_size // 64
                        if oph_buckets <= 0 or (oph_buckets & (oph_buckets - 1)) != 0:
                            raise SystemExit(
                                f"memory-matched OPH bucket count must be a power of two: "
                                f"n={sketch_size}, L={oph_buckets}"
                            )
                        command.append(f"--oph-buckets={oph_buckets}")
                    run(command, stdout_path=eval_path)
                    estimates = eval_rows(eval_path)
                    if len(estimates) != len(pairs):
                        raise SystemExit(
                            f"eval row count mismatch: expected {len(pairs)}, got {len(estimates)}"
                        )

                    for pair, estimate in zip(pairs, estimates):
                        pair_id = int(pair["pair_id"])
                        truth = truths.get(pair_id)
                        if truth is None:
                            raise SystemExit(f"missing truth for pair {pair_id} in {truth_path}")
                        if estimate["file1"] != pair["file1"] or estimate["file2"] != pair["file2"]:
                            raise SystemExit(f"pair ordering mismatch for pair {pair_id}")
                        writer.writerow({
                            "experiment": args.experiment,
                            "replicate": replicate,
                            "genome_seed": genome_seed,
                            "hash_seed": hash_seed,
                            "pair_id": pair_id,
                            "mutation_count": pair["mutation_count"],
                            "genome_length": pair["genome_length"],
                            "kmerlen": kmerlen,
                            "sketch_size": sketch_size,
                            "j0": j0,
                            "jaccard_true": truth["jaccard_true"],
                            "jaccard_oddsketch": estimate["jaccard_oddsketch"],
                            "jaccard_oph": estimate["jaccard_oph"],
                            "odd_raw_estimate": estimate["odd_raw_estimate"],
                            "hamming_distance": estimate["hamming_distance"],
                            "num_buckets": estimate["num_buckets"],
                            "oph_mismatches": estimate["oph_mismatches"],
                            "empty_buckets_left": estimate["empty_buckets_left"],
                            "empty_buckets_right": estimate["empty_buckets_right"],
                            "oph_num_buckets": estimate["oph_num_buckets"],
                            "oph_storage_bits": estimate["oph_storage_bits"],
                            "oph_empty_buckets_left": estimate["oph_empty_buckets_left"],
                            "oph_empty_buckets_right": estimate["oph_empty_buckets_right"],
                            "clipped": estimate["clipped"],
                        })
                    observations_handle.flush()

    analysis_script = task_root / "scripts" / "analysis" / "validation.py"
    bins = ",".join(str(value) for value in analysis_cfg["bins"])
    run([
        sys.executable,
        str(analysis_script),
        "--experiment",
        args.experiment,
        "--input",
        str(observations_path),
        "--outdir",
        str(run_dir / "summary"),
        "--bins",
        bins,
        "--bootstrap",
        str(bootstrap),
    ])
    print(f"[validation] complete: {run_dir}")


if __name__ == "__main__":
    main()
