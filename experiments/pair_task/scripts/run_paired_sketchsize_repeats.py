#!/usr/bin/env python3
"""Run a fresh paired OddSketch/BinDash sketch-size experiment."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

from common import allocate_run_dir, resolve_repo_root, resolve_task_root
from run_bindash_validation_repeats import (
    command_output,
    resolve_bindash,
    run_chunk,
)


EXACT_FIELDS = (
    "replicate",
    "genome_seed",
    "pair_id",
    "mutation_count",
    "genome_length",
    "kmerlen",
    "jaccard_true",
)

PAIRED_FIELDS = (
    "replicate",
    "genome_seed",
    "odd_hash_seed",
    "bindash_randseed",
    "pair_id",
    "mutation_count",
    "genome_length",
    "kmerlen",
    "sketch_size",
    "bindash_sketchsize64",
    "bindash_bbits",
    "bindash_payload_bits",
    "jaccard_true",
    "jaccard_oddsketch",
    "jaccard_bindash",
    "odd_raw_estimate",
    "hamming_distance",
    "num_buckets",
    "clipped",
)


def resolve_path(raw: str, base: Path) -> Path:
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def resolve_executable(raw: str | None, environment: str, default: Path) -> Path:
    candidate = os.environ.get(environment) or raw
    path = resolve_path(candidate, resolve_repo_root()) if candidate else default.resolve()
    if not path.exists() or not os.access(path, os.X_OK):
        raise SystemExit(f"required executable not found: {path}")
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_gzip_tsv_atomic(
    path: Path,
    rows: list[dict],
    fields: tuple[str, ...],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with gzip.open(temporary, "wt", newline="", compresslevel=6) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def parse_int(row: dict[str, str], field: str) -> int:
    return int(row[field])


def parse_finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"{field} is not finite")
    return value


def validate_exact_rows(
    rows: list[dict[str, str]],
    *,
    replicate: int,
    genome_seed: int,
    pairs_per_replicate: int,
    kmerlen: int,
) -> dict[int, dict[str, str]]:
    if len(rows) != pairs_per_replicate:
        raise ValueError(f"expected {pairs_per_replicate} exact rows, found {len(rows)}")
    if not rows or set(rows[0]) != set(EXACT_FIELDS):
        raise ValueError("exact intermediate file has an unexpected schema")

    keyed: dict[int, dict[str, str]] = {}
    for row in rows:
        if parse_int(row, "replicate") != replicate:
            raise ValueError("exact intermediate file replicate differs")
        if parse_int(row, "genome_seed") != genome_seed:
            raise ValueError("exact intermediate file genome seed differs")
        if parse_int(row, "kmerlen") != kmerlen:
            raise ValueError("exact intermediate file k-mer length differs")
        pair_id = parse_int(row, "pair_id")
        if pair_id in keyed:
            raise ValueError(f"duplicate pair ID in exact intermediate file: {pair_id}")
        true_value = parse_finite(row, "jaccard_true")
        if not 0.0 <= true_value <= 1.0:
            raise ValueError("true Jaccard is outside [0, 1]")
        parse_int(row, "mutation_count")
        parse_int(row, "genome_length")
        keyed[pair_id] = row
    if set(keyed) != set(range(1, pairs_per_replicate + 1)):
        raise ValueError("exact intermediate file pair IDs are incomplete")
    return keyed


def validate_paired_rows(
    rows: list[dict[str, str]],
    *,
    exact_by_id: dict[int, dict[str, str]],
    replicate: int,
    genome_seed: int,
    odd_hash_seed: int,
    bindash_randseed: int,
    sketch_size: int,
    bbits: int,
    pairs_per_replicate: int,
    kmerlen: int,
) -> None:
    if len(rows) != pairs_per_replicate:
        raise ValueError(f"expected {pairs_per_replicate} paired rows, found {len(rows)}")
    if not rows or set(rows[0]) != set(PAIRED_FIELDS):
        raise ValueError("paired intermediate file has an unexpected schema")
    sketchsize64 = sketch_size // (64 * bbits)
    payload_bits = sketchsize64 * 64 * bbits
    if payload_bits != sketch_size:
        raise ValueError("requested sketch size is not exactly representable by BinDash")

    seen: set[int] = set()
    for row in rows:
        expected_values = {
            "replicate": replicate,
            "genome_seed": genome_seed,
            "odd_hash_seed": odd_hash_seed,
            "bindash_randseed": bindash_randseed,
            "kmerlen": kmerlen,
            "sketch_size": sketch_size,
            "bindash_sketchsize64": sketchsize64,
            "bindash_bbits": bbits,
            "bindash_payload_bits": payload_bits,
        }
        for field, expected in expected_values.items():
            if parse_int(row, field) != expected:
                raise ValueError(f"{field} differs from the configured value")
        pair_id = parse_int(row, "pair_id")
        if pair_id in seen or pair_id not in exact_by_id:
            raise ValueError(f"unexpected or duplicate pair ID: {pair_id}")
        seen.add(pair_id)
        truth = parse_finite(row, "jaccard_true")
        expected_truth = parse_finite(exact_by_id[pair_id], "jaccard_true")
        if not math.isclose(truth, expected_truth, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"true Jaccard differs for pair ID {pair_id}")
        for field in ("jaccard_oddsketch", "jaccard_bindash"):
            value = parse_finite(row, field)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{field} is outside [0, 1]")
        parse_finite(row, "odd_raw_estimate")
        parse_int(row, "hamming_distance")
        parse_int(row, "num_buckets")
        if parse_int(row, "clipped") not in (0, 1):
            raise ValueError("clipped must be 0 or 1")
    if seen != set(exact_by_id):
        raise ValueError(
            "OddSketch and BinDash pair IDs do not match the exact intermediate file"
        )


def read_validated_exact_intermediate(
    path: Path,
    **validation: int,
) -> dict[int, dict[str, str]] | None:
    if not path.exists():
        return None
    try:
        rows = read_gzip_tsv(path)
        return validate_exact_rows(rows, **validation)
    except Exception as error:
        print(f"[intermediate-invalid] {path}: {error}", flush=True)
        return None


def read_validated_paired_intermediate(
    path: Path,
    **validation,
) -> list[dict[str, str]] | None:
    if not path.exists():
        return None
    try:
        rows = read_gzip_tsv(path)
        validate_paired_rows(rows, **validation)
        return rows
    except Exception as error:
        print(f"[intermediate-invalid] {path}: {error}", flush=True)
        return None


def run(command: list[str], *, stdout: Path | None = None) -> None:
    print("[run]", " ".join(command), flush=True)
    if stdout is None:
        subprocess.run(command, check=True)
        return
    with stdout.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)


def write_dataset_config(
    base_config: dict,
    path: Path,
    dataset_dir: Path,
    *,
    genome_seed: int,
    genome_length: int,
    pairs_per_replicate: int,
    mutation_min: int,
    mutation_max: int,
) -> None:
    config = json.loads(json.dumps(base_config))
    config.setdefault("paths", {})["outdir"] = str(dataset_dir)
    make = config.setdefault("make_genomes", {})
    make.update(
        {
            "seed_base": genome_seed,
            "genome_length": genome_length,
            "num_pairs": pairs_per_replicate,
            "mutation_min": mutation_min,
            "mutation_max": mutation_max,
        }
    )
    path.write_text(json.dumps(config, indent=2) + "\n")


def write_pairlist(path: Path, pairs: list[dict[str, str]]) -> None:
    with path.open("w") as handle:
        for pair in pairs:
            handle.write(f"{pair['file1']}\t{pair['file2']}\n")


def build_exact_rows(
    raw_rows: list[dict[str, str]],
    *,
    replicate: int,
    genome_seed: int,
    kmerlen: int,
) -> list[dict]:
    return [
        {
            "replicate": replicate,
            "genome_seed": genome_seed,
            "pair_id": row["pair_id"],
            "mutation_count": row["mutation_count"],
            "genome_length": row["genome_length"],
            "kmerlen": kmerlen,
            "jaccard_true": row["jaccard_true"],
        }
        for row in raw_rows
    ]


def verify_regenerated_pairs(
    pairs: list[dict[str, str]],
    exact_by_id: dict[int, dict[str, str]],
) -> None:
    keyed = {int(row["pair_id"]): row for row in pairs}
    if set(keyed) != set(exact_by_id):
        raise SystemExit("regenerated dataset pair IDs differ from the exact intermediate file")
    for pair_id, exact in exact_by_id.items():
        pair = keyed[pair_id]
        for field in ("mutation_count", "genome_length"):
            if int(pair[field]) != int(exact[field]):
                raise SystemExit(
                    f"regenerated dataset differs from exact intermediate file: "
                    f"pair={pair_id}, field={field}"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="experiments/pair_task/configs/sketchsize_repeats/config.json",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Create this directory for a new run instead of using paths.outdir",
    )
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--num-pairs", type=int, default=None)
    parser.add_argument("--genome-length", type=int, default=None)
    parser.add_argument("--mutation-min", type=int, default=None)
    parser.add_argument("--mutation-max", type=int, default=None)
    parser.add_argument("--sketch-sizes", default=None, help="Comma-separated payload bits")
    parser.add_argument("--bootstrap", type=int, default=None)
    parser.add_argument("--chunk-size", type=int, default=None)
    args = parser.parse_args()

    repo_root = resolve_repo_root()
    task_root = resolve_task_root()
    config_path = resolve_path(args.config, repo_root)
    config = json.loads(config_path.read_text())
    experiment = config["paired_sketchsize_repeats"]
    make = config["make_genomes"]
    bindash_config = experiment["bindash"]

    replicates = (
        args.replicates if args.replicates is not None else int(experiment["replicates"])
    )
    pairs_per_replicate = (
        args.num_pairs if args.num_pairs is not None else int(make["num_pairs"])
    )
    genome_length = (
        args.genome_length if args.genome_length is not None else int(make["genome_length"])
    )
    mutation_min = args.mutation_min if args.mutation_min is not None else int(make["mutation_min"])
    mutation_max = args.mutation_max if args.mutation_max is not None else int(make["mutation_max"])
    sketch_sizes = (
        [int(value) for value in args.sketch_sizes.split(",") if value.strip()]
        if args.sketch_sizes
        else [int(value) for value in experiment["sketch_sizes"]]
    )
    bootstrap = args.bootstrap if args.bootstrap is not None else int(experiment["bootstrap"])
    chunk_size = (
        args.chunk_size if args.chunk_size is not None else int(bindash_config["chunk_size"])
    )
    threads = int(os.environ.get("NSLOTS", experiment["threads"]))
    kmerlen = int(experiment["kmerlen"])
    j0 = float(experiment["j0"])
    pos_mode = str(experiment["pos_mode"])
    canonical = bool(experiment["canonical"])
    bbits = int(bindash_config["bbits"])
    minhashtype = int(bindash_config["minhashtype"])
    dens = int(bindash_config["dens"])
    genome_seed_base = int(experiment["genome_seed_base"])
    genome_seed_stride = int(experiment.get("genome_seed_stride", 1))
    odd_hash_seed_base = int(experiment["odd_hash_seed_base"])
    odd_hash_seed_stride = int(experiment.get("odd_hash_seed_stride", 1))
    bindash_randseed_base = int(experiment["bindash_randseed_base"])
    bindash_randseed_stride = int(experiment.get("bindash_randseed_stride", 1))
    if (
        replicates <= 0
        or pairs_per_replicate <= 0
        or genome_length <= 0
        or bootstrap <= 0
        or chunk_size <= 0
        or threads <= 0
        or genome_seed_stride < pairs_per_replicate
        or odd_hash_seed_stride <= 0
        or bindash_randseed_stride <= 0
    ):
        raise SystemExit(
            "replicates, pairs, genome length, bootstrap, chunk size, threads, and seed "
            "strides must be positive; genome_seed_stride must be at least pairs_per_replicate"
        )
    if mutation_min < 0 or mutation_max < mutation_min or mutation_max > genome_length:
        raise SystemExit("require 0 <= mutation_min <= mutation_max <= genome_length")
    if len(set(sketch_sizes)) != len(sketch_sizes) or not sketch_sizes:
        raise SystemExit("sketch sizes must be non-empty and unique")
    for sketch_size in sketch_sizes:
        if sketch_size <= 0 or sketch_size % (64 * bbits) != 0:
            raise SystemExit(
                f"sketch size {sketch_size} is not exactly representable with b={bbits}"
            )
    output_raw = Path(config["paths"]["outdir"])
    output_root = (
        output_raw.resolve() if output_raw.is_absolute() else (task_root / output_raw).resolve()
    )
    if args.output_dir:
        run_dir = resolve_path(args.output_dir, repo_root)
        if run_dir.exists():
            raise SystemExit(f"output directory already exists: {run_dir}")
        run_dir.mkdir(parents=True)
    else:
        output_root.mkdir(parents=True, exist_ok=True)
        run_dir = allocate_run_dir(output_root)
    work_root = run_dir / ".work"
    summary = run_dir / "summary"
    work_root.mkdir(parents=True)
    summary.mkdir(parents=True, exist_ok=True)

    runtime = {
        "replicates": replicates,
        "pairs_per_replicate": pairs_per_replicate,
        "genome_length": genome_length,
        "mutation_min": mutation_min,
        "mutation_max": mutation_max,
        "kmerlen": kmerlen,
        "sketch_sizes": sketch_sizes,
        "j0": j0,
        "pos_mode": pos_mode,
        "canonical": canonical,
        "threads": threads,
        "genome_seed_rule": "genome_seed_base + (replicate - 1) * genome_seed_stride",
        "genome_seed_base": genome_seed_base,
        "genome_seed_stride": genome_seed_stride,
        "odd_hash_seed_rule": "odd_hash_seed_base + (replicate - 1) * odd_hash_seed_stride",
        "odd_hash_seed_base": odd_hash_seed_base,
        "odd_hash_seed_stride": odd_hash_seed_stride,
        "bindash_randseed_rule": (
            "bindash_randseed_base + (replicate - 1) * bindash_randseed_stride"
        ),
        "bindash_randseed_base": bindash_randseed_base,
        "bindash_randseed_stride": bindash_randseed_stride,
        "bootstrap": bootstrap,
        "bootstrap_seed": int(experiment["bootstrap_seed"]),
        "high_jaccard_threshold": float(experiment["high_jaccard_threshold"]),
        "high_jaccard_operator": ">=",
        "jaccard_bins": [float(value) for value in experiment["jaccard_bins"]],
        "bindash_bbits": bbits,
        "bindash_chunk_size": chunk_size,
        "bindash_minhashtype": minhashtype,
        "bindash_dens": dens,
        "summary_figure": str(
            experiment.get(
                "summary_figure",
                "RMSE_by_true_jaccard_panels.png",
            )
        ),
    }
    if (
        Path(runtime["summary_figure"]).name != runtime["summary_figure"]
        or not runtime["summary_figure"].endswith(".png")
    ):
        raise SystemExit("summary_figure must be a PNG filename without a directory")
    resolved_config = json.loads(json.dumps(config))
    resolved_config["runtime"] = runtime
    used_config_path = run_dir / "used_config.json"
    used_config_path.write_text(json.dumps(resolved_config, indent=2) + "\n")

    oddsketch = resolve_executable(None, "ODDSKETCH_BIN", repo_root / "src" / "oddsketch")
    true_jaccard = resolve_executable(
        None,
        "TRUE_JACCARD_BIN",
        repo_root / "experiments" / "tools" / "bin" / "true_jaccard",
    )
    bindash_raw = os.environ.get("BINDASH_BIN") or str(
        bindash_config.get("bindash_bin", "")
    )
    bindash_candidate = resolve_path(bindash_raw, repo_root) if bindash_raw else None
    bindash = resolve_bindash(str(bindash_candidate) if bindash_candidate else None)
    make_genomes = task_root / "scripts" / "make_genomes.py"
    metadata_path = run_dir / "run_metadata.json"
    analyzer = (
        task_root
        / "analysis"
        / "aggregate"
        / "analyze_paired_sketchsize_repeats.py"
    )

    now = datetime.now().astimezone().isoformat()
    metadata = {
        "created_at": now,
        "status": "running",
        "config": str(config_path),
        "run_dir": str(run_dir),
        **runtime,
        "oddsketch_bin": str(oddsketch),
        "oddsketch_version": "not reported by executable",
        "oddsketch_sha256": sha256_file(oddsketch),
        "true_jaccard_bin": str(true_jaccard),
        "true_jaccard_sha256": sha256_file(true_jaccard),
        "bindash_bin": str(bindash),
        "bindash_version": command_output([str(bindash), "--version"]),
        "bindash_sha256": sha256_file(bindash),
        "rmse_definition": "sqrt(sum_r SSE_r / sum_r N_r)",
        "ci_definition": "paired bootstrap over replicate IDs",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    print(f"[paired-sketchsize] run_dir={run_dir}", flush=True)
    print(
        f"[paired-sketchsize] replicates={replicates} pairs={pairs_per_replicate} "
        f"sketch_sizes={','.join(str(value) for value in sketch_sizes)}",
        flush=True,
    )

    for replicate in range(1, replicates + 1):
        genome_seed = genome_seed_base + (replicate - 1) * genome_seed_stride
        odd_hash_seed = odd_hash_seed_base + (replicate - 1) * odd_hash_seed_stride
        bindash_randseed = (
            bindash_randseed_base + (replicate - 1) * bindash_randseed_stride
        )
        replicate_dir = work_root / f"replicate_{replicate:03d}"
        replicate_dir.mkdir(parents=True, exist_ok=True)
        exact_path = replicate_dir / "exact.tsv.gz"
        exact_validation = {
            "replicate": replicate,
            "genome_seed": genome_seed,
            "pairs_per_replicate": pairs_per_replicate,
            "kmerlen": kmerlen,
        }
        exact_by_id = read_validated_exact_intermediate(exact_path, **exact_validation)

        completed: dict[int, list[dict[str, str]]] = {}
        if exact_by_id is not None:
            for sketch_size in sketch_sizes:
                intermediate_path = replicate_dir / f"n{sketch_size}.tsv.gz"
                rows = read_validated_paired_intermediate(
                    intermediate_path,
                    exact_by_id=exact_by_id,
                    replicate=replicate,
                    genome_seed=genome_seed,
                    odd_hash_seed=odd_hash_seed,
                    bindash_randseed=bindash_randseed,
                    sketch_size=sketch_size,
                    bbits=bbits,
                    pairs_per_replicate=pairs_per_replicate,
                    kmerlen=kmerlen,
                )
                if rows is not None:
                    completed[sketch_size] = rows
                    print(
                        f"[intermediate-ok] replicate={replicate} n={sketch_size}",
                        flush=True,
                    )
        missing_sizes = [value for value in sketch_sizes if value not in completed]
        if exact_by_id is not None and not missing_sizes:
            continue

        temp_parent = Path(os.environ["TMPDIR"]) if os.environ.get("TMPDIR") else None
        with tempfile.TemporaryDirectory(
            prefix=f"paired-sketchsize-r{replicate:03d}-",
            dir=str(temp_parent) if temp_parent else None,
        ) as raw_work:
            work = Path(raw_work)
            dataset_dir = work / "dataset"
            dataset_config = work / "dataset_config.json"
            write_dataset_config(
                config,
                dataset_config,
                dataset_dir,
                genome_seed=genome_seed,
                genome_length=genome_length,
                pairs_per_replicate=pairs_per_replicate,
                mutation_min=mutation_min,
                mutation_max=mutation_max,
            )
            run([sys.executable, str(make_genomes), "--config", str(dataset_config)])
            pairs = read_tsv(dataset_dir / "pair_info.txt")
            if len(pairs) != pairs_per_replicate:
                raise SystemExit(
                    f"replicate={replicate}: expected {pairs_per_replicate} generated pairs, "
                    f"found {len(pairs)}"
                )
            pairlist = work / "fasta_pairs.tsv"
            write_pairlist(pairlist, pairs)

            if exact_by_id is None:
                raw_exact = work / "exact_raw.tsv"
                run(
                    [
                        str(true_jaccard),
                        f"--kmer={kmerlen}",
                        f"--pair-info={dataset_dir / 'pair_info.txt'}",
                        f"--out={raw_exact}",
                    ]
                )
                exact_rows = build_exact_rows(
                    read_tsv(raw_exact),
                    replicate=replicate,
                    genome_seed=genome_seed,
                    kmerlen=kmerlen,
                )
                exact_by_id = validate_exact_rows(exact_rows, **exact_validation)
                write_gzip_tsv_atomic(exact_path, exact_rows, EXACT_FIELDS)
                print(f"[intermediate-write] {exact_path}", flush=True)
            else:
                verify_regenerated_pairs(pairs, exact_by_id)

            command_log_root = work / "bindash_command_log"
            (command_log_root / "metadata").mkdir(parents=True)
            for sketch_size in missing_sizes:
                sketchsize64 = sketch_size // (64 * bbits)
                odd_output = work / f"odd_n{sketch_size}.tsv"
                odd_command = [
                    str(oddsketch),
                    "eval",
                    f"--pairlist={pairlist}",
                    f"--kmer={kmerlen}",
                    f"--sketch-size={sketch_size}",
                    f"--j0={j0}",
                    f"--pos-mode={pos_mode}",
                    f"--canonical={1 if canonical else 0}",
                    f"--hash-seed={odd_hash_seed}",
                    f"--threads={threads}",
                ]
                run(odd_command, stdout=odd_output)
                odd_rows = read_tsv(odd_output)
                if len(odd_rows) != pairs_per_replicate:
                    raise SystemExit(
                        f"replicate={replicate}, n={sketch_size}: "
                        f"OddSketch returned {len(odd_rows)} rows"
                    )

                bindash_estimates: dict[int, float] = {}
                for start in range(0, pairs_per_replicate, chunk_size):
                    chunk = pairs[start : start + chunk_size]
                    bindash_estimates.update(
                        run_chunk(
                            chunk,
                            bindash=bindash,
                            run_dir=command_log_root,
                            replicate=replicate,
                            sketch_size=sketch_size,
                            chunk_number=start // chunk_size + 1,
                            threads=threads,
                            kmerlen=kmerlen,
                            sketchsize64=sketchsize64,
                            bbits=bbits,
                            randseed=bindash_randseed,
                            minhashtype=minhashtype,
                            dens=dens,
                        )
                    )
                if set(bindash_estimates) != set(range(1, pairs_per_replicate + 1)):
                    raise SystemExit(
                        f"replicate={replicate}, n={sketch_size}: "
                        "BinDash did not return every requested pair"
                    )

                paired_rows = []
                for pair, odd in zip(pairs, odd_rows):
                    pair_id = int(pair["pair_id"])
                    if odd["file1"] != pair["file1"] or odd["file2"] != pair["file2"]:
                        raise SystemExit(
                            f"replicate={replicate}, n={sketch_size}, pair={pair_id}: "
                            "OddSketch pair order differs"
                        )
                    exact = exact_by_id[pair_id]
                    paired_rows.append(
                        {
                            "replicate": replicate,
                            "genome_seed": genome_seed,
                            "odd_hash_seed": odd_hash_seed,
                            "bindash_randseed": bindash_randseed,
                            "pair_id": pair_id,
                            "mutation_count": exact["mutation_count"],
                            "genome_length": exact["genome_length"],
                            "kmerlen": kmerlen,
                            "sketch_size": sketch_size,
                            "bindash_sketchsize64": sketchsize64,
                            "bindash_bbits": bbits,
                            "bindash_payload_bits": sketchsize64 * 64 * bbits,
                            "jaccard_true": exact["jaccard_true"],
                            "jaccard_oddsketch": odd["jaccard_oddsketch"],
                            "jaccard_bindash": f"{bindash_estimates[pair_id]:.17g}",
                            "odd_raw_estimate": odd["odd_raw_estimate"],
                            "hamming_distance": odd["hamming_distance"],
                            "num_buckets": odd["num_buckets"],
                            "clipped": odd["clipped"],
                        }
                    )
                validation = {
                    "exact_by_id": exact_by_id,
                    "replicate": replicate,
                    "genome_seed": genome_seed,
                    "odd_hash_seed": odd_hash_seed,
                    "bindash_randseed": bindash_randseed,
                    "sketch_size": sketch_size,
                    "bbits": bbits,
                    "pairs_per_replicate": pairs_per_replicate,
                    "kmerlen": kmerlen,
                }
                validate_paired_rows(paired_rows, **validation)
                intermediate_path = replicate_dir / f"n{sketch_size}.tsv.gz"
                write_gzip_tsv_atomic(intermediate_path, paired_rows, PAIRED_FIELDS)
                validated = read_validated_paired_intermediate(
                    intermediate_path, **validation
                )
                if validated is None:
                    raise SystemExit(
                        f"new intermediate file failed validation: {intermediate_path}"
                    )
                print(f"[intermediate-write] {intermediate_path}", flush=True)

    all_rows: list[dict[str, str]] = []
    for replicate in range(1, replicates + 1):
        replicate_dir = work_root / f"replicate_{replicate:03d}"
        exact_by_id = read_validated_exact_intermediate(
            replicate_dir / "exact.tsv.gz",
            replicate=replicate,
            genome_seed=genome_seed_base + (replicate - 1) * genome_seed_stride,
            pairs_per_replicate=pairs_per_replicate,
            kmerlen=kmerlen,
        )
        if exact_by_id is None:
            raise SystemExit(
                f"final exact intermediate validation failed for replicate={replicate}"
            )
        for sketch_size in sketch_sizes:
            intermediate_path = replicate_dir / f"n{sketch_size}.tsv.gz"
            rows = read_validated_paired_intermediate(
                intermediate_path,
                exact_by_id=exact_by_id,
                replicate=replicate,
                genome_seed=genome_seed_base + (replicate - 1) * genome_seed_stride,
                odd_hash_seed=(
                    odd_hash_seed_base + (replicate - 1) * odd_hash_seed_stride
                ),
                bindash_randseed=(
                    bindash_randseed_base
                    + (replicate - 1) * bindash_randseed_stride
                ),
                sketch_size=sketch_size,
                bbits=bbits,
                pairs_per_replicate=pairs_per_replicate,
                kmerlen=kmerlen,
            )
            if rows is None:
                raise SystemExit(
                    f"final paired intermediate validation failed: replicate={replicate}, "
                    f"n={sketch_size}"
                )
            all_rows.extend(rows)
    all_rows.sort(
        key=lambda row: (
            int(row["replicate"]),
            int(row["sketch_size"]),
            int(row["pair_id"]),
        )
    )
    paired_output = run_dir / "paired_observations.tsv.gz"
    write_gzip_tsv_atomic(paired_output, all_rows, PAIRED_FIELDS)

    metadata["paired_observation_count"] = len(all_rows)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    run(
        [
            sys.executable,
            str(analyzer),
            "--input",
            str(paired_output),
            "--metadata",
            str(metadata_path),
            "--outdir",
            str(summary),
        ]
    )
    metadata["status"] = "complete"
    metadata["completed_at"] = datetime.now().astimezone().isoformat()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    shutil.rmtree(work_root)
    print(f"[paired-sketchsize] complete: {run_dir}", flush=True)


if __name__ == "__main__":
    main()
