#!/usr/bin/env python3
"""Run BinDash on the saved 20-replicate OddSketch validation datasets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


BINDASH_FIELDS = [
    "replicate",
    "genome_seed",
    "bindash_randseed",
    "pair_id",
    "mutation_count",
    "genome_length",
    "kmerlen",
    "sketch_size",
    "bindash_sketchsize64",
    "bindash_bbits",
    "jaccard_true",
    "jaccard_bindash",
    "file1",
    "file2",
]

PAIRED_FIELDS = [
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
    "jaccard_true",
    "jaccard_oddsketch",
    "jaccard_bindash",
    "odd_raw_estimate",
    "hamming_distance",
    "num_buckets",
    "clipped",
]


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def repo_root() -> Path:
    return task_root().parents[1]


def resolve_path(raw: str, base: Path) -> Path:
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str:
    process = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return process.stdout.strip()


def resolve_bindash(raw: str | None) -> Path:
    candidates = [
        raw or "",
        os.environ.get("BINDASH_BIN", ""),
        str(repo_root() / "experiments" / "tools" / "bin" / "bindash"),
        shutil.which("bindash") or "",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.exists() and os.access(path, os.X_OK):
            return path.resolve()
    raise SystemExit("BinDash executable not found; set BINDASH_BIN")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def allocate_run_dir(output_root: Path) -> Path:
    root = output_root / "bindash_repeats"
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    run_dir = root / f"{stamp}_{os.getpid()}"
    run_dir.mkdir()
    (root / "latest_run.txt").write_text(str(run_dir.resolve()) + "\n")
    return run_dir


def resolve_odd_run(raw: str | None, output_root: Path) -> Path:
    if raw:
        path = resolve_path(raw, repo_root())
    else:
        latest = output_root / "repeats" / "latest_run.txt"
        if not latest.exists():
            raise SystemExit(f"OddSketch latest-run pointer not found: {latest}")
        path = Path(latest.read_text().strip()).resolve()
    if not (path / "observations.tsv").exists():
        raise SystemExit(f"OddSketch observations not found: {path / 'observations.tsv'}")
    return path


def parse_jaccard(raw: str) -> float:
    if "/" in raw:
        numerator, denominator = raw.split("/", 1)
        return float(numerator) / float(denominator)
    return float(raw)


def append_command_log(path: Path, values: list[str]) -> None:
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if not exists:
            writer.writerow(["replicate", "sketch_size", "chunk", "phase", "command"])
        writer.writerow(values)


def run_command(command: list[str], *, capture: bool = False) -> str:
    process = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed ({process.returncode}): {shlex.join(command)}\n"
            f"stdout:\n{process.stdout or ''}\nstderr:\n{process.stderr or ''}"
        )
    return process.stdout or ""


def sketch_command(
    bindash: Path,
    input_list: Path,
    out_prefix: Path,
    *,
    threads: int,
    kmerlen: int,
    sketchsize64: int,
    bbits: int,
    randseed: int,
    minhashtype: int,
    dens: int,
) -> list[str]:
    return [
        str(bindash),
        "sketch",
        f"--listfname={input_list}",
        f"--nthreads={threads}",
        f"--kmerlen={kmerlen}",
        f"--sketchsize64={sketchsize64}",
        f"--bbits={bbits}",
        f"--randseed={randseed}",
        f"--minhashtype={minhashtype}",
        f"--dens={dens}",
        "--isstrandpreserved=false",
        "--iscasepreserved=false",
        f"--outfname={out_prefix}",
    ]


def run_chunk(
    pairs: list[dict[str, str]],
    *,
    bindash: Path,
    run_dir: Path,
    replicate: int,
    sketch_size: int,
    chunk_number: int,
    threads: int,
    kmerlen: int,
    sketchsize64: int,
    bbits: int,
    randseed: int,
    minhashtype: int,
    dens: int,
) -> dict[int, float]:
    temp_parent_raw = os.environ.get("TMPDIR")
    temp_parent = Path(temp_parent_raw) if temp_parent_raw else None
    with tempfile.TemporaryDirectory(
        prefix=f"bindash-r{replicate:03d}-n{sketch_size}-c{chunk_number:03d}-",
        dir=str(temp_parent) if temp_parent else None,
    ) as raw_work:
        work = Path(raw_work)
        query_list = work / "queries.txt"
        target_list = work / "targets.txt"
        query_list.write_text("".join(f"{row['file1']}\n" for row in pairs))
        target_list.write_text("".join(f"{row['file2']}\n" for row in pairs))
        query_prefix = work / "query_sketch"
        target_prefix = work / "target_sketch"

        query_cmd = sketch_command(
            bindash,
            query_list,
            query_prefix,
            threads=threads,
            kmerlen=kmerlen,
            sketchsize64=sketchsize64,
            bbits=bbits,
            randseed=randseed,
            minhashtype=minhashtype,
            dens=dens,
        )
        target_cmd = sketch_command(
            bindash,
            target_list,
            target_prefix,
            threads=threads,
            kmerlen=kmerlen,
            sketchsize64=sketchsize64,
            bbits=bbits,
            randseed=randseed,
            minhashtype=minhashtype,
            dens=dens,
        )
        dist_cmd = [
            str(bindash),
            "dist",
            f"--nthreads={threads}",
            "--ithres=0",
            "--mthres=1000000",
            "--pthres=1.0001",
            "--nneighbors=0",
            "--outfname=-",
            str(query_prefix),
            str(target_prefix),
        ]
        command_log = run_dir / "metadata" / "commands.tsv"
        for phase, command in (
            ("query_sketch", query_cmd),
            ("target_sketch", target_cmd),
            ("dist", dist_cmd),
        ):
            append_command_log(
                command_log,
                [str(replicate), str(sketch_size), str(chunk_number), phase, shlex.join(command)],
            )

        run_command(query_cmd)
        run_command(target_cmd)
        output = run_command(dist_cmd, capture=True)

    expected = {(row["file1"], row["file2"]): int(row["pair_id"]) for row in pairs}
    estimates: dict[int, float] = {}
    for line in output.splitlines():
        fields = line.rstrip().split("\t")
        if len(fields) < 5:
            continue
        pair_id = expected.get((fields[0], fields[1]))
        if pair_id is not None:
            estimates[pair_id] = parse_jaccard(fields[4])
    missing = sorted(set(expected.values()) - set(estimates))
    if missing:
        raise RuntimeError(
            f"BinDash omitted {len(missing)} requested pairs for replicate={replicate}, "
            f"sketch_size={sketch_size}, chunk={chunk_number}; first={missing[:5]}"
        )
    return estimates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/validation/config.json")
    parser.add_argument("--odd-run", default=None)
    parser.add_argument("--run-dir", default=None, help="Existing/new output directory; enables resume")
    parser.add_argument("--bindash-bin", default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--num-pairs", type=int, default=None)
    parser.add_argument("--sketch-sizes", default=None, help="Comma-separated payload bits")
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--bootstrap", type=int, default=None)
    args = parser.parse_args()

    config_path = resolve_path(args.config, task_root())
    config = json.loads(config_path.read_text())
    validation = config["validation"]
    repeats_cfg = validation["repeats"]
    bindash_cfg = validation["bindash_repeats"]
    output_raw = Path(config["paths"]["outdir"])
    output_root = output_raw.resolve() if output_raw.is_absolute() else (task_root() / output_raw).resolve()
    odd_run = resolve_odd_run(args.odd_run, output_root)

    requested_replicates = args.replicates or int(repeats_cfg["replicates"])
    pairs_per_replicate = args.num_pairs or int(config["make_genomes"]["num_pairs"])
    kmerlen = int(repeats_cfg["kmerlen"])
    sketch_sizes = (
        [int(value) for value in args.sketch_sizes.split(",") if value.strip()]
        if args.sketch_sizes
        else [int(value) for value in repeats_cfg["sketch_sizes"]]
    )
    bbits = int(bindash_cfg["bbits"])
    threads = int(os.environ.get("NSLOTS", bindash_cfg["threads"]))
    chunk_size = args.chunk_size or int(bindash_cfg["chunk_size"])
    randseed_base = int(bindash_cfg["randseed_base"])
    minhashtype = int(bindash_cfg["minhashtype"])
    dens = int(bindash_cfg["dens"])
    bootstrap = args.bootstrap or int(validation["bootstrap"])
    if requested_replicates <= 0 or pairs_per_replicate <= 0 or chunk_size <= 0 or threads <= 0:
        raise SystemExit("replicates, num-pairs, chunk-size, and threads must be positive")

    bindash = resolve_bindash(args.bindash_bin)
    if args.run_dir:
        run_dir = resolve_path(args.run_dir, repo_root())
        run_dir.mkdir(parents=True, exist_ok=True)
    else:
        run_dir = allocate_run_dir(output_root)
    for directory in (run_dir / "metadata", run_dir / "partials", run_dir / "summary"):
        directory.mkdir(parents=True, exist_ok=True)

    odd_rows = read_tsv(odd_run / "observations.tsv")
    odd_by_key = {
        (int(row["replicate"]), int(row["sketch_size"]), int(row["pair_id"])): row
        for row in odd_rows
        if int(row["replicate"]) <= requested_replicates
        and int(row["sketch_size"]) in sketch_sizes
        and int(row["kmerlen"]) == kmerlen
    }
    expected_odd = requested_replicates * len(sketch_sizes) * pairs_per_replicate
    if len(odd_by_key) != expected_odd:
        raise SystemExit(
            f"expected {expected_odd} OddSketch observations, found {len(odd_by_key)}; "
            "the saved run is incomplete or does not contain the requested pairs per replicate"
        )
    source_seeds = {}
    truth_by_pair = {}
    for (replicate, _, pair_id), row in odd_by_key.items():
        source_seeds.setdefault(
            replicate,
            (int(row["genome_seed"]), int(row["hash_seed"])),
        )
        truth_key = (replicate, pair_id)
        previous_truth = truth_by_pair.setdefault(truth_key, row["jaccard_true"])
        if previous_truth != row["jaccard_true"]:
            raise SystemExit(f"true Jaccard differs across sketch sizes for key={truth_key}")
    genome_seeds = [source_seeds[index][0] for index in sorted(source_seeds)]
    odd_hash_seeds = [source_seeds[index][1] for index in sorted(source_seeds)]
    if len(set(genome_seeds)) != requested_replicates:
        raise SystemExit("saved OddSketch datasets do not have unique genome seeds per replicate")
    if len(set(odd_hash_seeds)) != requested_replicates:
        raise SystemExit("saved OddSketch observations do not have unique hash seeds per replicate")

    metadata = {
        "created_at": datetime.now().astimezone().isoformat(),
        "source_odd_run": str(odd_run),
        "config": str(config_path),
        "replicates": requested_replicates,
        "pairs_per_replicate": pairs_per_replicate,
        "kmerlen": kmerlen,
        "sketch_sizes": sketch_sizes,
        "threads": threads,
        "chunk_size": chunk_size,
        "bindash_bin": str(bindash),
        "bindash_version": command_output([str(bindash), "--version"]),
        "bindash_sha256": sha256_file(bindash),
        "bbits": bbits,
        "randseed_rule": "randseed_base + replicate - 1",
        "randseed_base": randseed_base,
        "minhashtype": minhashtype,
        "dens": dens,
        "source_genome_seeds": genome_seeds,
        "source_odd_hash_seeds": odd_hash_seeds,
    }
    (run_dir / "metadata" / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"[bindash-repeats] run_dir={run_dir}")
    print(f"[bindash-repeats] source_odd_run={odd_run}")
    print(f"[bindash-repeats] version={metadata['bindash_version']}")
    print(f"[bindash-repeats] sha256={metadata['bindash_sha256']}")

    all_bindash_rows: list[dict] = []
    for replicate in range(1, requested_replicates + 1):
        pair_info = odd_run / "datasets" / f"replicate_{replicate:03d}" / "pair_info.txt"
        pairs = read_tsv(pair_info)
        if len(pairs) != pairs_per_replicate:
            raise SystemExit(
                f"expected {pairs_per_replicate} pairs in {pair_info}, found {len(pairs)}"
            )
        randseed = randseed_base + replicate - 1
        genome_seed = int(odd_by_key[(replicate, sketch_sizes[0], 1)]["genome_seed"])

        for sketch_size in sketch_sizes:
            if sketch_size % (64 * bbits) != 0:
                raise SystemExit(
                    f"sketch size {sketch_size} is not exactly representable with bbits={bbits}"
                )
            sketchsize64 = sketch_size // (64 * bbits)
            partial = run_dir / "partials" / f"replicate_{replicate:03d}_n{sketch_size}.tsv"
            if partial.exists():
                rows = read_tsv(partial)
                if len(rows) != pairs_per_replicate:
                    raise SystemExit(f"incomplete existing partial: {partial}")
                print(f"[resume] replicate={replicate} n={sketch_size}")
                all_bindash_rows.extend(rows)
                continue

            print(
                f"[run] replicate={replicate}/{requested_replicates} n={sketch_size} "
                f"sketchsize64={sketchsize64} randseed={randseed}",
                flush=True,
            )
            estimates: dict[int, float] = {}
            for start in range(0, len(pairs), chunk_size):
                chunk = pairs[start : start + chunk_size]
                chunk_number = start // chunk_size + 1
                estimates.update(
                    run_chunk(
                        chunk,
                        bindash=bindash,
                        run_dir=run_dir,
                        replicate=replicate,
                        sketch_size=sketch_size,
                        chunk_number=chunk_number,
                        threads=threads,
                        kmerlen=kmerlen,
                        sketchsize64=sketchsize64,
                        bbits=bbits,
                        randseed=randseed,
                        minhashtype=minhashtype,
                        dens=dens,
                    )
                )
            if len(estimates) != pairs_per_replicate:
                raise RuntimeError(
                    f"expected {pairs_per_replicate} BinDash estimates for replicate={replicate}, "
                    f"n={sketch_size}; found {len(estimates)}"
                )

            rows = []
            for pair in pairs:
                pair_id = int(pair["pair_id"])
                odd = odd_by_key[(replicate, sketch_size, pair_id)]
                rows.append(
                    {
                        "replicate": replicate,
                        "genome_seed": genome_seed,
                        "bindash_randseed": randseed,
                        "pair_id": pair_id,
                        "mutation_count": pair["mutation_count"],
                        "genome_length": pair["genome_length"],
                        "kmerlen": kmerlen,
                        "sketch_size": sketch_size,
                        "bindash_sketchsize64": sketchsize64,
                        "bindash_bbits": bbits,
                        "jaccard_true": odd["jaccard_true"],
                        "jaccard_bindash": f"{estimates[pair_id]:.17g}",
                        "file1": pair["file1"],
                        "file2": pair["file2"],
                    }
                )
            write_tsv(partial, rows, BINDASH_FIELDS)
            all_bindash_rows.extend(rows)

    all_bindash_rows.sort(
        key=lambda row: (int(row["replicate"]), int(row["sketch_size"]), int(row["pair_id"]))
    )
    bindash_output = run_dir / "bindash_observations.tsv"
    write_tsv(bindash_output, all_bindash_rows, BINDASH_FIELDS)

    paired_rows = []
    for bindash_row in all_bindash_rows:
        key = (
            int(bindash_row["replicate"]),
            int(bindash_row["sketch_size"]),
            int(bindash_row["pair_id"]),
        )
        odd = odd_by_key[key]
        paired_rows.append(
            {
                "replicate": key[0],
                "genome_seed": odd["genome_seed"],
                "odd_hash_seed": odd["hash_seed"],
                "bindash_randseed": bindash_row["bindash_randseed"],
                "pair_id": key[2],
                "mutation_count": odd["mutation_count"],
                "genome_length": odd["genome_length"],
                "kmerlen": odd["kmerlen"],
                "sketch_size": key[1],
                "bindash_sketchsize64": bindash_row["bindash_sketchsize64"],
                "bindash_bbits": bindash_row["bindash_bbits"],
                "jaccard_true": odd["jaccard_true"],
                "jaccard_oddsketch": odd["jaccard_oddsketch"],
                "jaccard_bindash": bindash_row["jaccard_bindash"],
                "odd_raw_estimate": odd["odd_raw_estimate"],
                "hamming_distance": odd["hamming_distance"],
                "num_buckets": odd["num_buckets"],
                "clipped": odd["clipped"],
            }
        )
    paired_output = run_dir / "paired_observations.tsv"
    write_tsv(paired_output, paired_rows, PAIRED_FIELDS)

    analyzer = task_root() / "analysis" / "validation" / "analyze_paired_bindash.py"
    bins = ",".join(str(value) for value in validation["bins"])
    analyze_command = [
        sys.executable,
        str(analyzer),
        "--input",
        str(paired_output),
        "--metadata",
        str(run_dir / "metadata" / "run_metadata.json"),
        "--outdir",
        str(run_dir / "summary"),
        "--bins",
        bins,
        "--bootstrap",
        str(bootstrap),
    ]
    print(f"[analyze] {shlex.join(analyze_command)}", flush=True)
    subprocess.run(analyze_command, check=True)
    print(f"[complete] {run_dir}")


if __name__ == "__main__":
    main()
