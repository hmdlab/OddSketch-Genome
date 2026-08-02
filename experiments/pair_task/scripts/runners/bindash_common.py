#!/usr/bin/env python3
"""Shared BinDash command helpers for pair-task experiments."""

from __future__ import annotations

import csv
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from common import resolve_repo_root


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
        str(resolve_repo_root() / "experiments" / "tools" / "bin" / "bindash"),
        shutil.which("bindash") or "",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.exists() and os.access(path, os.X_OK):
            return path.resolve()
    raise SystemExit("BinDash executable not found; set BINDASH_BIN")


def parse_jaccard(raw: str) -> float:
    if "/" in raw:
        numerator, denominator = raw.split("/", 1)
        return float(numerator) / float(denominator)
    return float(raw)


def append_command_log(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
                [
                    str(replicate),
                    str(sketch_size),
                    str(chunk_number),
                    phase,
                    shlex.join(command),
                ],
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
