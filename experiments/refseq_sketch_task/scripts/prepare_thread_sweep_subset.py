#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import shutil
from datetime import datetime, timezone
from pathlib import Path


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(raw: str) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (task_root() / path).resolve()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_paths(path: Path) -> list[Path]:
    result: list[Path] = []
    with path.open() as f:
        for line in f:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            genome = Path(raw).expanduser()
            if not genome.is_absolute():
                genome = (path.parent / genome).resolve()
            if not genome.exists():
                raise SystemExit(f"genome not found: {genome}")
            result.append(genome)
    return result


def write_paths(path: Path, values: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value.resolve()}\n" for value in values))


def decompress_gzip(gz_path: Path, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with gzip.open(gz_path, "rb") as inf, tmp.open("wb") as outf:
        shutil.copyfileobj(inf, outf)
    tmp.replace(out_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source-list",
        default="data/assembly/manifests/gzip_paths.txt",
        help="List of downloaded .fna.gz files.",
    )
    ap.add_argument(
        "--outdir",
        default="data/thread_sweep_1024",
        help="Output directory for the fixed FASTA subset.",
    )
    ap.add_argument("--count", type=int, default=1024)
    ap.add_argument("--seed", type=int, default=20260517)
    args = ap.parse_args()

    source_list = resolve_path(args.source_list)
    outdir = resolve_path(args.outdir)
    fasta_dir = outdir / "fasta"
    manifests_dir = outdir / "manifests"
    metadata_dir = outdir / "metadata"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    genomes = read_paths(source_list)
    if args.count <= 0:
        raise SystemExit("--count must be positive")
    if args.count > len(genomes):
        raise SystemExit(f"--count={args.count} exceeds available genomes={len(genomes)}")

    rng = random.Random(args.seed)
    selected_indices = sorted(rng.sample(range(len(genomes)), args.count))
    selected_gz = [genomes[i] for i in selected_indices]

    selected_fasta: list[Path] = []
    for idx, gz_path in enumerate(selected_gz, 1):
        if not gz_path.name.endswith(".gz"):
            raise SystemExit(f"expected .gz input: {gz_path}")
        fasta_path = fasta_dir / gz_path.name[:-3]
        decompress_gzip(gz_path, fasta_path)
        selected_fasta.append(fasta_path.resolve())
        if idx % 128 == 0 or idx == len(selected_gz):
            print(f"[prepare] decompressed {idx}/{len(selected_gz)}")

    write_paths(manifests_dir / "gzip_paths.txt", selected_gz)
    write_paths(manifests_dir / "fasta_paths.txt", selected_fasta)

    metadata = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_list": str(source_list),
        "source_list_sha256": sha256_file(source_list),
        "count": len(selected_gz),
        "seed": args.seed,
        "outdir": str(outdir),
        "gzip_manifest": str(manifests_dir / "gzip_paths.txt"),
        "fasta_manifest": str(manifests_dir / "fasta_paths.txt"),
        "total_gzip_bytes": sum(path.stat().st_size for path in selected_gz),
        "total_fasta_bytes": sum(path.stat().st_size for path in selected_fasta),
    }
    (metadata_dir / "subset_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"[done] subset={len(selected_fasta)} outdir={outdir}")


if __name__ == "__main__":
    main()
