#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter


ASSEMBLY_COLUMNS = [
    "assembly_accession",
    "bioproject",
    "biosample",
    "wgs_master",
    "refseq_category",
    "taxid",
    "species_taxid",
    "organism_name",
    "infraspecific_name",
    "isolate",
    "version_status",
    "assembly_level",
    "release_type",
    "genome_rep",
    "seq_rel_date",
    "asm_name",
    "submitter",
    "gbrs_paired_asm",
    "paired_asm_comp",
    "ftp_path",
    "excluded_from_refseq",
    "relation_to_type_material",
]


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def repo_root() -> Path:
    return task_root().parents[1]


def resolve_config(path_arg: str) -> Path:
    candidates = [Path(path_arg), task_root() / path_arg, task_root() / "config.json"]
    for path in candidates:
        if path.exists():
            return path.resolve()
    raise SystemExit(f"config not found: {path_arg}")


def resolve_path(base: Path, raw: str | None) -> Path | None:
    if raw is None:
        return None
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def allocate_run_dir(data_root: Path, run_id: str | None) -> Path:
    if run_id:
        return data_root / "runs" / run_id
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    candidate = data_root / "runs" / stamp
    suffix = 1
    while candidate.exists():
        candidate = data_root / "runs" / f"{stamp}_{suffix:02d}"
        suffix += 1
    return candidate


def refseq_sketch_cfg(cfg: dict) -> dict:
    return cfg.get("refseq_sketch", {})


def copy_assembly_summary(cfg: dict, metadata_dir: Path) -> Path:
    paths = cfg.get("paths", {})
    source = resolve_path(task_root(), paths.get("assembly_summary"))
    saved = metadata_dir / "assembly_summary_refseq.txt"

    if not source or not source.exists():
        raise SystemExit("assembly_summary not found. Set paths.assembly_summary.")
    shutil.copy2(source, saved)
    return saved


def assembly_comments(path: Path) -> list[str]:
    comments: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
            elif line.strip():
                break
    return comments


def read_assembly_summary(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(ASSEMBLY_COLUMNS):
                continue
            rows.append(dict(zip(ASSEMBLY_COLUMNS, parts)))
    return rows


def select_assemblies(summary_path: Path, cfg: dict) -> list[dict[str, str]]:
    limit = refseq_sketch_cfg(cfg).get("limit")
    excluded = {
        str(value)
        for value in cfg.get("download", {}).get("excluded_accessions", [])
    }
    rows = [
        row
        for row in read_assembly_summary(summary_path)
        if row.get("ftp_path", "") not in ("", "na")
        and row.get("assembly_accession", "") not in excluded
    ]
    if limit is not None:
        rows = rows[: int(limit)]
    expected = cfg.get("download", {}).get("expected_genome_count")
    if limit is None and expected is not None and len(rows) != int(expected):
        raise SystemExit(
            f"selected assembly count mismatch: expected={expected}, actual={len(rows)}"
        )
    return rows


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def read_local_genome_list(list_path: Path) -> list[tuple[str, Path]]:
    result: list[tuple[str, Path]] = []
    with list_path.open() as f:
        for idx, line in enumerate(f, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            path = Path(raw).expanduser()
            if not path.is_absolute():
                path = (list_path.parent / path).resolve()
            if not path.exists():
                raise SystemExit(f"genome not found in local_genome_list: {path}")
            result.append((f"local_{idx:06d}", path))
    return result


def is_gzip_path(path: Path) -> bool:
    return path.name.endswith(".gz")


def ungzip_name(path: Path) -> str:
    return path.name[:-3] if path.name.endswith(".gz") else path.name


def make_input_links(genomes: list[tuple[str, Path]], input_dir: Path) -> list[Path]:
    input_dir.mkdir(parents=True, exist_ok=True)
    linked: list[Path] = []
    for idx, (name, source) in enumerate(genomes, 1):
        if is_gzip_path(source):
            link = input_dir / f"{idx:08d}_{safe_name(ungzip_name(source))}.gz"
        else:
            suffix = "".join(source.suffixes) or ".fna"
            link = input_dir / f"{idx:06d}_{safe_name(name)}{suffix}"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(source)
        linked.append(link)
    return linked


def write_list(path: Path, values: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value.absolute()}\n" for value in values))


def resolve_oddsketch_bin() -> str:
    candidates = [
        os.environ.get("ODDSKETCH_BIN", ""),
        str(repo_root() / "src" / "oddsketch"),
        "oddsketch",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    return "oddsketch"


def sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        p = subprocess.run(
            ["git", "-C", str(repo_root()), "rev-parse", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return p.stdout.strip()
    except Exception:
        return None


def oddsketch_cmd(
    cfg: dict,
    input_list: Path | None = None,
    out_dir: Path | None = None,
    sketch_paths_out: Path | None = None,
    skip_existing: bool = False,
) -> list[str]:
    odd = cfg.get("oddsketch", {})
    cmd = [
        resolve_oddsketch_bin(),
        "sketch",
        f"--kmer={odd.get('kmerlen', 64)}",
        f"--sketch-size={odd.get('sketch_size', 8192)}",
        f"--j0={odd.get('j0', 0.9)}",
        f"--pos-mode={odd.get('pos_mode', 'mix')}",
        f"--canonical={1 if odd.get('canonical', True) else 0}",
        f"--threads={odd.get('threads', 1)}",
    ]
    if input_list is not None:
        cmd.extend(["--input-paths", str(input_list)])
    if out_dir is not None:
        cmd.extend(["--out-dir", str(out_dir)])
    if sketch_paths_out is not None:
        cmd.extend(["--sketch-paths-out", str(sketch_paths_out)])
    if skip_existing:
        cmd.append("--skip-existing")
    return cmd


def parse_time_v(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    if not path.exists():
        return metrics
    for line in path.read_text(errors="replace").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        metrics[key.strip()] = value.strip()
    return metrics


def run_sketch(
    input_list: Path,
    stdout_path: Path,
    time_path: Path,
    cfg: dict,
    out_dir: Path,
    sketch_paths_out: Path,
    skip_existing: bool = False,
) -> tuple[float, dict[str, str]]:
    cmd = oddsketch_cmd(
        cfg,
        input_list=input_list,
        out_dir=out_dir,
        sketch_paths_out=sketch_paths_out,
        skip_existing=skip_existing,
    )
    timed_cmd = ["/usr/bin/time", "-v", "-o", str(time_path), *cmd]
    if not Path("/usr/bin/time").exists():
        timed_cmd = cmd
    print("[run]", " ".join(timed_cmd), flush=True)
    t0 = perf_counter()
    with stdout_path.open("w") as stdout:
        subprocess.run(timed_cmd, stdout=stdout, check=True)
    elapsed = perf_counter() - t0
    return elapsed, parse_time_v(time_path)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--resume", action="store_true", help="Resume an existing --run-id from completed batches.")
    args = ap.parse_args()

    cfg_path = resolve_config(args.config)
    cfg = json.loads(cfg_path.read_text())
    data_root = resolve_path(task_root(), cfg.get("paths", {}).get("data_root")) or (task_root() / "outputs")
    if args.resume and not args.run_id:
        raise SystemExit("--resume requires --run-id")
    run_dir = allocate_run_dir(data_root, args.run_id)
    if args.resume and not run_dir.exists():
        raise SystemExit(f"resume run not found: {run_dir}")
    if not args.resume and args.run_id and run_dir.exists():
        raise SystemExit(f"run already exists: {run_dir}")
    metadata_dir = run_dir / "metadata"
    manifests_dir = run_dir / "manifests"
    input_dir = run_dir / "genome_inputs"
    sketches_dir = run_dir / "sketches"
    logs_dir = run_dir / "logs"
    for d in (metadata_dir, manifests_dir, input_dir, sketches_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    used_config = metadata_dir / "used_config.json"
    used_config.write_text(json.dumps(cfg, indent=2) + "\n")

    summary_path = copy_assembly_summary(cfg, metadata_dir)
    selected_rows = select_assemblies(summary_path, cfg)
    write_tsv(metadata_dir / "selected_assemblies.tsv", selected_rows, ASSEMBLY_COLUMNS)

    local_list = resolve_path(task_root(), cfg.get("paths", {}).get("local_genome_list"))
    if not local_list:
        raise SystemExit("paths.local_genome_list is required for sketch. Run the download job first.")
    genomes = read_local_genome_list(local_list)
    limit = refseq_sketch_cfg(cfg).get("limit")
    if limit is not None:
        genomes = genomes[: int(limit)]

    input_list = manifests_dir / "genome_paths.txt"
    original_inputs = [path for _, path in genomes]
    uses_gzip_inputs = any(is_gzip_path(path) for path in original_inputs)
    linked_inputs = make_input_links(genomes, input_dir)
    write_list(input_list, linked_inputs)

    metadata = {
        "run_dir": str(run_dir),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "refseq_sketch_version_label": refseq_sketch_cfg(cfg).get("version_label"),
        "assembly_summary_source": cfg.get("paths", {}).get("assembly_summary"),
        "assembly_summary_saved": str(summary_path),
        "assembly_summary_comments": assembly_comments(summary_path),
        "selected_assemblies": len(selected_rows),
        "input_genomes": len(linked_inputs),
        "uses_gzip_inputs": uses_gzip_inputs,
        "git_commit": git_commit(),
        "oddsketch_bin": resolve_oddsketch_bin(),
        "oddsketch_bin_sha256": sha256_file(Path(resolve_oddsketch_bin())),
        "oddsketch_command": oddsketch_cmd(
            cfg,
            input_list=input_list,
            out_dir=sketches_dir,
            sketch_paths_out=manifests_dir / "sketch_paths.txt",
            skip_existing=args.resume,
        ),
    }
    (metadata_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    if args.prepare_only:
        print(f"[prepare] wrote {input_list}")
        print(f"[prepare] run_dir={run_dir}")
        return

    workflow_t0 = perf_counter()
    stdout_path = logs_dir / "oddsketch_sketch_stdout.txt"
    time_path = logs_dir / "oddsketch_sketch_time.txt"
    sketch_list = manifests_dir / "sketch_paths.txt"
    elapsed_sec, time_metrics = run_sketch(
        input_list,
        stdout_path,
        time_path,
        cfg,
        out_dir=sketches_dir,
        sketch_paths_out=sketch_list,
        skip_existing=args.resume,
    )
    relocated = read_local_genome_list(sketch_list)
    relocated_paths = [path for _, path in relocated]
    batch_count = 1
    workflow_elapsed_sec = perf_counter() - workflow_t0

    total_sketch_bytes = sum(path.stat().st_size for path in relocated_paths)
    timing = {
        "elapsed_sec": f"{elapsed_sec:.6f}",
        "workflow_elapsed_sec": f"{workflow_elapsed_sec:.6f}",
        "temporary_decompression_sec": "",
        "max_rss_kbytes": time_metrics.get("Maximum resident set size (kbytes)", ""),
        "user_sec": time_metrics.get("User time (seconds)", ""),
        "system_sec": time_metrics.get("System time (seconds)", ""),
        "cpu_percent": time_metrics.get("Percent of CPU this job got", ""),
        "exit_status": time_metrics.get("Exit status", ""),
        "input_genomes": str(len(linked_inputs)),
        "sketch_count": str(len(relocated_paths)),
        "batch_count": str(batch_count),
        "temporary_decompression": "False",
        "total_sketch_bytes": str(total_sketch_bytes),
    }
    write_tsv(run_dir / "results" / "oddsketch_sketch_metrics.tsv", [timing], list(timing.keys()))
    print(f"[done] sketches={len(relocated_paths)} bytes={total_sketch_bytes}")
    print(f"[done] metrics={run_dir / 'results' / 'oddsketch_sketch_metrics.tsv'}")
    print(f"[done] run_dir={run_dir}")


if __name__ == "__main__":
    main()
