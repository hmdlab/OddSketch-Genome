#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
from time import perf_counter

from refseq_sketch_runner import (
    ASSEMBLY_COLUMNS,
    allocate_run_dir,
    assembly_comments,
    copy_assembly_summary,
    git_commit,
    is_gzip_path,
    make_input_links,
    read_local_genome_list,
    refseq_sketch_cfg,
    repo_root,
    resolve_config,
    resolve_path,
    select_assemblies,
    sha256_file,
    task_root,
    write_list,
    write_tsv,
)


def resolve_bindash_bin(raw: str | None = None) -> str:
    candidates = [
        os.environ.get("BINDASH_BIN", ""),
        raw or "",
        str(repo_root() / "experiments" / "tools" / "bin" / "bindash"),
        "bindash",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.is_absolute() or path.parent != Path("."):
            if path.exists() and os.access(path, os.X_OK):
                return str(path.resolve())
            continue
        repo_relative = repo_root() / candidate
        if repo_relative.exists() and os.access(repo_relative, os.X_OK):
            return str(repo_relative.resolve())
        found = shutil.which(candidate)
        if found:
            return found
    return raw or os.environ.get("BINDASH_BIN", "") or "bindash"


def bindash_cfg(cfg: dict) -> dict:
    return cfg.get("bindash", {}) if isinstance(cfg.get("bindash"), dict) else {}


def resolve_bindash_sketch_params(cfg: dict) -> tuple[int, int, int]:
    odd = cfg.get("oddsketch", {}) if isinstance(cfg.get("oddsketch"), dict) else {}
    b = bindash_cfg(cfg)
    bbits = int(b.get("bbits", 16))
    if bbits <= 0:
        raise SystemExit("bindash.bbits must be positive")
    if "sketch_size" in b:
        target_bits = int(b["sketch_size"])
    elif "sketchsize64" in b:
        target_bits = 64 * int(b["sketchsize64"]) * bbits
    else:
        target_bits = int(odd.get("sketch_size", 16384))
    if target_bits <= 0:
        raise SystemExit("bindash.sketch_size must be positive")
    sketchsize64 = max(1, math.ceil(target_bits / (64 * bbits)))
    effective_bits = 64 * sketchsize64 * bbits
    return sketchsize64, effective_bits, bbits


def bindash_cmd(cfg: dict, input_list: Path, out_prefix: Path) -> list[str]:
    odd = cfg.get("oddsketch", {}) if isinstance(cfg.get("oddsketch"), dict) else {}
    b = bindash_cfg(cfg)
    sketchsize64, _, bbits = resolve_bindash_sketch_params(cfg)
    threads = int(b.get("threads", odd.get("threads", 1)))
    kmerlen = int(b.get("kmerlen", odd.get("kmerlen", 64)))
    cmd = [
        resolve_bindash_bin(b.get("bindash_bin")),
        "sketch",
        f"--listfname={input_list}",
        f"--nthreads={threads}",
        f"--kmerlen={kmerlen}",
        f"--sketchsize64={sketchsize64}",
        f"--bbits={bbits}",
        f"--outfname={out_prefix}",
    ]
    if "dens" in b:
        cmd.append(f"--dens={int(b['dens'])}")
    if "minhashtype" in b:
        cmd.append(f"--minhashtype={int(b['minhashtype'])}")
    if "randseed" in b:
        cmd.append(f"--randseed={int(b['randseed'])}")
    if bool(b.get("iscasepreserved", False)):
        cmd.append("--iscasepreserved=true")

    # BinDash preserves strand only when requested. The default reverse-complement
    # handling matches OddSketch canonical=true.
    if "isstrandpreserved" in b:
        if bool(b["isstrandpreserved"]):
            cmd.append("--isstrandpreserved=true")
    elif not bool(odd.get("canonical", True)):
        cmd.append("--isstrandpreserved=true")
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


def run_sketch(input_list: Path, stdout_path: Path, time_path: Path, cfg: dict, out_prefix: Path) -> tuple[float, dict[str, str]]:
    cmd = bindash_cmd(cfg, input_list, out_prefix)
    timed_cmd = ["/usr/bin/time", "-v", "-o", str(time_path), *cmd]
    if not Path("/usr/bin/time").exists():
        timed_cmd = cmd
    print("[run]", " ".join(timed_cmd), flush=True)
    t0 = perf_counter()
    with stdout_path.open("w") as stdout:
        subprocess.run(timed_cmd, stdout=stdout, stderr=subprocess.STDOUT, check=True)
    elapsed = perf_counter() - t0
    return elapsed, parse_time_v(time_path)


def bindash_output_files(out_prefix: Path) -> list[Path]:
    parent = out_prefix.parent
    return sorted(path for path in parent.glob(out_prefix.name + "*") if path.is_file())


def write_bindash_sketch_manifest(path: Path, files: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["path", "bytes"])
        for sketch_file in files:
            writer.writerow([str(sketch_file.resolve()), sketch_file.stat().st_size])


def bindash_version(bindash_bin: str) -> str | None:
    try:
        proc = subprocess.run(
            [bindash_bin, "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except Exception:
        return None
    text = proc.stdout.strip()
    return text or None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--prepare-only", action="store_true")
    args = ap.parse_args()

    cfg_path = resolve_config(args.config)
    cfg = json.loads(cfg_path.read_text())
    data_root = resolve_path(task_root(), cfg.get("paths", {}).get("data_root")) or (task_root() / "outputs")
    run_dir = allocate_run_dir(data_root, args.run_id)
    if args.run_id and run_dir.exists():
        raise SystemExit(f"run already exists: {run_dir}")

    metadata_dir = run_dir / "metadata"
    manifests_dir = run_dir / "manifests"
    input_dir = run_dir / "genome_inputs"
    sketches_dir = run_dir / "bindash_sketches"
    logs_dir = run_dir / "logs"
    for d in (metadata_dir, manifests_dir, input_dir, sketches_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    (metadata_dir / "used_config.json").write_text(json.dumps(cfg, indent=2) + "\n")

    summary_path = copy_assembly_summary(cfg, metadata_dir)
    selected_rows = select_assemblies(summary_path, cfg)
    write_tsv(metadata_dir / "selected_assemblies.tsv", selected_rows, ASSEMBLY_COLUMNS)

    local_list = resolve_path(task_root(), cfg.get("paths", {}).get("local_genome_list"))
    if not local_list:
        raise SystemExit("paths.local_genome_list is required for BinDash sketch. Run the download job first.")
    genomes = read_local_genome_list(local_list)
    limit = refseq_sketch_cfg(cfg).get("limit")
    if limit is not None:
        genomes = genomes[: int(limit)]

    input_list = manifests_dir / "bindash_genome_paths.txt"
    original_inputs = [path for _, path in genomes]
    uses_gzip_inputs = any(is_gzip_path(path) for path in original_inputs)
    linked_inputs = make_input_links(genomes, input_dir)
    write_list(input_list, linked_inputs)

    out_prefix = sketches_dir / "bindash_refseq_sketch"
    sketchsize64, effective_bits, bbits = resolve_bindash_sketch_params(cfg)
    b = bindash_cfg(cfg)
    bindash_bin = resolve_bindash_bin(b.get("bindash_bin"))
    metadata = {
        "run_dir": str(run_dir),
        "refseq_sketch_version_label": refseq_sketch_cfg(cfg).get("version_label"),
        "assembly_summary_source": cfg.get("paths", {}).get("assembly_summary"),
        "assembly_summary_saved": str(summary_path),
        "assembly_summary_comments": assembly_comments(summary_path),
        "selected_assemblies": len(selected_rows),
        "input_genomes": len(linked_inputs),
        "uses_gzip_inputs": uses_gzip_inputs,
        "git_commit": git_commit(),
        "bindash_bin": bindash_bin,
        "bindash_version": bindash_version(bindash_bin),
        "bindash_bin_sha256": sha256_file(Path(bindash_bin)),
        "bindash_sketchsize64": sketchsize64,
        "bindash_effective_sketch_bits": effective_bits,
        "bindash_bbits": bbits,
        "bindash_command": bindash_cmd(cfg, input_list, out_prefix),
    }
    (metadata_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    if args.prepare_only:
        print(f"[prepare] wrote {input_list}")
        print(f"[prepare] run_dir={run_dir}")
        return

    workflow_t0 = perf_counter()
    stdout_path = logs_dir / "bindash_sketch_stdout.txt"
    time_path = logs_dir / "bindash_sketch_time.txt"
    elapsed_sec, time_metrics = run_sketch(input_list, stdout_path, time_path, cfg, out_prefix)
    workflow_elapsed_sec = perf_counter() - workflow_t0

    sketch_files = bindash_output_files(out_prefix)
    write_bindash_sketch_manifest(manifests_dir / "bindash_sketch_files.tsv", sketch_files)
    total_sketch_bytes = sum(path.stat().st_size for path in sketch_files)
    timing = {
        "elapsed_sec": f"{elapsed_sec:.6f}",
        "workflow_elapsed_sec": f"{workflow_elapsed_sec:.6f}",
        "max_rss_kbytes": time_metrics.get("Maximum resident set size (kbytes)", ""),
        "user_sec": time_metrics.get("User time (seconds)", ""),
        "system_sec": time_metrics.get("System time (seconds)", ""),
        "cpu_percent": time_metrics.get("Percent of CPU this job got", ""),
        "exit_status": time_metrics.get("Exit status", ""),
        "input_genomes": str(len(linked_inputs)),
        "sketch_file_count": str(len(sketch_files)),
        "sketchsize64": str(sketchsize64),
        "effective_sketch_bits": str(effective_bits),
        "bbits": str(bbits),
        "total_sketch_bytes": str(total_sketch_bytes),
    }
    write_tsv(run_dir / "results" / "bindash_sketch_metrics.tsv", [timing], list(timing.keys()))
    print(f"[done] bindash_sketch_files={len(sketch_files)} bytes={total_sketch_bytes}")
    print(f"[done] metrics={run_dir / 'results' / 'bindash_sketch_metrics.tsv'}")
    print(f"[done] run_dir={run_dir}")


if __name__ == "__main__":
    main()
