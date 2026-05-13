#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
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


def download_file(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    print(f"[download] {url} -> {out_path}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(out_path)


def copy_or_download_assembly_summary(cfg: dict, metadata_dir: Path) -> Path:
    paths = cfg.get("paths", {})
    refseq = cfg.get("refseq", {})
    source = resolve_path(task_root(), paths.get("assembly_summary"))
    saved = metadata_dir / "assembly_summary_refseq.txt"

    if refseq.get("download_assembly_summary", False):
        url = refseq.get("assembly_summary_url")
        if not url:
            raise SystemExit("refseq.assembly_summary_url is required when download_assembly_summary=true")
        download_file(url, saved)
        return saved

    if not source or not source.exists():
        raise SystemExit(
            "assembly_summary not found. Set paths.assembly_summary, or set "
            "refseq.download_assembly_summary=true."
        )
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


def passes_filters(row: dict[str, str], filters: dict) -> bool:
    taxid = filters.get("taxid")
    if taxid is not None and str(row.get("taxid", "")) != str(taxid):
        return False

    for key in ("assembly_level", "refseq_category"):
        allowed = filters.get(key) or []
        if allowed and row.get(key, "") not in set(map(str, allowed)):
            return False

    return row.get("ftp_path", "") not in ("", "na")


def select_assemblies(summary_path: Path, cfg: dict) -> list[dict[str, str]]:
    refseq = cfg.get("refseq", {})
    filters = refseq.get("filters", {})
    limit = refseq.get("limit")
    rows = [row for row in read_assembly_summary(summary_path) if passes_filters(row, filters)]
    if limit is not None:
        rows = rows[: int(limit)]
    return rows


def genomic_fna_url(ftp_path: str) -> str:
    base = ftp_path.rstrip("/").split("/")[-1]
    url = ftp_path.rstrip("/") + f"/{base}_genomic.fna.gz"
    if url.startswith("ftp://"):
        url = "https://" + url[len("ftp://") :]
    return url


def decompress_gzip(gz_path: Path, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with gzip.open(gz_path, "rb") as inf, tmp.open("wb") as outf:
        shutil.copyfileobj(inf, outf)
    tmp.replace(out_path)


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def materialize_refseq_genomes(rows: list[dict[str, str]], cfg: dict, genomes_dir: Path) -> list[tuple[str, Path]]:
    if not cfg.get("refseq", {}).get("download_genomes", False):
        raise SystemExit(
            "RefSeq assemblies were selected, but refseq.download_genomes=false. "
            "Enable it, or provide paths.local_genome_list."
        )

    downloads_dir = genomes_dir / "downloads"
    fasta_dir = genomes_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    result: list[tuple[str, Path]] = []
    for row in rows:
        accession = row["assembly_accession"]
        url = genomic_fna_url(row["ftp_path"])
        gz_path = downloads_dir / Path(url).name
        fna_path = fasta_dir / f"{safe_name(accession)}.fna"
        if not gz_path.exists():
            download_file(url, gz_path)
        decompress_gzip(gz_path, fna_path)
        result.append((accession, fna_path))
    return result


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
        suffix = "".join(source.suffixes) or ".fna"
        link = input_dir / f"{idx:06d}_{safe_name(name)}{suffix}"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(source)
        linked.append(link)
    return linked


def write_list(path: Path, values: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value.resolve()}\n" for value in values))


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


def oddsketch_cmd(cfg: dict) -> list[str]:
    odd = cfg.get("oddsketch", {})
    return [
        resolve_oddsketch_bin(),
        "sketch",
        f"--kmer={odd.get('kmerlen', 64)}",
        f"--sketch-size={odd.get('sketch_size', 8192)}",
        f"--j0={odd.get('j0', 0.9)}",
        f"--pos-mode={odd.get('pos_mode', 'mix')}",
        f"--canonical={1 if odd.get('canonical', True) else 0}",
        f"--threads={odd.get('threads', 1)}",
    ]


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


def run_sketch(input_list: Path, stdout_path: Path, time_path: Path, cfg: dict) -> tuple[float, dict[str, str]]:
    cmd = oddsketch_cmd(cfg)
    timed_cmd = ["/usr/bin/time", "-v", "-o", str(time_path), *cmd]
    if not Path("/usr/bin/time").exists():
        timed_cmd = cmd
    print("[run]", " ".join(timed_cmd))
    t0 = perf_counter()
    with input_list.open() as stdin, stdout_path.open("w") as stdout:
        subprocess.run(timed_cmd, stdin=stdin, stdout=stdout, check=True)
    elapsed = perf_counter() - t0
    return elapsed, parse_time_v(time_path)


def relocate_sketches(stdout_path: Path, sketches_dir: Path) -> list[Path]:
    sketches_dir.mkdir(parents=True, exist_ok=True)
    relocated: list[Path] = []
    for line in stdout_path.read_text().splitlines():
        raw = line.strip()
        if not raw:
            continue
        src = Path(raw)
        if not src.exists():
            raise SystemExit(f"oddsketch reported missing sketch: {src}")
        dst = sketches_dir / src.name
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), dst)
        relocated.append(dst.resolve())
    return relocated


def metric_float(metrics: dict[str, str], key: str) -> float:
    try:
        return float(metrics.get(key, "") or 0.0)
    except ValueError:
        return 0.0


def metric_int(metrics: dict[str, str], key: str) -> int:
    try:
        return int(metrics.get(key, "") or 0)
    except ValueError:
        return 0


def run_sketch_with_temporary_inputs(
    genomes: list[tuple[str, Path]],
    temp_root: Path,
    manifests_dir: Path,
    logs_dir: Path,
    sketches_dir: Path,
    cfg: dict,
) -> tuple[float, dict[str, str], list[Path], int]:
    odd = cfg.get("oddsketch", {})
    threads = int(odd.get("threads", 1))
    batch_size = max(1, int(odd.get("gzip_batch_size", max(1, threads * 16))))

    total_elapsed = 0.0
    max_rss = 0
    user_sec = 0.0
    system_sec = 0.0
    relocated_all: list[Path] = []
    batch_count = 0

    for batch_start in range(0, len(genomes), batch_size):
        batch_count += 1
        batch = genomes[batch_start : batch_start + batch_size]
        batch_dir = temp_root / f"batch_{batch_count:06d}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_inputs: list[Path] = []
        try:
            for offset, (name, source) in enumerate(batch):
                idx = batch_start + offset + 1
                if is_gzip_path(source):
                    target_name = f"{idx:08d}_{safe_name(ungzip_name(source))}"
                    target = batch_dir / target_name
                    decompress_gzip(source, target)
                else:
                    suffix = "".join(source.suffixes) or ".fna"
                    target = batch_dir / f"{idx:08d}_{safe_name(name)}{suffix}"
                    if target.exists() or target.is_symlink():
                        target.unlink()
                    target.symlink_to(source)
                batch_inputs.append(target)

            batch_input_list = manifests_dir / f"genome_paths_batch_{batch_count:06d}.txt"
            batch_stdout = logs_dir / f"oddsketch_sketch_stdout_batch_{batch_count:06d}.txt"
            batch_time = logs_dir / f"oddsketch_sketch_time_batch_{batch_count:06d}.txt"
            write_list(batch_input_list, batch_inputs)
            elapsed, metrics = run_sketch(batch_input_list, batch_stdout, batch_time, cfg)
            total_elapsed += elapsed
            max_rss = max(max_rss, metric_int(metrics, "Maximum resident set size (kbytes)"))
            user_sec += metric_float(metrics, "User time (seconds)")
            system_sec += metric_float(metrics, "System time (seconds)")
            relocated_all.extend(relocate_sketches(batch_stdout, sketches_dir))
            print(f"[sketch] batch {batch_count}: {len(batch)} genomes")
        finally:
            shutil.rmtree(batch_dir, ignore_errors=True)

    time_metrics = {
        "Maximum resident set size (kbytes)": str(max_rss),
        "User time (seconds)": f"{user_sec:.2f}",
        "System time (seconds)": f"{system_sec:.2f}",
        "Percent of CPU this job got": "",
        "Exit status": "0",
    }
    return total_elapsed, time_metrics, relocated_all, batch_count


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
    args = ap.parse_args()

    cfg_path = resolve_config(args.config)
    cfg = json.loads(cfg_path.read_text())
    data_root = resolve_path(task_root(), cfg.get("paths", {}).get("data_root")) or (task_root() / "outputs")
    run_dir = allocate_run_dir(data_root, args.run_id)
    metadata_dir = run_dir / "metadata"
    manifests_dir = run_dir / "manifests"
    input_dir = run_dir / "genome_inputs"
    genomes_dir = run_dir / "genomes"
    temp_fasta_dir = run_dir / "temporary_fasta"
    sketches_dir = run_dir / "sketches"
    logs_dir = run_dir / "logs"
    for d in (metadata_dir, manifests_dir, input_dir, genomes_dir, temp_fasta_dir, sketches_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    used_config = metadata_dir / "used_config.json"
    used_config.write_text(json.dumps(cfg, indent=2) + "\n")

    summary_path = copy_or_download_assembly_summary(cfg, metadata_dir)
    selected_rows = select_assemblies(summary_path, cfg)
    write_tsv(metadata_dir / "selected_assemblies.tsv", selected_rows, ASSEMBLY_COLUMNS)

    local_list = resolve_path(task_root(), cfg.get("paths", {}).get("local_genome_list"))
    if local_list:
        genomes = read_local_genome_list(local_list)
    else:
        genomes = materialize_refseq_genomes(selected_rows, cfg, genomes_dir)

    input_list = manifests_dir / "genome_paths.txt"
    original_inputs = [path for _, path in genomes]
    uses_gzip_inputs = any(is_gzip_path(path) for path in original_inputs)
    if uses_gzip_inputs:
        write_list(input_list, original_inputs)
        linked_inputs = original_inputs
    else:
        linked_inputs = make_input_links(genomes, input_dir)
        write_list(input_list, linked_inputs)

    metadata = {
        "run_dir": str(run_dir),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "refseq_version_label": cfg.get("refseq", {}).get("version_label"),
        "assembly_summary_source": cfg.get("paths", {}).get("assembly_summary"),
        "assembly_summary_url": cfg.get("refseq", {}).get("assembly_summary_url"),
        "assembly_summary_saved": str(summary_path),
        "assembly_summary_comments": assembly_comments(summary_path),
        "selected_assemblies": len(selected_rows),
        "input_genomes": len(linked_inputs),
        "uses_gzip_inputs": uses_gzip_inputs,
        "gzip_batch_size": cfg.get("oddsketch", {}).get("gzip_batch_size"),
        "git_commit": git_commit(),
        "oddsketch_bin": resolve_oddsketch_bin(),
        "oddsketch_bin_sha256": sha256_file(Path(resolve_oddsketch_bin())),
        "oddsketch_command": oddsketch_cmd(cfg),
    }
    (metadata_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    if args.prepare_only:
        print(f"[prepare] wrote {input_list}")
        print(f"[prepare] run_dir={run_dir}")
        return

    if uses_gzip_inputs:
        elapsed_sec, time_metrics, relocated, batch_count = run_sketch_with_temporary_inputs(
            genomes,
            temp_fasta_dir,
            manifests_dir,
            logs_dir,
            sketches_dir,
            cfg,
        )
    else:
        stdout_path = logs_dir / "oddsketch_sketch_stdout.txt"
        time_path = logs_dir / "oddsketch_sketch_time.txt"
        elapsed_sec, time_metrics = run_sketch(input_list, stdout_path, time_path, cfg)
        relocated = relocate_sketches(stdout_path, sketches_dir)
        batch_count = 1
    sketch_list = manifests_dir / "sketch_paths.txt"
    write_list(sketch_list, relocated)

    total_sketch_bytes = sum(path.stat().st_size for path in relocated)
    timing = {
        "elapsed_sec": f"{elapsed_sec:.6f}",
        "max_rss_kbytes": time_metrics.get("Maximum resident set size (kbytes)", ""),
        "user_sec": time_metrics.get("User time (seconds)", ""),
        "system_sec": time_metrics.get("System time (seconds)", ""),
        "cpu_percent": time_metrics.get("Percent of CPU this job got", ""),
        "exit_status": time_metrics.get("Exit status", ""),
        "input_genomes": str(len(linked_inputs)),
        "sketch_count": str(len(relocated)),
        "batch_count": str(batch_count),
        "temporary_decompression": str(uses_gzip_inputs),
        "total_sketch_bytes": str(total_sketch_bytes),
    }
    write_tsv(run_dir / "results" / "oddsketch_sketch_metrics.tsv", [timing], list(timing.keys()))
    print(f"[done] sketches={len(relocated)} bytes={total_sketch_bytes}")
    print(f"[done] metrics={run_dir / 'results' / 'oddsketch_sketch_metrics.tsv'}")
    print(f"[done] run_dir={run_dir}")


if __name__ == "__main__":
    main()
