#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import time
import urllib.request
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


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


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_assembly_summary(path: Path) -> tuple[list[str], list[str], list[dict[str, str]]]:
    comments: list[str] = []
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open() as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                if line.startswith("#assembly_accession"):
                    header = line[1:].split("\t")
                continue
            if header is None:
                raise SystemExit(f"assembly_summary header not found in {path}")
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    if header is None:
        raise SystemExit(f"assembly_summary header not found in {path}")
    return comments, header, rows


def genomic_fna_url(ftp_path: str) -> str:
    base = ftp_path.rstrip("/").split("/")[-1]
    url = ftp_path.rstrip("/") + f"/{base}_genomic.fna.gz"
    if url.startswith("ftp://"):
        url = "https://" + url[len("ftp://") :]
    return url


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def download_file(url: str, out_path: Path, retries: int, timeout: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response, tmp.open("wb") as outf:
                shutil.copyfileobj(response, outf)
            tmp.replace(out_path)
            return
        except Exception as exc:
            last_error = exc
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                time.sleep(min(60, 2**attempt))
    raise RuntimeError(f"download failed after {retries} attempts: {url}: {last_error}")


def decompress_gzip(gz_path: Path, fna_path: Path) -> None:
    if fna_path.exists() and fna_path.stat().st_size > 0:
        return
    fna_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = fna_path.with_suffix(fna_path.suffix + ".tmp")
    with gzip.open(gz_path, "rb") as inf, tmp.open("wb") as outf:
        shutil.copyfileobj(inf, outf)
    tmp.replace(fna_path)


def genome_paths(row: dict[str, str], gzip_dir: Path, fasta_dir: Path) -> tuple[str, Path, Path]:
    accession = row.get("assembly_accession", "")
    url = genomic_fna_url(row.get("ftp_path", ""))
    gz_name = Path(url).name
    if not gz_name:
        gz_name = f"{safe_name(accession)}_genomic.fna.gz"
    gz_path = gzip_dir / gz_name
    fna_name = gz_name[:-3] if gz_name.endswith(".gz") else f"{gz_name}.fna"
    fna_path = fasta_dir / fna_name
    return url, gz_path, fna_path


def materialize_one(
    row: dict[str, str],
    gzip_dir: Path,
    fasta_dir: Path,
    retries: int,
    timeout: int,
    decompress: bool,
    dry_run: bool,
) -> dict[str, str]:
    started = time.perf_counter()
    accession = row.get("assembly_accession", "")
    ftp_path = row.get("ftp_path", "")
    organism = row.get("organism_name", "")
    url, gz_path, fna_path = genome_paths(row, gzip_dir, fasta_dir)
    status = "ok"
    error = ""
    downloaded = "0"
    decompressed = "0"

    try:
        if dry_run:
            status = "dry_run"
        else:
            if not gz_path.exists() or gz_path.stat().st_size == 0:
                download_file(url, gz_path, retries=retries, timeout=timeout)
                downloaded = "1"
            if decompress:
                before = fna_path.exists() and fna_path.stat().st_size > 0
                decompress_gzip(gz_path, fna_path)
                decompressed = "0" if before else "1"
    except Exception as exc:
        status = "error"
        error = str(exc)

    elapsed = time.perf_counter() - started
    gz_bytes = gz_path.stat().st_size if gz_path.exists() else 0
    fna_bytes = fna_path.stat().st_size if fna_path.exists() else 0
    return {
        "assembly_accession": accession,
        "organism_name": organism,
        "ftp_path": ftp_path,
        "genomic_fna_url": url,
        "gzip_path": str(gz_path),
        "fasta_path": str(fna_path) if decompress else "",
        "status": status,
        "downloaded": downloaded,
        "decompressed": decompressed,
        "gzip_bytes": str(gz_bytes),
        "fasta_bytes": str(fna_bytes),
        "elapsed_sec": f"{elapsed:.3f}",
        "error": error,
    }


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def select_rows(rows: list[dict[str, str]], limit: int | None) -> list[dict[str, str]]:
    selected = [row for row in rows if row.get("ftp_path", "") not in ("", "na")]
    if limit is not None:
        selected = selected[:limit]
    return selected


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--assembly-summary", default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg_path = resolve_config(args.config)
    cfg = json.loads(cfg_path.read_text())
    download_cfg = cfg.get("download", {})

    assembly_summary = resolve_path(
        task_root(),
        args.assembly_summary
        or download_cfg.get("assembly_summary")
        or "data/refseq_bacteria/assembly_summary.txt",
    )
    if assembly_summary is None or not assembly_summary.exists():
        raise SystemExit(f"assembly_summary not found: {assembly_summary}")

    outdir = resolve_path(task_root(), args.outdir or download_cfg.get("outdir") or "data/assembly")
    if outdir is None:
        raise SystemExit("download output directory could not be resolved")

    threads = max(1, int(args.threads or download_cfg.get("threads", 4)))
    retries = max(1, int(download_cfg.get("retries", 3)))
    timeout = max(1, int(download_cfg.get("timeout_sec", 120)))
    decompress = bool(download_cfg.get("decompress", True))
    limit = args.limit
    if limit is None and download_cfg.get("limit") is not None:
        limit = int(download_cfg["limit"])

    metadata_dir = outdir / "metadata"
    manifests_dir = outdir / "manifests"
    gzip_dir = outdir / "gzip"
    fasta_dir = outdir / "fasta"
    for path in (metadata_dir, manifests_dir, gzip_dir, fasta_dir):
        path.mkdir(parents=True, exist_ok=True)

    started_at = utc_now()
    saved_summary = metadata_dir / "assembly_summary.txt"
    shutil.copy2(assembly_summary, saved_summary)
    comments, header, rows = read_assembly_summary(saved_summary)
    selected = select_rows(rows, limit)

    write_json(
        metadata_dir / "download_metadata.json",
        {
            "status": "running",
            "started_at_utc": started_at,
            "version_label": download_cfg.get("version_label", "RefSeq bacteria assembly_summary.txt"),
            "assembly_summary_source": str(assembly_summary),
            "assembly_summary_saved": str(saved_summary),
            "assembly_summary_sha256": sha256_file(saved_summary),
            "assembly_summary_comments": comments,
            "assembly_summary_columns": header,
            "total_rows": len(rows),
            "selected_rows": len(selected),
            "outdir": str(outdir),
            "threads": threads,
            "decompress": decompress,
            "dry_run": args.dry_run,
        },
    )

    manifest_path = manifests_dir / "assembly_download_manifest.tsv"
    fasta_list_path = manifests_dir / "fasta_paths.txt"
    failed_path = manifests_dir / "failed_assemblies.tsv"
    fieldnames = [
        "assembly_accession",
        "organism_name",
        "ftp_path",
        "genomic_fna_url",
        "gzip_path",
        "fasta_path",
        "status",
        "downloaded",
        "decompressed",
        "gzip_bytes",
        "fasta_bytes",
        "elapsed_sec",
        "error",
    ]

    print(f"[download] assembly_summary={saved_summary}")
    print(f"[download] selected assemblies={len(selected)}")
    print(f"[download] outdir={outdir}")
    print(f"[download] threads={threads} decompress={decompress} dry_run={args.dry_run}")

    ok = 0
    failed = 0
    total_gzip_bytes = 0
    total_fasta_bytes = 0
    selected_iter = iter(selected)

    with manifest_path.open("w", newline="") as mf, failed_path.open("w", newline="") as ff:
        manifest_writer = csv.DictWriter(mf, fieldnames=fieldnames, delimiter="\t")
        failed_writer = csv.DictWriter(ff, fieldnames=fieldnames, delimiter="\t")
        manifest_writer.writeheader()
        failed_writer.writeheader()

        with ThreadPoolExecutor(max_workers=threads) as executor:
            pending = set()
            for _ in range(threads * 2):
                try:
                    row = next(selected_iter)
                except StopIteration:
                    break
                pending.add(executor.submit(materialize_one, row, gzip_dir, fasta_dir, retries, timeout, decompress, args.dry_run))

            done_count = 0
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    result = future.result()
                    manifest_writer.writerow(result)
                    mf.flush()
                    done_count += 1
                    if result["status"] in ("ok", "dry_run"):
                        ok += 1
                        total_gzip_bytes += int(result["gzip_bytes"])
                        total_fasta_bytes += int(result["fasta_bytes"])
                    else:
                        failed += 1
                        failed_writer.writerow(result)
                        ff.flush()

                    if done_count % 100 == 0 or done_count == len(selected):
                        print(f"[download] progress {done_count}/{len(selected)} ok={ok} failed={failed}")

                    try:
                        row = next(selected_iter)
                    except StopIteration:
                        continue
                    pending.add(executor.submit(materialize_one, row, gzip_dir, fasta_dir, retries, timeout, decompress, args.dry_run))

    if decompress and not args.dry_run:
        fasta_paths = sorted(path for path in fasta_dir.glob("*.fna") if path.stat().st_size > 0)
        fasta_list_path.write_text("".join(f"{path.resolve()}\n" for path in fasta_paths))
    else:
        fasta_list_path.write_text("")

    finished_at = utc_now()
    write_json(
        metadata_dir / "download_metadata.json",
        {
            "status": "completed" if failed == 0 else "completed_with_errors",
            "started_at_utc": started_at,
            "finished_at_utc": finished_at,
            "version_label": download_cfg.get("version_label", "RefSeq bacteria assembly_summary.txt"),
            "assembly_summary_source": str(assembly_summary),
            "assembly_summary_saved": str(saved_summary),
            "assembly_summary_sha256": sha256_file(saved_summary),
            "assembly_summary_comments": comments,
            "assembly_summary_columns": header,
            "total_rows": len(rows),
            "selected_rows": len(selected),
            "ok": ok,
            "failed": failed,
            "total_gzip_bytes": total_gzip_bytes,
            "total_fasta_bytes": total_fasta_bytes,
            "outdir": str(outdir),
            "manifest": str(manifest_path),
            "failed_manifest": str(failed_path),
            "fasta_list": str(fasta_list_path),
            "threads": threads,
            "retries": retries,
            "timeout_sec": timeout,
            "decompress": decompress,
            "dry_run": args.dry_run,
        },
    )
    print(f"[download] done ok={ok} failed={failed}")
    print(f"[download] manifest={manifest_path}")
    print(f"[download] fasta_list={fasta_list_path}")
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
