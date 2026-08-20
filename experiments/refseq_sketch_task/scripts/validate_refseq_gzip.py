#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import csv
import gzip
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
    candidates = [
        Path(path_arg).expanduser(),
        task_root() / path_arg,
        task_root() / "config.json",
    ]
    for path in candidates:
        if path.exists():
            return path.resolve()
    raise SystemExit(f"config not found: {path_arg}")


def resolve_path(raw: str | Path) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (task_root() / path).resolve()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_gzip_paths(path: Path) -> list[Path]:
    result: list[Path] = []
    with path.open() as f:
        for raw in f:
            value = raw.strip()
            if not value or value.startswith("#"):
                continue
            genome = Path(value).expanduser()
            if not genome.is_absolute():
                genome = (path.parent / genome).resolve()
            result.append(genome)
    return result


def read_download_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as f:
        rows = csv.DictReader(f, delimiter="\t")
        return {str(Path(row["gzip_path"]).resolve()): row for row in rows if row.get("gzip_path")}


def gzip_ok(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as inf:
            while inf.read(1024 * 1024):
                pass
        return True, ""
    except Exception as exc:
        return False, str(exc)


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


def validate_one(
    path: Path,
    manifest_rows: dict[str, dict[str, str]],
    repair: bool,
    retries: int,
    timeout: int,
) -> dict[str, str]:
    started = time.perf_counter()
    row = manifest_rows.get(str(path.resolve()), {})
    before_bytes = path.stat().st_size if path.exists() else 0
    ok_before, error_before = gzip_ok(path) if path.exists() else (False, "missing")
    repaired = "0"
    redownloaded = "0"
    status = "ok" if ok_before else "invalid"
    error_after = ""

    if not ok_before and repair:
        url = row.get("genomic_fna_url", "")
        if not url:
            status = "repair_error"
            error_after = "genomic_fna_url not found in download manifest"
        else:
            try:
                if path.exists():
                    path.unlink()
                download_file(url, path, retries=retries, timeout=timeout)
                redownloaded = "1"
                ok_after, error_after = gzip_ok(path)
                repaired = "1" if ok_after else "0"
                status = "repaired" if ok_after else "repair_error"
            except Exception as exc:
                status = "repair_error"
                error_after = str(exc)

    elapsed = time.perf_counter() - started
    return {
        "assembly_accession": row.get("assembly_accession", ""),
        "gzip_path": str(path),
        "genomic_fna_url": row.get("genomic_fna_url", ""),
        "before_bytes": str(before_bytes),
        "after_bytes": str(path.stat().st_size if path.exists() else 0),
        "status": status,
        "valid_before": "1" if ok_before else "0",
        "repaired": repaired,
        "redownloaded": redownloaded,
        "error_before": error_before,
        "error_after": error_after,
        "elapsed_sec": f"{elapsed:.3f}",
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--gzip-list", default=None)
    ap.add_argument("--download-manifest", default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--retries", type=int, default=None)
    ap.add_argument("--timeout-sec", type=int, default=None)
    ap.add_argument("--repair", action="store_true")
    args = ap.parse_args()

    config_path = resolve_config(args.config)
    config = json.loads(config_path.read_text())
    download_config = config.get("download", {})
    download_outdir = resolve_path(
        download_config.get("outdir") or "outputs/assembly"
    )
    manifests_dir = download_outdir / "manifests"

    gzip_list = (
        resolve_path(args.gzip_list)
        if args.gzip_list
        else manifests_dir / "gzip_paths.txt"
    )
    download_manifest = (
        resolve_path(args.download_manifest)
        if args.download_manifest
        else manifests_dir / "assembly_download_manifest.tsv"
    )
    outdir = resolve_path(args.outdir) if args.outdir else manifests_dir
    if not gzip_list.is_file():
        raise SystemExit(f"gzip list not found: {gzip_list}")
    if not download_manifest.is_file():
        raise SystemExit(f"download manifest not found: {download_manifest}")

    paths = read_gzip_paths(gzip_list)
    manifest_rows = read_download_manifest(download_manifest)
    threads = max(1, int(args.threads or download_config.get("threads", 4)))
    retries = max(1, int(args.retries or download_config.get("retries", 3)))
    timeout = max(1, int(args.timeout_sec or download_config.get("timeout_sec", 120)))
    fieldnames = [
        "assembly_accession",
        "gzip_path",
        "genomic_fna_url",
        "before_bytes",
        "after_bytes",
        "status",
        "valid_before",
        "repaired",
        "redownloaded",
        "error_before",
        "error_after",
        "elapsed_sec",
    ]

    started_at = utc_now()
    rows: list[dict[str, str]] = []
    pending_paths = iter(paths)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        pending = set()
        for _ in range(threads * 2):
            try:
                path = next(pending_paths)
            except StopIteration:
                break
            pending.add(
                executor.submit(
                    validate_one,
                    path,
                    manifest_rows,
                    args.repair,
                    retries,
                    timeout,
                )
            )

        done_count = 0
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                rows.append(future.result())
                done_count += 1
                if done_count % 1000 == 0 or done_count == len(paths):
                    invalid = sum(row["status"] not in ("ok", "repaired") for row in rows)
                    repaired = sum(row["status"] == "repaired" for row in rows)
                    print(f"[gzip-check] progress {done_count}/{len(paths)} invalid={invalid} repaired={repaired}")
                try:
                    path = next(pending_paths)
                except StopIteration:
                    continue
                pending.add(
                    executor.submit(
                        validate_one,
                        path,
                        manifest_rows,
                        args.repair,
                        retries,
                        timeout,
                    )
                )

    rows.sort(key=lambda row: row["gzip_path"])
    report_path = outdir / "gzip_integrity_manifest.tsv"
    invalid_path = outdir / "invalid_gzip_files.tsv"
    write_tsv(report_path, rows, fieldnames)
    write_tsv(
        invalid_path,
        [row for row in rows if row["status"] not in ("ok", "repaired")],
        fieldnames,
    )
    metadata = {
        "started_at_utc": started_at,
        "finished_at_utc": utc_now(),
        "config": str(config_path),
        "download_outdir": str(download_outdir),
        "gzip_list": str(gzip_list),
        "download_manifest": str(download_manifest),
        "checked": len(rows),
        "ok_before": sum(row["valid_before"] == "1" for row in rows),
        "invalid_before": sum(row["valid_before"] != "1" for row in rows),
        "repaired": sum(row["status"] == "repaired" for row in rows),
        "remaining_invalid": sum(row["status"] not in ("ok", "repaired") for row in rows),
        "repair": args.repair,
        "threads": threads,
        "retries": retries,
        "timeout_sec": timeout,
        "report": str(report_path),
        "invalid_report": str(invalid_path),
    }
    (outdir / "gzip_integrity_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"[gzip-check] report={report_path}")
    print(f"[gzip-check] invalid_report={invalid_path}")
    print(
        "[gzip-check] done "
        f"checked={metadata['checked']} invalid_before={metadata['invalid_before']} "
        f"repaired={metadata['repaired']} remaining_invalid={metadata['remaining_invalid']}"
    )
    if metadata["remaining_invalid"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
