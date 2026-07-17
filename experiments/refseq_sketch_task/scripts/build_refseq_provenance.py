#!/usr/bin/env python3
"""Build the public provenance files for the paper RefSeq dataset."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

from assembly_summary_io import content_sha256, content_size, open_text


MANIFEST_COLUMNS = (
    "assembly_accession",
    "ftp_path",
    "genomic_fna_url",
    "local_filename",
    "file_size",
)


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def repo_root() -> Path:
    return task_root().parents[1]


def resolve_path(raw: str | Path) -> Path:
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (task_root() / path).resolve()


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root()))
    except ValueError:
        return str(path.resolve())


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_assembly_summary(path: Path) -> tuple[dict[str, str], int, int]:
    header: list[str] | None = None
    ftp_paths: dict[str, str] = {}
    total_rows = 0
    eligible_rows = 0

    with open_text(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("#assembly_accession"):
                header = line[1:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise SystemExit(f"assembly summary header not found: {path}")

            total_rows += 1
            values = line.split("\t")
            row = dict(zip(header, values))
            accession = row.get("assembly_accession", "")
            ftp_path = row.get("ftp_path", "")
            if accession and ftp_path not in ("", "na"):
                ftp_paths[accession] = ftp_path
                eligible_rows += 1

    return ftp_paths, total_rows, eligible_rows


def expected_genomic_url(ftp_path: str) -> str:
    directory = ftp_path.rstrip("/")
    basename = directory.split("/")[-1]
    url = f"{directory}/{basename}_genomic.fna.gz"
    if url.startswith("ftp://"):
        url = "https://" + url[len("ftp://") :]
    return url


def read_failed_assemblies(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    fields = ("assembly_accession", "ftp_path", "genomic_fna_url", "error")
    with path.open(newline="") as f:
        return [
            {field: row.get(field, "") for field in fields}
            for row in csv.DictReader(f, delimiter="\t")
        ]


def write_manifest(
    integrity_manifest: Path,
    ftp_paths: dict[str, str],
    output_path: Path,
    verify_files: bool,
) -> tuple[int, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    rows_written = 0
    total_file_size = 0

    with (
        integrity_manifest.open(newline="") as input_file,
        output_path.open("wb") as raw_output,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as gzip_output,
        io.TextIOWrapper(gzip_output, encoding="utf-8", newline="") as text_output,
    ):
        reader = csv.DictReader(input_file, delimiter="\t")
        writer = csv.DictWriter(
            text_output,
            fieldnames=MANIFEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            if row.get("status") != "ok":
                raise SystemExit(
                    f"non-ok integrity row for {row.get('assembly_accession', '')}: "
                    f"{row.get('status', '')}"
                )

            accession = row.get("assembly_accession", "")
            if not accession or accession in seen:
                raise SystemExit(f"missing or duplicate assembly accession: {accession}")
            seen.add(accession)

            ftp_path = ftp_paths.get(accession)
            if not ftp_path:
                raise SystemExit(f"accession not found in assembly summary: {accession}")

            genomic_url = row.get("genomic_fna_url", "")
            expected_url = expected_genomic_url(ftp_path)
            if genomic_url != expected_url:
                raise SystemExit(
                    f"genomic URL mismatch for {accession}: "
                    f"manifest={genomic_url}, expected={expected_url}"
                )

            local_path = Path(row.get("gzip_path", ""))
            recorded_size = int(row.get("after_bytes", "0"))
            if verify_files:
                if not local_path.is_file():
                    raise SystemExit(f"local genome file not found: {local_path}")
                actual_size = local_path.stat().st_size
                if actual_size != recorded_size:
                    raise SystemExit(
                        f"file size mismatch for {accession}: "
                        f"integrity manifest={recorded_size}, actual={actual_size}"
                    )

            writer.writerow(
                {
                    "assembly_accession": accession,
                    "ftp_path": ftp_path,
                    "genomic_fna_url": genomic_url,
                    "local_filename": local_path.name,
                    "file_size": recorded_size,
                }
            )
            rows_written += 1
            total_file_size += recorded_size

    return rows_written, total_file_size


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument(
        "--integrity-manifest",
        default="data/assembly/manifests/gzip_integrity_manifest.tsv",
    )
    ap.add_argument(
        "--integrity-metadata",
        default="data/assembly/manifests/gzip_integrity_metadata.json",
    )
    ap.add_argument(
        "--download-metadata",
        default="data/assembly/metadata/download_metadata.json",
    )
    ap.add_argument(
        "--failed-manifest",
        default="data/assembly/manifests/failed_assemblies.tsv",
    )
    ap.add_argument(
        "--manifest-out",
        default="provenance/refseq_bacteria_genomes.tsv.gz",
    )
    ap.add_argument(
        "--metadata-out",
        default="provenance/refseq_bacteria_dataset.json",
    )
    ap.add_argument("--no-verify-files", action="store_true")
    args = ap.parse_args()

    config_path = resolve_path(args.config)
    cfg = json.loads(config_path.read_text())
    download_cfg = cfg.get("download", {})

    assembly_summary = resolve_path(download_cfg["assembly_summary"])
    integrity_manifest = resolve_path(args.integrity_manifest)
    integrity_metadata_path = resolve_path(args.integrity_metadata)
    download_metadata_path = resolve_path(args.download_metadata)
    failed_manifest = resolve_path(args.failed_manifest)
    manifest_out = resolve_path(args.manifest_out)
    metadata_out = resolve_path(args.metadata_out)

    summary_sha256 = content_sha256(assembly_summary)
    expected_summary_sha256 = download_cfg.get("assembly_summary_sha256")
    if expected_summary_sha256 and summary_sha256 != expected_summary_sha256:
        raise SystemExit(
            f"assembly summary SHA256 mismatch: "
            f"expected={expected_summary_sha256}, actual={summary_sha256}"
        )

    ftp_paths, total_rows, eligible_rows = read_assembly_summary(assembly_summary)
    rows_written, total_file_size = write_manifest(
        integrity_manifest,
        ftp_paths,
        manifest_out,
        verify_files=not args.no_verify_files,
    )

    download_metadata = json.loads(download_metadata_path.read_text())
    integrity_metadata = json.loads(integrity_metadata_path.read_text())
    failed_assemblies = read_failed_assemblies(failed_manifest)
    excluded_accessions = {
        str(value) for value in download_cfg.get("excluded_accessions", [])
    }
    expected_genome_count = download_cfg.get("expected_genome_count")
    if expected_genome_count is not None and rows_written != int(expected_genome_count):
        raise SystemExit(
            f"provenance manifest row count mismatch: "
            f"expected={expected_genome_count}, actual={rows_written}"
        )
    failed_accessions = {
        row["assembly_accession"] for row in failed_assemblies
    }
    if failed_accessions != excluded_accessions:
        raise SystemExit(
            f"excluded accessions do not match the recorded download failures: "
            f"excluded={sorted(excluded_accessions)}, "
            f"failed={sorted(failed_accessions)}"
        )

    metadata = {
        "schema_version": 1,
        "dataset": "NCBI RefSeq bacteria genomes used in the OddSketch-Genome paper experiments",
        "assembly_summary": {
            "source_url": download_cfg.get("assembly_summary_source_url"),
            "acquired_on": download_cfg.get("assembly_summary_acquired_on"),
            "source_last_modified_at": download_cfg.get(
                "assembly_summary_source_last_modified_at"
            ),
            "local_path": download_cfg.get("assembly_summary"),
            "compression": "gzip" if assembly_summary.suffix == ".gz" else "none",
            "file_size": content_size(assembly_summary),
            "sha256": summary_sha256,
            "compressed_file_size": (
                assembly_summary.stat().st_size
                if assembly_summary.suffix == ".gz"
                else None
            ),
            "compressed_sha256": (
                sha256_file(assembly_summary)
                if assembly_summary.suffix == ".gz"
                else None
            ),
            "total_rows": total_rows,
            "eligible_rows": eligible_rows,
            "excluded_accessions": sorted(excluded_accessions),
            "expected_genome_count": expected_genome_count,
        },
        "genome_download": {
            "started_at_utc": download_metadata.get("started_at_utc"),
            "finished_at_utc": download_metadata.get("finished_at_utc"),
            "successful_files": rows_written,
            "failed_assemblies": failed_assemblies,
        },
        "gzip_integrity_check": {
            "started_at_utc": integrity_metadata.get("started_at_utc"),
            "finished_at_utc": integrity_metadata.get("finished_at_utc"),
            "checked": integrity_metadata.get("checked"),
            "remaining_invalid": integrity_metadata.get("remaining_invalid"),
        },
        "genome_manifest": {
            "path": repo_relative(manifest_out),
            "compression": "gzip",
            "columns": list(MANIFEST_COLUMNS),
            "rows": rows_written,
            "total_file_size": total_file_size,
            "compressed_file_size": manifest_out.stat().st_size,
            "sha256": sha256_file(manifest_out),
        },
    }

    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    metadata_out.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote: {manifest_out}")
    print(f"wrote: {metadata_out}")


if __name__ == "__main__":
    main()
