#!/usr/bin/env python3

import argparse
import csv
import math
import json
import subprocess
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(base: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    return path if path.is_absolute() else (base / path).resolve()


def resolve_config_path(config_arg: str) -> Path:
    task_root = resolve_task_root()
    candidates = [
        Path(config_arg),
        task_root / config_arg,
        Path(__file__).resolve().parent / config_arg,
        task_root / "config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (task_root / "config.json").resolve()


def load_cfg(cpath: Path) -> dict:
    try:
        return json.loads(cpath.read_text())
    except Exception:
        return {}


def resolve_output_root(task_root: Path, cfg: dict) -> Path:
    paths = cfg.get("paths", {}) if isinstance(cfg, dict) else {}
    if isinstance(paths, dict) and paths.get("outdir"):
        return resolve_path(task_root, paths["outdir"])
    legacy = cfg.get("make_genomes", {}).get("outdir") if isinstance(cfg, dict) else None
    if legacy:
        legacy_path = resolve_path(task_root, legacy)
        return legacy_path.parent if legacy_path.name == "genomes" else legacy_path
    return (task_root / "outputs" / "default").resolve()


def read_pair_info(pair_info_path: Path) -> list[dict]:
    pairs = []
    with pair_info_path.open() as f:
        next(f, None)
        for line in f:
            pid, f1, f2, mut, glen = line.strip().split("\t")
            pairs.append({
                "pair_id": int(pid),
                "file1": f1,
                "file2": f2,
                "mutation_count": int(mut),
                "genome_length": int(glen),
            })
    return pairs


def run_cmd(cmd: str, capture: bool = True) -> str:
    print(f"[run] {cmd}")
    process = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
        text=True,
    )
    if process.stderr and process.stderr.strip():
        print(process.stderr.strip())
    return process.stdout if capture else ""


def resolve_bindash_sketch_params(cfg: dict) -> tuple[int, int, int]:
    bbits = int(cfg.get("bbits", 16))
    if bbits <= 0:
        raise SystemExit("bindash.bbits must be positive")
    if "sketch_size" in cfg:
        target_bits = int(cfg["sketch_size"])
    else:
        target_bits = 64 * int(cfg.get("sketchsize64", 32)) * bbits
    if target_bits <= 0:
        raise SystemExit("bindash.sketch_size must be positive")
    sketchsize64 = max(1, math.ceil(target_bits / (64 * bbits)))
    effective_bits = 64 * sketchsize64 * bbits
    return sketchsize64, effective_bits, bbits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json", help="Path to task config JSON")
    args = ap.parse_args()

    task_root = resolve_task_root()
    config_path = resolve_config_path(args.config)
    full_cfg = load_cfg(config_path)
    cfg = full_cfg.get("bindash", {}) if isinstance(full_cfg.get("bindash"), dict) else {}
    output_root = resolve_output_root(task_root, full_cfg)
    pair_info = output_root / "pair_info.txt"
    results_dir = output_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    if not pair_info.exists():
        raise SystemExit(f"pair_info が見つかりません: {pair_info}")

    bindash_bin = cfg.get("bindash_bin", "bindash")
    threads = int(cfg.get("threads", 8))
    mode = cfg.get("mode", "sketch_dist")
    kmerlen = int(cfg.get("kmerlen", 64))
    sketchsize64, _, bbits = resolve_bindash_sketch_params(cfg)
    pair_cmd = cfg.get("pair_cmd", "{bin} exact --nthreads={threads} {f1} {f2}")

    pairs = read_pair_info(pair_info)
    results = []

    if mode == "pairwise":
        for pair in pairs:
            cmd = pair_cmd.format(bin=bindash_bin, f1=pair["file1"], f2=pair["file2"], threads=threads)
            out = run_cmd(cmd, capture=True)
            line = next((ln for ln in reversed(out.splitlines()) if ln.strip()), None)
            if not line:
                continue
            parts = line.strip().split("\t")
            if len(parts) == 3:
                jaccard = float(parts[2])
            elif len(parts) >= 5:
                jac = parts[4]
                jaccard = float(jac.split("/")[0]) / float(jac.split("/")[1]) if "/" in jac else float(jac)
            else:
                continue
            results.append({
                **pair,
                "jaccard_bindash": jaccard,
            })
    elif mode == "sketch_dist":
        listfile = results_dir / "bindash_sketch_list.txt"
        with listfile.open("w") as f:
            for pair in pairs:
                f.write(pair["file1"] + "\n")
                f.write(pair["file2"] + "\n")

        sketch_path = results_dir / "bindash_sketch"
        sketch_cmd = (
            f"{bindash_bin} sketch --listfname={listfile} "
            f"--nthreads={threads} --kmerlen={kmerlen} "
            f"--sketchsize64={sketchsize64} --bbits={bbits} "
            f"--outfname={sketch_path}"
        )
        run_cmd(sketch_cmd, capture=False)

        dist_cmd = f"{bindash_bin} dist --nthreads={threads} --outfname=- {sketch_path}"
        dist_out = run_cmd(dist_cmd, capture=True)

        pair_index = {(pair["file1"], pair["file2"]): pair for pair in pairs}
        pair_index_rev = {(pair["file2"], pair["file1"]): pair for pair in pairs}
        for line in dist_out.splitlines():
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            q, t, _, _, jac = parts[:5]
            pair = pair_index.get((q, t)) or pair_index_rev.get((t, q))
            if not pair:
                continue
            jaccard = float(jac.split("/")[0]) / float(jac.split("/")[1]) if "/" in jac else float(jac)
            results.append({
                **pair,
                "jaccard_bindash": jaccard,
            })
    else:
        raise SystemExit("Unsupported mode. Use 'pairwise' or 'sketch_dist'.")

    out_txt = results_dir / "jaccard_bindash_results.txt"
    with out_txt.open("w") as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_bindash\tfile1\tfile2\n")
        for result in results:
            f.write(
                f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                f"{result['jaccard_bindash']:.10f}\t{result['file1']}\t{result['file2']}\n"
            )
    print(f"wrote {out_txt}")

    true_path = results_dir / "jaccard_true_results.txt"
    if true_path.exists():
        true = {}
        with true_path.open() as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if row:
                    true[int(row[0])] = {
                        "mutation_count": int(row[1]),
                        "jaccard_true": float(row[4]),
                    }
        out_csv = results_dir / "comparison_results_bindash.csv"
        with out_csv.open("w") as f:
            writer = csv.writer(f)
            writer.writerow(["pair_id", "mutation_count", "jaccard_true", "jaccard_bindash"])
            for result in results:
                truth = true.get(result["pair_id"])
                if truth:
                    writer.writerow([
                        result["pair_id"],
                        truth["mutation_count"],
                        truth["jaccard_true"],
                        result["jaccard_bindash"],
                    ])
        print(f"wrote {out_csv}")
    else:
        print(f"note: true results not found ({true_path}); skipped CSV merge.")


if __name__ == "__main__":
    main()
