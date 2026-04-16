#!/usr/bin/env python3

import argparse
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    task_root = resolve_task_root()
    cfg_path = resolve_config_path(args.config)
    cfg = json.loads(cfg_path.read_text())
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))
    outdir.mkdir(parents=True, exist_ok=True)

    db_list = outdir / "db_genomes.list"
    q_list = outdir / "queries.list"
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    combined = outdir / "all_genomes.list"
    with combined.open("w") as cf:
        cf.write(db_list.read_text())
        cf.write(q_list.read_text())

    k = int(cfg.get("true_jaccard", {}).get("kmerlen", 64))
    db_paths = [line.strip() for line in db_list.read_text().splitlines() if line.strip()]
    q_paths = [line.strip() for line in q_list.read_text().splitlines() if line.strip()]
    for fasta in set(db_paths + q_paths):
        for suffix in (f".k{k}.bin", f".k{k}.idx"):
            idx_path = Path(fasta + suffix)
            if idx_path.exists():
                try:
                    idx_path.unlink()
                except Exception:
                    pass

    cpp = task_root.parents[1] / "experiments" / "tools" / "bin" / "true_index_pairs"
    if not cpp.exists():
        raise SystemExit(f"binary not found: {cpp}")

    cmd1 = [str(cpp), "preprocess", "--list", str(combined), "--k", str(k)]
    print("[run]", " ".join(cmd1))
    subprocess.run(cmd1, check=True)

    true_pairs = outdir / "true_pairs.tsv"
    true_nn = outdir / "true_nn.tsv"
    cmd2 = [
        str(cpp),
        "pairs",
        "--qlist", str(q_list),
        "--dblist", str(db_list),
        "--out-pairs", str(true_pairs),
        "--out-nn", str(true_nn),
        "--k", str(k),
    ]
    print("[run]", " ".join(cmd2))
    subprocess.run(cmd2, check=True)
    print(f"[true] wrote {true_pairs}")
    print(f"[true] wrote {true_nn}")


if __name__ == "__main__":
    main()
