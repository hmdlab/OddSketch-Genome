#!/usr/bin/env python3

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_repo_root() -> Path:
    return resolve_task_root().parents[1]


def display_path(raw: str | Path) -> str:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.resolve().relative_to(resolve_repo_root()))
    except Exception:
        return str(path)


def display_cmd(cmd: list[str]) -> str:
    return " ".join(display_path(part) for part in cmd)


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

    manifests_dir = outdir / "data" / "manifests"
    db_list = manifests_dir / "db_genome_paths.txt"
    q_list = manifests_dir / "query_genome_paths.txt"
    if not db_list.exists() or not q_list.exists():
        raise SystemExit(f"missing inputs: {db_list} / {q_list}")

    k = int(cfg.get("true_jaccard", {}).get("kmerlen", 64))
    db_paths = [line.strip() for line in db_list.read_text().splitlines() if line.strip()]
    q_paths = [line.strip() for line in q_list.read_text().splitlines() if line.strip()]

    print("[true] start exact Jaccard preprocessing and pair search")
    print(f"[true] config={display_path(cfg_path)}")
    print(f"[true] output_root={display_path(outdir)}")
    print(
        "[true] inputs="
        f"queries={len(q_paths)}, db={len(db_paths)}, total_pairs={len(q_paths) * len(db_paths)}, k={k}"
    )
    print(f"[true] db_list={display_path(db_list)}")
    print(f"[true] query_list={display_path(q_list)}")

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

    truth_dir = outdir / "results" / "truth"
    truth_dir.mkdir(parents=True, exist_ok=True)
    exact_pairs = truth_dir / "exact_query_db_jaccard.tsv"
    exact_top1 = truth_dir / "exact_top1_neighbors.tsv"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".list") as tmp:
        tmp.write(db_list.read_text())
        tmp.write(q_list.read_text())
        combined_path = Path(tmp.name)

    try:
        print("[true] stage 1/2 preprocess indices")
        cmd1 = [str(cpp), "preprocess", "--list", str(combined_path), "--k", str(k)]
        print("[run]", display_cmd(cmd1))
        subprocess.run(cmd1, check=True)

        print("[true] stage 2/2 compute exact query-db pairs")
        cmd2 = [
            str(cpp),
            "pairs",
            "--qlist", str(q_list),
            "--dblist", str(db_list),
            "--out-pairs", str(exact_pairs),
            "--out-nn", str(exact_top1),
            "--k", str(k),
        ]
        print("[run]", display_cmd(cmd2))
        subprocess.run(cmd2, check=True)
    finally:
        try:
            combined_path.unlink()
        except Exception:
            pass
        for fasta in set(db_paths + q_paths):
            for suffix in (f".k{k}.bin", f".k{k}.idx"):
                idx_path = Path(fasta + suffix)
                if idx_path.exists():
                    try:
                        idx_path.unlink()
                    except Exception:
                        pass

    print(f"[true] wrote exact query-db Jaccard -> {display_path(exact_pairs)}")
    print(f"[true] wrote exact top1 neighbors   -> {display_path(exact_top1)}")


if __name__ == "__main__":
    main()
