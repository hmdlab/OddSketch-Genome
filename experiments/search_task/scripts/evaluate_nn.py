#!/usr/bin/env python3

import argparse
import json
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


def load_map(tsv_path: Path, q_col: int, nn_col: int) -> dict:
    mapping = {}
    with tsv_path.open() as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > max(q_col, nn_col):
                mapping[parts[q_col]] = parts[nn_col]
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--skip-bindash", action="store_true", help="Evaluate only OddSketch columns")
    args = ap.parse_args()

    task_root = resolve_task_root()
    cfg = json.loads(resolve_config_path(args.config).read_text())
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))

    odd_tsv = outdir / "oddsketch_nn.tsv"
    label_path = outdir / "true_nn.tsv"
    bds_tsv = outdir / "bindash_nn.tsv"
    need_bindash = not args.skip_bindash
    if not odd_tsv.exists() or not label_path.exists() or (need_bindash and not bds_tsv.exists()):
        raise SystemExit(f"missing inputs under {outdir}")

    true_map = load_map(label_path, 0, 1)
    odd_map = load_map(odd_tsv, 0, 1)
    bds_map = load_map(bds_tsv, 0, 1) if need_bindash else {}

    ok_odd = 0
    ok_bds = 0
    eval_path = outdir / "nn_eval.tsv"
    with eval_path.open("w") as f:
        f.write("query\tnn_true\tnn_oddsketch\tcorrect_oddsketch\tnn_bindash\tcorrect_bindash\n")
        for query, truth in true_map.items():
            odd_pred = odd_map.get(query, "")
            bds_pred = bds_map.get(query, "") if need_bindash else ""
            c_odd = int(odd_pred == truth)
            c_bds = int(bds_pred == truth) if need_bindash else 0
            ok_odd += c_odd
            ok_bds += c_bds
            f.write(f"{query}\t{truth}\t{odd_pred}\t{c_odd}\t{bds_pred}\t{c_bds}\n")

    n = len(true_map)
    print(f"[eval] queries={n}")
    print(f"[eval] oddsketch top1 accuracy = {ok_odd}/{n} ({(ok_odd / n * 100):.2f}%)")
    if need_bindash:
        print(f"[eval] bindash   top1 accuracy = {ok_bds}/{n} ({(ok_bds / n * 100):.2f}%)")
    print(f"[eval] wrote {eval_path}")


if __name__ == "__main__":
    main()
