#!/usr/bin/env python3
"""Run the small local paired-comparison smoke configuration."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> None:
    task_root = Path(__file__).resolve().parents[2]
    command = [
        sys.executable,
        str(task_root / "scripts" / "runners" / "run_paired.py"),
        "--config",
        str(task_root / "configs" / "smoke.json"),
        "--experiment",
        "paired",
        *sys.argv[1:],
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
