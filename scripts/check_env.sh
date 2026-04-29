#!/usr/bin/env bash
set -euo pipefail

echo "[check] uv"
command -v uv >/dev/null
uv --version

echo "[check] bindash"
REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
REPO_BINDASH="${REPO_ROOT}/experiments/tools/bin/bindash"
if [[ -x "${REPO_BINDASH}" ]]; then
  "${REPO_BINDASH}" --version >/dev/null 2>&1 || "${REPO_BINDASH}" --help >/dev/null 2>&1 || true
  echo "bindash: ok (${REPO_BINDASH})"
elif command -v bindash >/dev/null 2>&1; then
  bindash --version >/dev/null 2>&1 || bindash --help >/dev/null 2>&1 || true
  echo "bindash: ok"
else
  echo "bindash: not found (optional unless BinDash baseline is used)"
fi

echo "[check] compiler"
if command -v clang++ >/dev/null 2>&1; then
  clang++ --version | head -n 1
elif command -v g++ >/dev/null 2>&1; then
  g++ --version | head -n 1
else
  echo "No C++ compiler found (need C++17 for src build)"
  exit 1
fi

echo "[check] python deps in uv environment"
uv run python -c "import numpy, pandas, matplotlib; print('python deps: ok')"

echo "[check] done"
