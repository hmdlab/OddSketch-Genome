#!/bin/bash
#$ -cwd
#$ -V
#$ -N refseq_gzip_check
#$ -q tsmall
#$ -pe OpenMP 4
#$ -l mem_req=8g
#$ -l h_vmem=8g

set -euo pipefail

REPO_ROOT=$(pwd)
if [[ ! -f "${REPO_ROOT}/README.md" || ! -d "${REPO_ROOT}/experiments/refseq_sketch_task" ]]; then
  echo "Run this job from the repository root:" >&2
  echo "  cd /path/to/genome-oddsketch && qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh" >&2
  exit 1
fi

TASK_DIR="${REPO_ROOT}/experiments/refseq_sketch_task"

UV_BIN=${UV_BIN:-"${HOME}/.local/bin/uv"}
if [[ ! -x "${UV_BIN}" ]]; then
  UV_BIN=$(command -v uv || true)
fi
if [[ -z "${UV_BIN}" || ! -x "${UV_BIN}" ]]; then
  echo "uv not found. Set UV_BIN or add uv to PATH." >&2
  exit 1
fi

if [[ -x /usr/bin/python3.11 ]]; then
  PYTHON_BIN_DEFAULT=/usr/bin/python3.11
else
  PYTHON_BIN_DEFAULT=$(command -v python3.11 || true)
fi
PYTHON_BIN=${PYTHON_BIN:-"${PYTHON_BIN_DEFAULT}"}

export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${USER}"}
mkdir -p "${UV_CACHE_DIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

CHECK_THREADS=${CHECK_THREADS:-${NSLOTS:-4}}

echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] uv=${UV_BIN}"
if [[ -n "${PYTHON_BIN}" ]]; then
  echo "[job] python=${PYTHON_BIN}"
fi
echo "[job] uv_cache_dir=${UV_CACHE_DIR}"
echo "[job] check_threads=${CHECK_THREADS}"

"${UV_BIN}" sync
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python "${TASK_DIR}/scripts/validate_refseq_gzip.py" \
  --threads "${CHECK_THREADS}" \
  --repair

echo "[job] end=$(date)"
