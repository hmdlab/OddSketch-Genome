#!/bin/bash
#$ -cwd
#$ -V
#$ -N refseq_thread_sweep
#$ -q tsmall
#$ -pe OpenMP 16
#$ -l mem_req=32g
#$ -l h_vmem=32g

set -euo pipefail

REPO_ROOT=$(pwd)
if [[ ! -f "${REPO_ROOT}/README.md" || ! -d "${REPO_ROOT}/experiments/refseq_sketch_task" ]]; then
  echo "Run this job from the repository root:" >&2
  echo "  cd /path/to/genome-oddsketch && qsub experiments/refseq_sketch_task/jobs/qsub_refseq_thread_sweep_1024.sh" >&2
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

ODDSKETCH_BIN=${ODDSKETCH_BIN:-"${REPO_ROOT}/src/oddsketch"}
if [[ ! -x "${ODDSKETCH_BIN}" ]]; then
  echo "warning: oddsketch binary not found at ${ODDSKETCH_BIN}" >&2
fi

export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${USER}"}
export ODDSKETCH_BIN
mkdir -p "${UV_CACHE_DIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] uv=${UV_BIN}"
if [[ -n "${PYTHON_BIN}" ]]; then
  echo "[job] python=${PYTHON_BIN}"
fi
echo "[job] oddsketch=${ODDSKETCH_BIN}"
echo "[job] uv_cache_dir=${UV_CACHE_DIR}"

"${UV_BIN}" sync
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python "${TASK_DIR}/scripts/prepare_thread_sweep_subset.py"
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python "${TASK_DIR}/scripts/run_thread_sweep.py"

echo "[job] end=$(date)"
