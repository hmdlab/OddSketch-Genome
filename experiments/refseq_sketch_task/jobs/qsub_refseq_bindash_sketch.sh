#!/bin/bash
#$ -cwd
#$ -V
#$ -N refseq_bindash
#$ -q h.q
#$ -pe OpenMP 8
#$ -l mem_req=16g
#$ -l h_vmem=16g

set -euo pipefail

REPO_ROOT=$(pwd)
if [[ ! -f "${REPO_ROOT}/README.md" || ! -d "${REPO_ROOT}/experiments/refseq_sketch_task" ]]; then
  echo "Run this job from the repository root:" >&2
  echo "  cd /path/to/genome-oddsketch && qsub experiments/refseq_sketch_task/jobs/qsub_refseq_bindash_sketch.sh" >&2
  exit 1
fi

TASK_DIR="${REPO_ROOT}/experiments/refseq_sketch_task"
CONFIG_PATH=${1:-"${TASK_DIR}/config.json"}
if [[ $# -gt 0 ]]; then
  shift
fi
if [[ "${CONFIG_PATH}" != /* ]]; then
  CONFIG_PATH="${REPO_ROOT}/${CONFIG_PATH}"
fi
if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "config not found: ${CONFIG_PATH}" >&2
  exit 1
fi

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
  PYTHON_BIN_DEFAULT=$(command -v python3.11 || command -v python3 || true)
fi
PYTHON_BIN=${PYTHON_BIN:-"${PYTHON_BIN_DEFAULT}"}

BINDASH_BIN=${BINDASH_BIN:-"${REPO_ROOT}/experiments/tools/bin/bindash"}
if [[ ! -x "${BINDASH_BIN}" ]]; then
  BINDASH_BIN=$(command -v bindash || true)
fi
if [[ -z "${BINDASH_BIN}" || ! -x "${BINDASH_BIN}" ]]; then
  echo "bindash not found. Set BINDASH_BIN or run scripts/bootstrap.sh first." >&2
  exit 1
fi

export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${USER}"}
export BINDASH_BIN
mkdir -p "${UV_CACHE_DIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] config=${CONFIG_PATH}"
if [[ $# -gt 0 ]]; then
  echo "[job] runner_args=$*"
fi
echo "[job] uv=${UV_BIN}"
if [[ -n "${PYTHON_BIN}" ]]; then
  echo "[job] python=${PYTHON_BIN}"
fi
echo "[job] bindash=${BINDASH_BIN}"
echo "[job] uv_cache_dir=${UV_CACHE_DIR}"

"${UV_BIN}" sync
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python "${TASK_DIR}/scripts/refseq_bindash_sketch_runner.py" --config "${CONFIG_PATH}" "$@"

echo "[job] end=$(date)"
