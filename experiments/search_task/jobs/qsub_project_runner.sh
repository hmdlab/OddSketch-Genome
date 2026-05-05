#!/bin/bash
#$ -cwd
#$ -V
#$ -N search_project_runner
#$ -q tsmall
#$ -l mem_req=16g
#$ -l h_vmem=16g

set -euo pipefail

TASK_DIR=$(pwd)
if [[ ! -f "${TASK_DIR}/config.json" || ! -d "${TASK_DIR}/scripts" || ! -d "${TASK_DIR}/jobs" ]]; then
  echo "Run this script from experiments/search_task:" >&2
  echo "  cd /path/to/genome-oddsketch/experiments/search_task && qsub jobs/qsub_project_runner.sh" >&2
  exit 1
fi

REPO_ROOT=$(cd "${TASK_DIR}/../.." && pwd)

CONFIG_PATH=${1:-"${TASK_DIR}/config.json"}
if [[ "${CONFIG_PATH}" != /* ]]; then
  CONFIG_PATH="${TASK_DIR}/${CONFIG_PATH}"
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
  PYTHON_BIN_DEFAULT=$(command -v python3.11 || true)
fi
PYTHON_BIN=${PYTHON_BIN:-"${PYTHON_BIN_DEFAULT}"}

ODDSKETCH_BIN=${ODDSKETCH_BIN:-"${REPO_ROOT}/src/oddsketch"}
if [[ ! -x "${ODDSKETCH_BIN}" ]]; then
  echo "warning: oddsketch binary not found at ${ODDSKETCH_BIN}" >&2
fi

export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${USER}"}
export MPLBACKEND=${MPLBACKEND:-Agg}
export ODDSKETCH_BIN

mkdir -p "${UV_CACHE_DIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] task_dir=${TASK_DIR}"
echo "[job] config=${CONFIG_PATH}"
echo "[job] uv=${UV_BIN}"
if [[ -n "${PYTHON_BIN}" ]]; then
  echo "[job] python=${PYTHON_BIN}"
fi
echo "[job] oddsketch=${ODDSKETCH_BIN}"
echo "[job] uv_cache_dir=${UV_CACHE_DIR}"

cd "${REPO_ROOT}"
echo "[job] syncing uv environment"
"${UV_BIN}" sync
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python "${TASK_DIR}/scripts/project_runner.py" --config "${CONFIG_PATH}"

echo "[job] end=$(date)"
