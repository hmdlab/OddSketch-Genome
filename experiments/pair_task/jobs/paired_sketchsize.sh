#!/bin/bash
#$ -cwd
#$ -N paired_sketchsize
#$ -q h.q
#$ -pe OpenMP 8
#$ -l mem_req=16g
#$ -l h_vmem=16g

set -euo pipefail

export ODDSKETCH_FIGURE_LABEL=${ODDSKETCH_FIGURE_LABEL:-OddSketch-Genome}

SUBMIT_DIR=${SGE_O_WORKDIR:-$(pwd)}
if [[ -f "${SUBMIT_DIR}/experiments/pair_task/configs/paired.json" ]]; then
  REPO_ROOT=${SUBMIT_DIR}
elif [[ -f "${SUBMIT_DIR}/configs/paired.json" ]]; then
  REPO_ROOT=$(cd "${SUBMIT_DIR}/../.." && pwd)
else
  echo "Submit from the repository root or experiments/pair_task." >&2
  exit 1
fi
TASK_DIR="${REPO_ROOT}/experiments/pair_task"

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
TRUE_JACCARD_BIN=${TRUE_JACCARD_BIN:-"${REPO_ROOT}/experiments/tools/bin/true_jaccard"}
BINDASH_BIN=${BINDASH_BIN:-"${REPO_ROOT}/experiments/tools/bin/bindash"}
export ODDSKETCH_BIN TRUE_JACCARD_BIN BINDASH_BIN
JOB_USER=${USER:-$(id -un)}
export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${JOB_USER}"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-"${TMPDIR:-/tmp}/matplotlib-${JOB_USER}"}
export MPLBACKEND=${MPLBACKEND:-Agg}
mkdir -p "${UV_CACHE_DIR}" "${MPLCONFIGDIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] task_dir=${TASK_DIR}"
echo "[job] slots=${NSLOTS:-8}"
echo "[job] uv=${UV_BIN}"
echo "[job] uv_cache_dir=${UV_CACHE_DIR}"

cd "${REPO_ROOT}"
make -C src CXX="${CXX:-g++}" LDFLAGS="${LDFLAGS:--lstdc++fs}"
"${UV_BIN}" sync

if [[ ! -x "${BINDASH_BIN}" ]]; then
  echo "BinDash executable not found: ${BINDASH_BIN}" >&2
  exit 1
fi

echo "[job] oddsketch=${ODDSKETCH_BIN}"
echo "[job] true_jaccard=${TRUE_JACCARD_BIN}"
echo "[job] bindash=$("${BINDASH_BIN}" --version 2>&1 | head -n 1)"
"${UV_BIN}" "${UV_RUN_ARGS[@]}" python \
  "${TASK_DIR}/scripts/runners/run_paired.py" \
  --config "${TASK_DIR}/configs/paired.json" \
  --experiment paired \
  "$@"

echo "[job] end=$(date)"
