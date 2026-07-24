#!/bin/bash
#$ -cwd
#$ -V
#$ -N pair_project_runner
#$ -q h.q
#$ -pe OpenMP 8
#$ -l mem_req=16g
#$ -l h_vmem=16g

set -euo pipefail

SUBMIT_DIR=${SGE_O_WORKDIR:-$(pwd)}
if [[ -n "${GENOME_ODDSKETCH_ROOT:-}" ]]; then
  REPO_ROOT=${GENOME_ODDSKETCH_ROOT}
elif [[ -f "${SUBMIT_DIR}/experiments/pair_task/configs/sketchsize_repeats/config.json" ]]; then
  REPO_ROOT=${SUBMIT_DIR}
elif [[ -f "${SUBMIT_DIR}/configs/sketchsize_repeats/config.json" ]]; then
  REPO_ROOT=$(cd "${SUBMIT_DIR}/../.." && pwd)
else
  REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
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
export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${USER}"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-"${TMPDIR:-/tmp}/matplotlib-${USER}"}
export MPLBACKEND=${MPLBACKEND:-Agg}
mkdir -p "${UV_CACHE_DIR}" "${MPLCONFIGDIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

MODE=paired
if [[ "${1:-}" == "--legacy-batch" ]]; then
  MODE=legacy
  shift
fi

echo "[job] mode=${MODE}"
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

if [[ "${MODE}" == "legacy" ]]; then
  LEGACY_ARGS=("$@")
  if [[ ${#LEGACY_ARGS[@]} -eq 0 ]]; then
    LEGACY_ARGS=(--config-dir configs/default --jobs "${PAIR_TASK_JOBS:-${NSLOTS:-1}}")
  fi
  "${UV_BIN}" "${UV_RUN_ARGS[@]}" python \
    "${TASK_DIR}/scripts/batch_project_runner.py" \
    "${LEGACY_ARGS[@]}"
else
  if [[ ! -x "${BINDASH_BIN}" ]]; then
    echo "BinDash executable not found: ${BINDASH_BIN}" >&2
    exit 1
  fi
  PAIRED_ARGS=(
    --config "${TASK_DIR}/configs/sketchsize_repeats/config.json"
  )
  if [[ -n "${PAIRED_SKETCHSIZE_RUN_DIR:-}" ]]; then
    PAIRED_ARGS+=(--run-dir "${PAIRED_SKETCHSIZE_RUN_DIR}")
  fi
  PAIRED_ARGS+=("$@")
  echo "[job] oddsketch=${ODDSKETCH_BIN}"
  echo "[job] true_jaccard=${TRUE_JACCARD_BIN}"
  echo "[job] bindash=$("${BINDASH_BIN}" --version 2>&1 | head -n 1)"
  echo "[job] resume_dir=${PAIRED_SKETCHSIZE_RUN_DIR:-new run}"
  # The paired runner calls analysis/aggregate/analyze_paired_sketchsize_repeats.py
  # after aggregation and writes the final figure to:
  #   <run_dir>/summary/RMSE_by_true_jaccard_panels.png
  "${UV_BIN}" "${UV_RUN_ARGS[@]}" python \
    "${TASK_DIR}/scripts/run_paired_sketchsize_repeats.py" \
    "${PAIRED_ARGS[@]}"
fi

echo "[job] end=$(date)"
