#!/bin/bash
#$ -cwd
#$ -N pair_supplementary
#$ -q h.q
#$ -pe OpenMP 8
#$ -l mem_req=16g
#$ -l h_vmem=16g

set -euo pipefail

usage() {
  echo "usage: supplementary.sh [bindash_recommended|k_sensitivity|oph] [runner options...]" >&2
}

EXPERIMENT=all
if [[ $# -gt 0 ]]; then
  case "$1" in
    bindash_recommended|k_sensitivity|oph)
      EXPERIMENT=$1
      shift
      ;;
    -*)
      ;;
    *)
      usage
      exit 2
      ;;
  esac
fi

SUBMIT_DIR=${SGE_O_WORKDIR:-$(pwd)}
if [[ -f "${SUBMIT_DIR}/experiments/pair_task/configs/supplementary.json" ]]; then
  REPO_ROOT=${SUBMIT_DIR}
elif [[ -f "${SUBMIT_DIR}/configs/supplementary.json" ]]; then
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
export ODDSKETCH_FIGURE_LABEL=${ODDSKETCH_FIGURE_LABEL:-OddSketch-Genome}
JOB_USER=${USER:-$(id -un)}
export UV_CACHE_DIR=${UV_CACHE_DIR:-"${TMPDIR:-/tmp}/uv-cache-${JOB_USER}"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-"${TMPDIR:-/tmp}/matplotlib-${JOB_USER}"}
export MPLBACKEND=${MPLBACKEND:-Agg}
mkdir -p "${UV_CACHE_DIR}" "${MPLCONFIGDIR}"

UV_RUN_ARGS=(run --no-sync)
if [[ -n "${PYTHON_BIN}" ]]; then
  UV_RUN_ARGS+=(--python "${PYTHON_BIN}")
fi

echo "[job] experiment=${EXPERIMENT}"
echo "[job] host=$(hostname)"
echo "[job] start=$(date)"
echo "[job] repo_root=${REPO_ROOT}"
echo "[job] slots=${NSLOTS:-8}"
echo "[job] uv=${UV_BIN}"

cd "${REPO_ROOT}"
make -C src CXX="${CXX:-g++}" LDFLAGS="${LDFLAGS:--lstdc++fs}"
"${UV_BIN}" sync

if [[ "${EXPERIMENT}" == "all" || "${EXPERIMENT}" == "bindash_recommended" ]]; then
  if [[ ! -x "${BINDASH_BIN}" ]]; then
    echo "BinDash executable not found: ${BINDASH_BIN}" >&2
    exit 1
  fi
fi

RUN_STAMP=$(date +%Y%m%d_%H%M%S)
RUN_SUFFIX=${JOB_ID:-$$}
VALIDATION_ROOT="${TASK_DIR}/outputs/validation"
RUN_DIR="${VALIDATION_ROOT}/run_${RUN_STAMP}_${RUN_SUFFIX}"
mkdir -p "${VALIDATION_ROOT}"
mkdir "${RUN_DIR}"
echo "[job] run_dir=${RUN_DIR}"

run_bindash_recommended() {
  "${UV_BIN}" "${UV_RUN_ARGS[@]}" python \
    "${TASK_DIR}/scripts/runners/run_paired.py" \
    --config "${TASK_DIR}/configs/supplementary.json" \
    --experiment bindash_recommended \
    --output-dir "${RUN_DIR}/bindash_recommended" \
    "$@"
}

run_validation() {
  local experiment_name=$1
  local output_name=$2
  shift 2
  "${UV_BIN}" "${UV_RUN_ARGS[@]}" python \
    "${TASK_DIR}/scripts/runners/run_validation.py" \
    --experiment "${experiment_name}" \
    --config "${TASK_DIR}/configs/supplementary.json" \
    --output-dir "${RUN_DIR}/${output_name}" \
    "$@"
}

case "${EXPERIMENT}" in
  all)
    run_bindash_recommended "$@"
    run_validation k_sensitivity k_sensitivity "$@"
    run_validation oph oph "$@"
    ;;
  bindash_recommended)
    run_bindash_recommended "$@"
    ;;
  k_sensitivity)
    run_validation k_sensitivity k_sensitivity "$@"
    ;;
  oph)
    run_validation oph oph "$@"
    ;;
esac

echo "[job] end=$(date)"
