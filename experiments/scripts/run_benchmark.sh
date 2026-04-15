#!/usr/bin/env bash
set -euo pipefail

MODE="all"
SKIP_BINDASH=0
PAIR_CONFIG=""
SEARCH_CONFIG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    --skip-bindash)
      SKIP_BINDASH=1
      shift
      ;;
    --pair-config)
      PAIR_CONFIG="$2"
      shift 2
      ;;
    --search-config)
      SEARCH_CONFIG="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${EXP_DIR}/.." && pwd)"

PAIR_DIR="${EXP_DIR}/pair_task"
SEARCH_DIR="${EXP_DIR}/search_task"
OUT_DIR="${EXP_DIR}/outputs"

if [[ -z "${PAIR_CONFIG}" ]]; then
  PAIR_CONFIG="${EXP_DIR}/configs/pair_task.pipeline_config.json"
fi
if [[ -z "${SEARCH_CONFIG}" ]]; then
  SEARCH_CONFIG="${EXP_DIR}/configs/search_task.config.json"
fi

mkdir -p "${OUT_DIR}" "${OUT_DIR}/pair_task" "${OUT_DIR}/search_task"

echo "[benchmark] sync python environment"
uv sync

echo "[benchmark] build oddsketch binaries"
make -C "${ROOT_DIR}/src"

export ODDSKETCH_BIN="${ODDSKETCH_BIN:-${ROOT_DIR}/src/oddsketch}"

echo "[benchmark] ODDSKETCH_BIN=${ODDSKETCH_BIN}"

if [[ "${MODE}" == "pair" || "${MODE}" == "all" ]]; then
  echo "[benchmark] running pair_task"
  pushd "${PAIR_DIR}" >/dev/null
  uv run python make_genomes/make_diverse_genomes.py --config "${PAIR_CONFIG}"
  uv run python cal/cal_jaccard_true.py --config "${PAIR_CONFIG}"
  uv run python cal/cal_jaccard_oddsketch.py --config "${PAIR_CONFIG}"
  if [[ "${SKIP_BINDASH}" -eq 0 ]]; then
    uv run python cal/cal_jaccard_bindash.py --config "${PAIR_CONFIG}"
  else
    echo "[benchmark] skip bindash for pair_task"
  fi
  popd >/dev/null

  uv run python "${EXP_DIR}/scripts/make_figures.py" --task pair --exp-root "${EXP_DIR}"

  cp -f "${PAIR_DIR}/data/test_genomes/jaccard_true_results.txt" "${OUT_DIR}/pair_task/" 2>/dev/null || true
  cp -f "${PAIR_DIR}/data/test_genomes/jaccard_oddsketch_results.txt" "${OUT_DIR}/pair_task/" 2>/dev/null || true
  cp -f "${PAIR_DIR}/data/test_genomes/jaccard_bindash_results.txt" "${OUT_DIR}/pair_task/" 2>/dev/null || true
  cp -f "${PAIR_DIR}/data/test_genomes/comparison_results_oddsketch.csv" "${OUT_DIR}/pair_task/" 2>/dev/null || true
  cp -f "${PAIR_DIR}/data/test_genomes/comparison_results_bindash.csv" "${OUT_DIR}/pair_task/" 2>/dev/null || true
fi

if [[ "${MODE}" == "search" || "${MODE}" == "all" ]]; then
  echo "[benchmark] running search_task"
  uv run python "${SEARCH_DIR}/project_runner.py" --config "${SEARCH_CONFIG}"
  uv run python "${EXP_DIR}/scripts/make_figures.py" --task search --exp-root "${EXP_DIR}"
fi

echo "[benchmark] done"
echo "[benchmark] outputs: ${OUT_DIR}"
