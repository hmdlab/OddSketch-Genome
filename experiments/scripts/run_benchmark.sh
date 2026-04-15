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
if [[ -z "${PAIR_CONFIG}" ]]; then
  PAIR_CONFIG="${PAIR_DIR}/config.json"
fi
if [[ -z "${SEARCH_CONFIG}" ]]; then
  SEARCH_CONFIG="${SEARCH_DIR}/config.json"
fi

echo "[benchmark] sync python environment"
uv sync

echo "[benchmark] build oddsketch binaries"
make -C "${ROOT_DIR}/src"

export ODDSKETCH_BIN="${ODDSKETCH_BIN:-${ROOT_DIR}/src/oddsketch}"

echo "[benchmark] ODDSKETCH_BIN=${ODDSKETCH_BIN}"

if [[ "${MODE}" == "pair" || "${MODE}" == "all" ]]; then
  echo "[benchmark] running pair_task"
  pushd "${PAIR_DIR}" >/dev/null
  uv run python scripts/make_genomes.py --config "${PAIR_CONFIG}"
  uv run python scripts/cal_jaccard_true.py --config "${PAIR_CONFIG}"
  uv run python scripts/cal_jaccard_oddsketch.py --config "${PAIR_CONFIG}"
  if [[ "${SKIP_BINDASH}" -eq 0 ]]; then
    uv run python scripts/cal_jaccard_bindash.py --config "${PAIR_CONFIG}"
  else
    echo "[benchmark] skip bindash for pair_task"
  fi
  popd >/dev/null

  uv run python "${EXP_DIR}/scripts/make_figures.py" --task pair --exp-root "${EXP_DIR}"
fi

if [[ "${MODE}" == "search" || "${MODE}" == "all" ]]; then
  echo "[benchmark] running search_task"
  if [[ "${SKIP_BINDASH}" -eq 0 ]]; then
    uv run python "${SEARCH_DIR}/scripts/project_runner.py" --config "${SEARCH_CONFIG}"
  else
    uv run python "${SEARCH_DIR}/scripts/project_runner.py" --config "${SEARCH_CONFIG}" --skip-bindash
  fi
  uv run python "${EXP_DIR}/scripts/make_figures.py" --task search --exp-root "${EXP_DIR}"
fi

echo "[benchmark] done"
echo "[benchmark] pair outputs: ${PAIR_DIR}/outputs"
echo "[benchmark] search outputs: ${SEARCH_DIR}/outputs"
