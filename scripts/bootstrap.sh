#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
TOOLS_BIN_DIR="${REPO_ROOT}/experiments/tools/bin"
INSTALL_PATH="${TOOLS_BIN_DIR}/bindash"

BINDASH_REPO=${BINDASH_REPO:-https://github.com/zhaoxiaofei/bindash.git}
BINDASH_TAG=${BINDASH_TAG:-v2.6}

echo "[bootstrap] repo_root=${REPO_ROOT}"
echo "[bootstrap] bindash_repo=${BINDASH_REPO}"
echo "[bootstrap] bindash_tag=${BINDASH_TAG}"

mkdir -p "${TOOLS_BIN_DIR}"

if [[ -x "${INSTALL_PATH}" ]]; then
  echo "[bootstrap] bindash already installed: ${INSTALL_PATH}"
  "${INSTALL_PATH}" --help >/dev/null 2>&1 || true
  exit 0
fi

for cmd in git cmake make; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "[bootstrap] required command not found: ${cmd}" >&2
    exit 1
  fi
done

if command -v g++ >/dev/null 2>&1; then
  export CXX=${CXX:-g++}
elif command -v clang++ >/dev/null 2>&1; then
  export CXX=${CXX:-clang++}
else
  echo "[bootstrap] no C++ compiler found (need g++ or clang++)" >&2
  exit 1
fi

BUILD_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/bindash-bootstrap.XXXXXX")
cleanup() {
  rm -rf "${BUILD_ROOT}"
}
trap cleanup EXIT

echo "[bootstrap] build_root=${BUILD_ROOT}"
git clone --depth 1 --branch "${BINDASH_TAG}" "${BINDASH_REPO}" "${BUILD_ROOT}/src"

cmake -S "${BUILD_ROOT}/src" -B "${BUILD_ROOT}/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_ROOT}/build" -j "${CMAKE_BUILD_PARALLEL_LEVEL:-1}"

BIN_CANDIDATE=$(find "${BUILD_ROOT}/build" -maxdepth 3 -type f -name bindash | head -n 1 || true)
if [[ -z "${BIN_CANDIDATE}" || ! -f "${BIN_CANDIDATE}" ]]; then
  echo "[bootstrap] built bindash binary not found" >&2
  exit 1
fi

cp "${BIN_CANDIDATE}" "${INSTALL_PATH}"
chmod +x "${INSTALL_PATH}"

echo "[bootstrap] installed bindash -> ${INSTALL_PATH}"
"${INSTALL_PATH}" --help >/dev/null 2>&1 || true
echo "[bootstrap] done"
