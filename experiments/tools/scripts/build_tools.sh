#!/usr/bin/env bash
set -euo pipefail

METHOD="auto"
BINDASH_REF="v2.2.0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method)
      METHOD="$2"
      shift 2
      ;;
    --ref)
      BINDASH_REF="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${TOOLS_DIR}/bin"
SRC_DIR="${TOOLS_DIR}/bindash"
BUILD_DIR="${SRC_DIR}/build"

mkdir -p "${BIN_DIR}" "${TOOLS_DIR}"

echo "[build_tools] BinDash is an external tool by its original authors."
echo "[build_tools] This repository does not vendor BinDash source or binaries."

if command -v bindash >/dev/null 2>&1; then
  echo "[build_tools] bindash found on PATH: $(command -v bindash)"
  exit 0
fi

install_from_source() {
  command -v git >/dev/null
  command -v cmake >/dev/null
  command -v make >/dev/null

  if [[ ! -d "${SRC_DIR}/.git" ]]; then
    git clone https://github.com/jianshu93/bindash.git "${SRC_DIR}"
  fi

  git -C "${SRC_DIR}" fetch --tags --all
  git -C "${SRC_DIR}" checkout "${BINDASH_REF}"

  cmake -S "${SRC_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" -j

  BUILT_BIN="$(find "${BUILD_DIR}" -type f -name bindash | head -n 1 || true)"
  if [[ -z "${BUILT_BIN}" ]]; then
    echo "[build_tools] Could not find built bindash binary in ${BUILD_DIR}" >&2
    exit 1
  fi

  cp "${BUILT_BIN}" "${BIN_DIR}/bindash"
  chmod +x "${BIN_DIR}/bindash"
  echo "[build_tools] Installed bindash to ${BIN_DIR}/bindash"
}

if [[ "${METHOD}" == "auto" ]]; then
  if [[ "$(uname -s)" == "Darwin" ]] && command -v brew >/dev/null 2>&1; then
    echo "[build_tools] Using Homebrew path"
    brew tap jianshu93/bindash
    brew update
    brew install --cc=gcc-13 bindash
    exit 0
  fi

  if command -v conda >/dev/null 2>&1; then
    echo "[build_tools] Using conda path"
    conda install -y bindash -c bioconda
    exit 0
  fi

  echo "[build_tools] Falling back to source build"
  install_from_source
  exit 0
fi

case "${METHOD}" in
  brew)
    command -v brew >/dev/null
    brew tap jianshu93/bindash
    brew update
    brew install --cc=gcc-13 bindash
    ;;
  conda)
    command -v conda >/dev/null
    conda install -y bindash -c bioconda
    ;;
  source)
    install_from_source
    ;;
  none)
    echo "[build_tools] Skipped bindash install by request"
    ;;
  *)
    echo "Unsupported method: ${METHOD}" >&2
    exit 1
    ;;
esac
