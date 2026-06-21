#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
OUT_PATH=${1:-/tmp/genome-oddsketch-public.tar.gz}
PREFIX=${PUBLIC_ARCHIVE_PREFIX:-genome-oddsketch/}

TMP_LIST=$(mktemp "${TMPDIR:-/tmp}/genome-oddsketch-public-files.XXXXXX")
cleanup() {
  rm -f "${TMP_LIST}"
}
trap cleanup EXIT

cd "${REPO_ROOT}"
mkdir -p "$(dirname "${OUT_PATH}")"

git ls-files -z --cached --others --exclude-standard > "${TMP_LIST}"
tar --null -czf "${OUT_PATH}" --transform "s|^|${PREFIX}|" -T "${TMP_LIST}"

echo "[archive] wrote ${OUT_PATH}"
