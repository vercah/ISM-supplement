#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   bash install_concorde_here.sh
#
# Installs into:
#   ./concorde
#   ./linkern
#
# Supported:
#   - macOS arm64
#   - macOS x86_64
#   - Linux x86_64

OUT_DIR="$(pwd)"
KEEP_BUILD="${KEEP_BUILD:-0}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

fetch() {
  local url="$1"
  local out="$2"

  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "$url" -o "$out"
  elif command -v wget >/dev/null 2>&1; then
    wget -qO "$out" "$url"
  else
    echo "Need either curl or wget" >&2
    exit 1
  fi
}

cpu_count() {
  if command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null && return 0
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc 2>/dev/null && return 0
  fi
  if command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.logicalcpu 2>/dev/null && return 0
  fi
  echo 2
}

OS="$(uname -s)"
ARCH="$(uname -m)"

CONCORDE_SRC_URL="https://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz"

CONFIGURE_EXTRA=()
case "${OS}:${ARCH}" in
  Darwin:arm64)
    QSOPT_A_URL="https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/m1/qsopt.a"
    QSOPT_H_URL="https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/m1/qsopt.h"
    CONFIGURE_EXTRA=(--host=darwin)
    ;;
  Darwin:x86_64)
    QSOPT_A_URL="https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/mac64/qsopt.a"
    QSOPT_H_URL="https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/mac64/qsopt.h"
    CONFIGURE_EXTRA=(--host=darwin)
    ;;
  Linux:x86_64)
    QSOPT_A_URL="https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/PIC/qsopt.PIC.a"
    QSOPT_H_URL="https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/PIC/qsopt.h"
    ;;
  *)
    echo "Unsupported platform: ${OS} ${ARCH}" >&2
    echo "Supported platforms: macOS arm64, macOS x86_64, Linux x86_64" >&2
    exit 1
    ;;
esac

need_cmd tar
need_cmd make
need_cmd cp
need_cmd uname
need_cmd mktemp

WORK_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/concorde-build.XXXXXX")"
cleanup() {
  if [[ "$KEEP_BUILD" == "1" ]]; then
    echo "Keeping build directory: $WORK_ROOT"
  else
    rm -rf "$WORK_ROOT"
  fi
}
trap cleanup EXIT

QSOPT_DIR="$WORK_ROOT/qsopt"
SRC_DIR="$WORK_ROOT/src"
mkdir -p "$QSOPT_DIR" "$SRC_DIR"

echo "Downloading QSopt for ${OS} ${ARCH} ..."
fetch "$QSOPT_A_URL" "$QSOPT_DIR/qsopt.a"
fetch "$QSOPT_H_URL" "$QSOPT_DIR/qsopt.h"

echo "Downloading Concorde source ..."
fetch "$CONCORDE_SRC_URL" "$WORK_ROOT/concorde.tgz"

echo "Extracting Concorde ..."
tar -xzf "$WORK_ROOT/concorde.tgz" -C "$SRC_DIR"

CONCORDE_DIR="$SRC_DIR/concorde"
if [[ ! -d "$CONCORDE_DIR" ]]; then
  echo "Expected source directory not found: $CONCORDE_DIR" >&2
  exit 1
fi

JOBS="$(cpu_count)"
export CFLAGS="${CFLAGS:--fPIC -O3 -g -std=gnu89}"

echo "Configuring Concorde ..."
(
  cd "$CONCORDE_DIR"
  ./configure --with-qsopt="$QSOPT_DIR" "${CONFIGURE_EXTRA[@]}"
)

echo "Building Concorde with $JOBS job(s) ..."
(
  cd "$CONCORDE_DIR"
  make -j"$JOBS"
)

if [[ ! -x "$CONCORDE_DIR/TSP/concorde" ]]; then
  echo "Build finished, but TSP/concorde was not produced." >&2
  exit 1
fi

cp "$CONCORDE_DIR/TSP/concorde" "$OUT_DIR/concorde"
chmod +x "$OUT_DIR/concorde"

if [[ -x "$CONCORDE_DIR/LINKERN/linkern" ]]; then
  cp "$CONCORDE_DIR/LINKERN/linkern" "$OUT_DIR/linkern"
  chmod +x "$OUT_DIR/linkern"
fi

echo
echo "Installed in current directory:"
echo "  $OUT_DIR/concorde"
if [[ -x "$OUT_DIR/linkern" ]]; then
  echo "  $OUT_DIR/linkern"
fi
echo
echo "Smoke test:"
echo "  ./concorde -s 99 -k 100"
