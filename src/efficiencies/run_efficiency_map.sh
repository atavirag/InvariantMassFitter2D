#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$(realpath -m "$SCRIPT_DIR/../../build/efficiency_map")"

mkdir -p "$BUILD_DIR"
rm -f "$BUILD_DIR"/*_cxx.* "$BUILD_DIR"/*_rdict.pcm

cd "$SCRIPT_DIR" || exit 1
root -b computeEfficiencyMap.cxx+