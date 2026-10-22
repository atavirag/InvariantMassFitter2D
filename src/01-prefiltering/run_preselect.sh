#!/bin/bash

set -e

# Directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SRC_DIR="$SCRIPT_DIR"
BUILD_DIR="$(realpath ../../build)"

# --- Prepare build directory ---
mkdir -p "$BUILD_DIR"
rm -f "$BUILD_DIR"/*_cxx.* "$BUILD_DIR"/*_rdict.pcm

cd "$SRC_DIR"

echo "Source dir : $SRC_DIR"
echo "Build dir  : $BUILD_DIR"

root -l -b <<EOF
gErrorIgnoreLevel = kWarning;
gSystem->SetBuildDir("$BUILD_DIR", true);
gSystem->AddIncludePath("-I$SRC_DIR");

.L TreeReader.cxx+
.L TreeSelector.cxx+
.x preselect.cpp+

.q
EOF