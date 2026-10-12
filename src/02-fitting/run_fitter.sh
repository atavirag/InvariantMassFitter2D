#!/bin/bash

# --- Resolve absolute paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR"
BUILD_DIR="$(realpath ../../build)"

# --- Prepare build directory ---
mkdir -p "$BUILD_DIR"
rm -f "$BUILD_DIR"/*_cxx.* "$BUILD_DIR"/*_rdict.pcm

# --- Run ROOT ---
cd "$SRC_DIR" || exit 1
root -b <<EOF
gErrorIgnoreLevel = kWarning;
gSystem->SetBuildDir("$BUILD_DIR", true);
gSystem->AddIncludePath("-I/Users/atavirag/Codigo/Analysis/Double-D/InvariantMassFitter2D/src/01-prefiltering");
.L functions/BifurcatedCB.cxx+
.L functions/SignalFunctions.cxx+
.L functions/BackgroundFunctions.cxx+
.L functions/ReflectedFunctions.cxx+
.L DatasetCreator.cxx+
.L InvariantMassFitter.cxx+
.L runFitter.cxx+

runFitter(0, nullptr);
.q
EOF
