#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$(realpath -m "$SCRIPT_DIR/../../build/splot")"   # separate dir: does not wipe the fitter build

mkdir -p "$BUILD_DIR"
rm -f "$BUILD_DIR"/*_cxx.* "$BUILD_DIR"/*_rdict.pcm

cd "$SCRIPT_DIR" || exit 1
root -b <<EOF
gErrorIgnoreLevel = kWarning;
gSystem->SetBuildDir("$BUILD_DIR", true);
.L SPlotter.cxx+
.L runSPlot.cxx+

runSPlot(0, nullptr);
.q
EOF