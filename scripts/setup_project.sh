#!/bin/bash
# setup_project.sh - Organize your analysis into proper directory structure

set -e

echo "========================================="
echo "  D Meson Analysis - Project Setup"
echo "========================================="
echo ""

# Create directory structure
echo "Creating directory structure..."
mkdir -p src
mkdir -p build
mkdir -p data
mkdir -p plots
mkdir -p python
mkdir -p scripts
mkdir -p docs

echo "✓ Directories created"
echo ""

# Move files to appropriate locations
echo "Organizing files..."

# Move C++ source files to src/
if [ -f "TreeSelector.h" ]; then
    mv TreeSelector.h src/ 2>/dev/null || true
    echo "  Moved TreeSelector.h → src/"
fi
if [ -f "TreeSelector.cxx" ]; then
    mv TreeSelector.cxx src/ 2>/dev/null || true
    echo "  Moved TreeSelector.cxx → src/"
fi
if [ -f "InvariantMassFitter.h" ]; then
    mv InvariantMassFitter.h src/ 2>/dev/null || true
    echo "  Moved InvariantMassFitter.h → src/"
fi
if [ -f "InvariantMassFitter.cxx" ]; then
    mv InvariantMassFitter.cxx src/ 2>/dev/null || true
    echo "  Moved InvariantMassFitter.cxx → src/"
fi
if [ -f "CandidatePair.h" ]; then
    mv CandidatePair.h src/ 2>/dev/null || true
    echo "  Moved CandidatePair.h → src/"
fi
if [ -f "runFit.cxx" ]; then
    mv runFit.cxx src/ 2>/dev/null || true
    echo "  Moved runFit.cxx → src/"
fi
if [ -f "advancedFitExample.cxx" ]; then
    mv advancedFitExample.cxx src/ 2>/dev/null || true
    echo "  Moved advancedFitExample.cxx → src/"
fi

# Move Python files
if [ -f "visualize_fits.py" ]; then
    mv visualize_fits.py python/ 2>/dev/null || true
    echo "  Moved visualize_fits.py → python/"
fi
if [ -f "visualize_with_fits.py" ]; then
    mv visualize_with_fits.py python/ 2>/dev/null || true
    echo "  Moved visualize_with_fits.py → python/"
fi
if [ -f "check_analysis.py" ]; then
    mv check_analysis.py python/ 2>/dev/null || true
    echo "  Moved check_analysis.py → python/"
fi
if [ -f "requirements.txt" ]; then
    cp requirements.txt python/ 2>/dev/null || true
fi

# Move scripts
if [ -f "run_complete_analysis.sh" ]; then
    mv run_complete_analysis.sh scripts/ 2>/dev/null || true
    chmod +x scripts/run_complete_analysis.sh
    echo "  Moved run_complete_analysis.sh → scripts/"
fi

# Move documentation
for doc in README*.md PROJECT_OVERVIEW.md; do
    if [ -f "$doc" ]; then
        mv "$doc" docs/ 2>/dev/null || true
        echo "  Moved $doc → docs/"
    fi
done

# Move build artifacts
echo ""
echo "Cleaning up build artifacts..."
mv *.so build/ 2>/dev/null || true
mv *.d build/ 2>/dev/null || true
mv *.o build/ 2>/dev/null || true
mv *_cxx.* build/ 2>/dev/null || true
mv *_cpp.* build/ 2>/dev/null || true
echo "✓ Build artifacts moved to build/"

# Move data files
echo ""
echo "Moving data files..."
mv *.root data/ 2>/dev/null || true
if [ "$(ls -A data/*.root 2>/dev/null)" ]; then
    echo "✓ ROOT files moved to data/"
fi

# Move plots
mv *.pdf plots/ 2>/dev/null || true
mv *.png plots/ 2>/dev/null || true
if [ "$(ls -A plots/ 2>/dev/null)" ]; then
    echo "✓ Plots moved to plots/"
fi

echo ""
echo "========================================="
echo "  Project organization complete!"
echo "========================================="
echo ""
echo "Directory structure:"
echo "  src/       - C++ source code"
echo "  build/     - Build artifacts (.so, .d, .o)"
echo "  data/      - ROOT data files"
echo "  plots/     - Generated plots"
echo "  python/    - Python visualization scripts"
echo "  scripts/   - Shell scripts"
echo "  docs/      - Documentation"
echo ""
echo "Next steps:"
echo "  1. Update Makefile paths if needed"
echo "  2. Review .gitignore"
echo "  3. Initialize git: git init && git add . && git commit -m 'Initial commit'"
echo ""
echo "To build: make"
echo "To run:   make run"
echo "To clean: make clean"
echo ""