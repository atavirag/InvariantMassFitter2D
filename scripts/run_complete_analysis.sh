#!/bin/bash
# run_complete_analysis.sh
# Complete D meson pair analysis pipeline: C++ fitting + Python visualization

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  D Meson Pair Analysis Pipeline${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Configuration
OS_FILE="os_pairs.root"
LS_FILE="ls_pairs.root"
WORKSPACE="fit_workspace.root"
PLOT_DIR="analysis_plots"
USE_PYROOT=true  # Set to false if PyROOT not available

# Step 1: Compile C++ code
echo -e "${YELLOW}[1/4] Compiling C++ code...${NC}"
if [ ! -f "Makefile" ]; then
    echo -e "${RED}Error: Makefile not found!${NC}"
    exit 1
fi

make clean
make

if [ $? -ne 0 ]; then
    echo -e "${RED}Compilation failed!${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Compilation successful${NC}"
echo ""

# Step 2: Run C++ fitting
echo -e "${YELLOW}[2/4] Running RooFit analysis...${NC}"
if [ ! -f "./runFit" ]; then
    echo -e "${RED}Error: runFit executable not found!${NC}"
    exit 1
fi

./runFit

if [ $? -ne 0 ]; then
    echo -e "${RED}Fit failed!${NC}"
    exit 1
fi

# Check if output files exist
if [ ! -f "$OS_FILE" ] || [ ! -f "$LS_FILE" ]; then
    echo -e "${RED}Error: Output ROOT files not created!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Fit completed successfully${NC}"
echo -e "  Created: $OS_FILE, $LS_FILE, $WORKSPACE"
echo ""

# Step 3: Install Python dependencies (if needed)
echo -e "${YELLOW}[3/4] Checking Python dependencies...${NC}"
python3 -c "import numpy, matplotlib, seaborn, uproot" 2>/dev/null
if [ $? -ne 0 ]; then
    echo -e "${YELLOW}Installing required packages...${NC}"
    pip install -r requirements.txt
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to install dependencies!${NC}"
        exit 1
    fi
fi
echo -e "${GREEN}✓ Dependencies satisfied${NC}"
echo ""

# Step 4: Generate visualizations
echo -e "${YELLOW}[4/4] Generating visualizations...${NC}"

# Create plot directory
mkdir -p "$PLOT_DIR"

if [ "$USE_PYROOT" = true ]; then
    # Check if PyROOT is available
    python3 -c "import ROOT" 2>/dev/null
    if [ $? -eq 0 ]; then
        echo -e "Using ${GREEN}advanced${NC} visualization (with PyROOT)"
        python3 visualize_with_fits.py \
            --workspace "$WORKSPACE" \
            --os-file "$OS_FILE" \
            --ls-file "$LS_FILE" \
            --output-dir "$PLOT_DIR"
    else
        echo -e "${YELLOW}PyROOT not available, using basic visualization${NC}"
        python3 visualize_fits.py \
            --os-file "$OS_FILE" \
            --ls-file "$LS_FILE" \
            --output-dir "$PLOT_DIR"
    fi
else
    echo "Using basic visualization (without PyROOT)"
    python3 visualize_fits.py \
        --os-file "$OS_FILE" \
        --ls-file "$LS_FILE" \
        --output-dir "$PLOT_DIR"
fi

if [ $? -ne 0 ]; then
    echo -e "${RED}Visualization failed!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Visualizations created${NC}"
echo ""

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Analysis Complete!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "Output files:"
echo -e "  📁 ROOT files:     $OS_FILE, $LS_FILE"
echo -e "  📁 Workspace:      $WORKSPACE"
echo -e "  📁 Plots:          $PLOT_DIR/"
echo ""
echo -e "Generated plots:"
if [ -d "$PLOT_DIR" ]; then
    ls -1 "$PLOT_DIR" | sed 's/^/    /'
fi
echo ""

# Optional: Open plots
echo -e "${YELLOW}Open plots? (y/n)${NC}"
read -t 10 -n 1 response
echo ""
if [[ "$response" =~ ^[Yy]$ ]]; then
    if command -v xdg-open &> /dev/null; then
        xdg-open "$PLOT_DIR"/*.pdf 2>/dev/null &
    elif command -v open &> /dev/null; then
        open "$PLOT_DIR"/*.pdf 2>/dev/null &
    else
        echo "Cannot auto-open files. Please check $PLOT_DIR manually."
    fi
fi

echo -e "${GREEN}All done! Happy analyzing! 🎉${NC}"