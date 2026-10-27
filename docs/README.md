# D Meson Pair Analysis

Complete C++/Python framework for 2D invariant mass fitting and visualization of D meson pairs.

## Project Structure

```
analysis/
├── src/              # C++ source code (version controlled)
├── build/            # Build artifacts (gitignored)
├── data/             # ROOT data files (gitignored)
├── plots/            # Generated plots (gitignored)
├── python/           # Python scripts (version controlled)
├── scripts/          # Automation scripts (version controlled)
├── docs/             # Documentation (version controlled)
├── Makefile          # Build system
├── requirements.txt  # Python dependencies
└── .gitignore        # Git ignore rules
```

## Analysis steps

1. Prefiltering
```
source src/01-prefiltering/run_preselect.sh
```
2. 2D fitting, with and without efficiency corrections
```
source src/efficiencies/run_efficiency_map.sh
source src/01-fitting/run_fitter.sh
```
3. Signal extraction with s-plot
```
source src/03-splot/run_splot.sh
```
4. Cross section extraction
```
python python.cross_sections.compute_xsec.py
```
