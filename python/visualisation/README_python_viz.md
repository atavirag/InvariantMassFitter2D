# Python Visualization for 2D Invariant Mass Fits

Beautiful, publication-quality visualizations for your D meson pair analysis using Python.

## Why Python for Visualization?

✅ **Beautiful plots** - matplotlib/seaborn produce cleaner, more customizable plots  
✅ **Faster iteration** - No recompilation needed  
✅ **Better aesthetics** - Modern color schemes, fonts, and styling  
✅ **Easy customization** - Simple to adjust colors, layouts, labels  
✅ **Interactive options** - Can use plotly for interactive plots  
✅ **Publication-ready** - Direct output to high-quality PDF/PNG/SVG  

## Workflow

```
C++ (Heavy lifting)          Python (Visualization)
├─ RooFit fitting       →    ├─ Load ROOT files with uproot
├─ TreeSelector         →    ├─ Extract fit parameters (PyROOT)
├─ Save workspace       →    ├─ Create beautiful plots
└─ Export trees         →    └─ Generate publication figures
```

## Installation

### Option 1: Without PyROOT (Basic, uses uproot)

```bash
pip install -r requirements.txt
```

### Option 2: With PyROOT (Advanced, includes fit overlay)

First install ROOT, then:
```bash
pip install -r requirements.txt
# PyROOT comes with ROOT installation
```

## Usage

### Basic Visualization (No fit parameters needed)

```bash
# Generate all plots from data files only
python visualize_fits.py \
    --os-file os_pairs.root \
    --ls-file ls_pairs.root \
    --output-dir plots/
```

This creates:
- `2d_mass_distributions.pdf` - Side-by-side OS/LS 2D plots
- `1d_projections.pdf` - All 1D mass projections
- `correlation_plot.pdf` - Correlation between m1 and m2
- `os_vs_ls_comparison.pdf` - Overlay comparison
- `kinematic_distributions.pdf` - pT distributions
- `summary.pdf` - Comprehensive summary page

### Advanced Visualization (With fit parameters)

```bash
# Includes fit function overlay and yield extraction
python visualize_with_fits.py \
    --workspace fit_workspace.root \
    --os-file os_pairs.root \
    --ls-file ls_pairs.root \
    --output-dir plots_advanced/
```

Additional plots:
- `os_fit_overlay.pdf` - Data with fit curves
- `yield_summary.pdf` - Bar and pie charts of yields
- `significance.pdf` - Signal quality metrics
- `publication_plot.pdf` - Journal-ready figure

## Example Output

### 2D Mass Distribution
High-resolution 2D histograms with proper colormaps, statistics boxes, and reference lines.

### Fit Overlay
Clean 1D projections with Gaussian signal overlays, shaded regions, and uncertainty bands.

### Yield Summary
Professional bar charts and pie charts showing fit component contributions.

## Customization

### Change Color Scheme

Edit the script:
```python
# In visualize_fits.py or visualize_with_fits.py
plt.style.use('seaborn-v0_8-paper')  # Try: 'ggplot', 'bmh', 'classic'
sns.set_palette("husl")              # Try: "deep", "muted", "bright"
```

### Use Science Plot Style

```python
import scienceplots
plt.style.use(['science', 'ieee'])  # IEEE style
# or
plt.style.use(['science', 'nature'])  # Nature style
```

### Adjust Figure Sizes

```python
fig = plt.figure(figsize=(14, 10))  # Width, height in inches
```

### Export to Different Formats

```python
plt.savefig(output, dpi=300, format='png')   # PNG
plt.savefig(output, dpi=300, format='pdf')   # PDF (vector)
plt.savefig(output, dpi=300, format='svg')   # SVG (vector)
```

## Interactive Plots with Plotly

Create an interactive version:

```python
import plotly.graph_objects as go
import plotly.express as px

# 2D histogram
fig = go.Figure(data=go.Histogram2d(x=m1_vals, y=m2_vals, 
                                    colorscale='Viridis'))
fig.update_layout(title='Interactive 2D Mass Distribution',
                  xaxis_title='m1 [GeV/c²]',
                  yaxis_title='m2 [GeV/c²]')
fig.show()
fig.write_html('interactive_2d.html')  # Save as HTML
```

## Integration with C++ Pipeline

### Complete Workflow

```bash
# Step 1: Run C++ fitting
./runFit
# Creates: os_pairs.root, ls_pairs.root, fit_workspace.root

# Step 2: Generate Python visualizations
python visualize_with_fits.py \
    --workspace fit_workspace.root \
    --os-file os_pairs.root \
    --ls-file ls_pairs.root

# Step 3: All plots are in plots_advanced/
```

### Automated Script

Create `run_analysis.sh`:

```bash
#!/bin/bash

echo "Running C++ fit..."
./runFit

if [ $? -eq 0 ]; then
    echo "Fit successful! Generating visualizations..."
    python visualize_with_fits.py \
        --workspace fit_workspace.root \
        --os-file os_pairs.root \
        --ls-file ls_pairs.root \
        --output-dir plots/
    
    echo "Analysis complete! Check plots/ directory"
else
    echo "Fit failed!"
    exit 1
fi
```

## Advanced Features

### 1. Reading RooWorkspace with PyROOT

If you need to extract specific fit parameters:

```python
import ROOT

# Open workspace
f = ROOT.TFile.Open("fit_workspace.root")
ws = f.Get("w")

# Get fit result
fitResult = ws.obj("osFitResult")
if fitResult:
    fitResult.Print("v")

# Get parameters
nSigSig = ws.var("nSigSig")
print(f"Signal yield: {nSigSig.getVal()} ± {nSigSig.getError()}")

# Get PDF and data
model = ws.pdf("model")
data = ws.data("osData")

# Create RooPlot and convert to matplotlib
frame = ws.var("m1").frame()
data.plotOn(frame)
model.plotOn(frame)

# Extract data from RooPlot for matplotlib
# (This requires some conversion, see visualize_with_fits.py)
```

### 2. Systematic Uncertainties Visualization

Add error bands:

```python
# Calculate systematic variations
signal_central = calculate_signal(params)
signal_up = calculate_signal(params_varied_up)
signal_down = calculate_signal(params_varied_down)

# Plot with uncertainty band
ax.plot(x, signal_central, 'r-', label='Central')
ax.fill_between(x, signal_down, signal_up, 
                alpha=0.3, color='red', label='Systematic')
```

### 3. Multiple Fit Comparison

Compare different fit scenarios:

```python
# Load multiple workspaces
fits = {
    'Standard': 'fit_standard.root',
    'Alternative': 'fit_alternative.root',
    'Background Only': 'fit_bkg_only.root'
}

fig, ax = plt.subplots(figsize=(12, 6))

for name, ws_file in fits.items():
    # Extract and plot each fit
    yields = extract_yields(ws_file)
    ax.bar(name, yields['nSigSig'], label=name)

ax.legend()
ax.set_ylabel('Signal Yield')
ax.set_title('Comparison of Fit Scenarios')
```

### 4. Contour Plots for Correlations

Show parameter correlations:

```python
from matplotlib.patches import Ellipse
from scipy.stats import chi2

# Get covariance matrix from fit result
cov = fitResult.covarianceMatrix()

# Create confidence ellipses
def confidence_ellipse(mean, cov, ax, n_std=2.0, **kwargs):
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    
    ellipse = Ellipse(mean, width=ell_radius_x * scale_x * 2,
                     height=ell_radius_y * scale_y * 2, **kwargs)
    ax.add_patch(ellipse)

# Plot
fig, ax = plt.subplots()
confidence_ellipse([mean1, mean2], cov_matrix, ax, 
                   n_std=1, facecolor='none', edgecolor='red')
confidence_ellipse([mean1, mean2], cov_matrix, ax,
                   n_std=2, facecolor='none', edgecolor='blue')
```

## Tips for Publication-Quality Plots

### 1. Font Sizes for Different Media

```python
# For presentations
plt.rcParams.update({'font.size': 14, 'axes.labelsize': 16})

# For posters
plt.rcParams.update({'font.size': 18, 'axes.labelsize': 20})

# For papers (default)
plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12})
```

### 2. High DPI for Publications

```python
# Nature, Science: 300-600 dpi
plt.savefig('figure.pdf', dpi=600, bbox_inches='tight')

# PRD, PRL: 300 dpi minimum
plt.savefig('figure.eps', dpi=300, format='eps')
```

### 3. Consistent Styling

```python
# Define color palette once
COLORS = {
    'os': '#2E86AB',
    'ls': '#F18F01',
    'signal': '#A23B72',
    'background': '#C73E1D'
}

# Use throughout
ax.plot(x, y_os, color=COLORS['os'], label='OS')
ax.plot(x, y_ls, color=COLORS['ls'], label='LS')
```

### 4. LaTeX Rendering (if available)

```python
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Now you can use:
ax.set_xlabel(r'$m_{D^0}$ [GeV/$c^2$]')
ax.set_title(r'$D^0 \bar{D}^0$ Invariant Mass')
```

## Troubleshooting

### "uproot not found"
```bash
pip install uproot awkward
```

### "PyROOT not available"
PyROOT comes with ROOT. Make sure ROOT is properly sourced:
```bash
source /path/to/root/bin/thisroot.sh
python -c "import ROOT; print(ROOT.__version__)"
```

### "Cannot read workspace"
Make sure the workspace was saved correctly:
```bash
root -l fit_workspace.root
# In ROOT:
.ls
w->Print()
```

### Plots look pixelated
Increase DPI:
```python
plt.rcParams['figure.dpi'] = 150  # For screen
plt.savefig('plot.pdf', dpi=300)  # For print
```

### Legend overlaps with data
Adjust position:
```python
ax.legend(loc='upper right')        # Corner
ax.legend(bbox_to_anchor=(1.05, 1)) # Outside plot
ax.legend(loc='best')                # Automatic
```

## Performance Tips

### For Large Datasets

```python
# Use 2D histogram instead of scatter
h, xedges, yedges = np.histogram2d(x, y, bins=100)
ax.imshow(h.T, ...)  # Much faster than scatter

# Or use hexbin
ax.hexbin(x, y, gridsize=50, cmap='viridis')
```

### Parallel Processing

```python
from multiprocessing import Pool

def plot_slice(data_slice):
    # Create plot for this slice
    ...

with Pool(4) as p:
    p.map(plot_slice, data_slices)
```

## Comparison: ROOT vs Python Plotting

| Feature | ROOT | Python |
|---------|------|--------|
| Speed (large data) | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| Aesthetics | ⭐⭐ | ⭐⭐⭐⭐⭐ |
| Customization | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| Learning curve | ⭐⭐ | ⭐⭐⭐⭐ |
| Interactive plots | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| Publication ready | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| Integration with RooFit | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |

## Example Gallery

Check `examples/` directory for:
- `example_2d.py` - 2D histogram variations
- `example_fit_overlay.py` - Fit function overlays
- `example_correlation.py` - Correlation plots
- `example_significance.py` - Significance calculations
- `example_systematic.py` - Systematic uncertainties

## Resources

- **matplotlib gallery**: https://matplotlib.org/stable/gallery/
- **seaborn tutorial**: https://seaborn.pydata.org/tutorial.html
- **uproot documentation**: https://uproot.readthedocs.io/
- **PyROOT manual**: https://root.cern/manual/python/
- **SciencePlots**: https://github.com/garrettj403/SciencePlots

## Contributing

Suggestions for new plot types or improvements? Feel free to extend the scripts!

## License

Free to use and modify for your analysis.