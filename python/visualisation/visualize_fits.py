import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse
import seaborn as sns
import uproot
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
import warnings

class FitVisualizer:
    def __init__(self, workspace_file: str, verbose: bool = True):
        """
        Initialize visualizer
    
        Parameters:
            workspace_file: Path to RooFit workspace ROOT file
            verbose: Print loading information
        """
        self.workspace_file = Path(workspace_file)
        self.verbose = verbose
    
        # Data containers
        self.workspace_contents = {}
        self.has_pyroot = False
        self.workspace = None
        self.root_file = None
    
        # Check for PyROOT
        try:
            import ROOT
            self.has_pyroot = True
        except ImportError:
            pass
    
        # Load everything
        self._load_workspace()
    
        # Print summary
        if verbose:
            self._print_summary()

    def _print_summary(self):
        """Print a summary of what was loaded"""
        print("\n" + "=" * 60)
        print("FitVisualizer Summary")
        print("=" * 60)
        print(f"Workspace: {'✓ Loaded' if self.workspace else '✗ Not available'}")
        print(f"PyROOT: {'✓ Available' if self.has_pyroot else '✗ Not available (limited workspace access)'}")
        if not self.has_pyroot:
            print("\nTip: Install PyROOT for full RooFit workspace access:")
            print("     pip install root")
        print("=" * 60 + "\n")

    def _print(self, msg: str):
        """Print if verbose mode is on"""
        if self.verbose:
            print(msg)

    def _load_workspace(self):
        """Load workspace contents - try both uproot and PyROOT"""
        if not self.workspace_file.exists():
            self._print(f"Warning: Workspace file {self.workspace_file} not found")
            return
    
        # Try PyROOT first (best for RooFit workspaces)
        try:
            import ROOT
            self._load_workspace_pyroot()
            return
        except ImportError:
            if self.verbose:
                self._print("PyROOT not available, using uproot (limited workspace access)")
    
        # Fall back to uproot for basic file inspection
        try:
            with uproot.open(self.workspace_file) as f:
                # Just list what's in the file without trying to read everything
                keys = f.keys()
                self.workspace_contents = {'_keys': keys}
                self._print(f"✓ Workspace file opened ({len(keys)} objects)")
                if self.verbose and len(keys) > 0:
                    self._print("  Available objects:")
                    for key in list(keys)[:10]:
                        self._print(f"    - {key}")
                    if len(keys) > 10:
                        self._print(f"    ... and {len(keys) - 10} more")
        except Exception as e:
            # Silently handle workspace parsing issues - we mainly need the data
            if self.verbose:
                self._print(f"Note: Workspace parsing limited (this is normal with uproot)")
                self._print(f"      For full workspace access, install PyROOT")

    def _load_workspace_pyroot(self):
        """Load workspace using PyROOT (if available)"""
        import ROOT
    
        self.root_file = ROOT.TFile.Open(str(self.workspace_file))
        if not self.root_file or self.root_file.IsZombie():
            self._print(f"Warning: Could not open {self.workspace_file}")
            return
    
        # Try to get the workspace
        self.workspace = self.root_file.Get("w")
        if self.workspace:
            self._print(f"✓ Loaded RooWorkspace 'w' using PyROOT")
            self.workspace_contents['w'] = self.workspace
        
            # List available objects
            if self.verbose:
                self._print("  Available in workspace:")
                # Variables
                for var in self.workspace.allVars():
                    self._print(f"    Variable: {var.GetName()}")

                # PDFs
                for pdf in self.workspace.allPdfs():
                    self._print(f"    PDF: {pdf.GetName()}")

                # Datasets
                for data in self.workspace.allData():
                    self._print(f"    Dataset: {data.GetName()}")
        else:
            # List all objects in file
            keys = [key.GetName() for key in self.root_file.GetListOfKeys()]
            self._print(f"✓ Opened ROOT file with {len(keys)} objects")
            if self.verbose:
                for key in keys[:10]:
                    self._print(f"    - {key}")
    
        self.workspace_contents['root_file'] = self.root_file

    @staticmethod
    def _plot_with_errors(ax, data, bins, range=None, weights=None,
                          color='steelblue', alpha=0.6, label=None,
                          histtype='step', linewidth=1.5,
                          edgecolor='black', capsize=3):
        """
        Plot a histogram with statistical error bars.

        For unweighted data:  σ_i = sqrt(N_i)
        For weighted data:    σ_i = sqrt(Σ w²_i)   (standard HEP convention)

        Parameters
        ----------
        ax       : matplotlib Axes
        data     : array-like of values
        bins     : int or array-like (bin edges)
        range    : (low, high) or None
        weights  : array-like or None
        color    : fill/line colour
        alpha    : transparency for the filled area
        label    : legend label
        histtype : 'stepfilled' or 'step' – controls how the histogram body is drawn
        linewidth: line width for the histogram outline
        edgecolor: edge colour for stepfilled histtype
        capsize  : cap size for error bars

        Returns
        -------
        counts, bin_edges
        """
        data = np.asarray(data)
        if weights is not None:
            weights = np.asarray(weights)

        # Compute bin counts and squared-weight sum for errors
        counts, edges = np.histogram(data, bins=bins, range=range, weights=weights)
        if weights is not None:
            counts_w2, _ = np.histogram(data, bins=bins, range=range,
                                        weights=weights ** 2)
            errors = np.sqrt(counts_w2)
        else:
            errors = np.sqrt(counts)

        centres = 0.5 * (edges[:-1] + edges[1:])
        width = edges[1:] - edges[:-1]

        # Draw histogram body
        mask = counts > 0
        if histtype == 'stepfilled':
            ax.bar(centres, counts, width=width, color=color, alpha=alpha,
                   edgecolor=edgecolor, linewidth=linewidth, label=label)
            ax.errorbar(centres[mask], counts[mask], yerr=errors[mask],
                        fmt='none', color=edgecolor,
                        linewidth=linewidth * 0.8, capsize=capsize, alpha=min(alpha + 0.2, 1.0))
        elif histtype == 'step':
            ax.step(edges[:-1], counts, where='post', color=color,
                    linewidth=linewidth, alpha=alpha, label=label)
            ax.step([edges[-2], edges[-1]], [counts[-1], counts[-1]],
                    where='post', color=color, linewidth=linewidth, alpha=alpha)
            ax.errorbar(centres[mask], counts[mask], yerr=errors[mask],
                        fmt='none', color=color,
                        linewidth=linewidth * 0.8, capsize=capsize, alpha=min(alpha + 0.2, 1.0))
        else:  # 'points' — data points with error bars, no histogram body
            ax.errorbar(centres[mask], counts[mask], yerr=errors[mask],
                        fmt='o', color=color, markersize=4,
                        linewidth=linewidth, capsize=capsize,
                        alpha=alpha, label=label)

        return counts, edges

    def list_workspace_contents(self):
        """List all contents of the workspace in a formatted way"""
        if not self.workspace:
            print("No workspace loaded. Make sure PyROOT is installed.")
            return
    
        print("\nWorkspace contents:")
        print("-" * 60)
    
        # List variables
        print("Variables:")
        for var in self.workspace.allVars():
            print(f"  {var.GetName()}: {var.getVal()}")
    
        # List PDFs
        print("\nPDFs:")
        for pdf in self.workspace.allPdfs():
            print(f"  {pdf.GetName()}")
    
        # List datasets
        print("\nDatasets:")
        for data in self.workspace.allData():
            print(f"  {data.GetName()}: {data.numEntries()} entries")

    def extract_data_to_numpy(self, dataset_name: str, variables: list) -> dict:
        """
        Extract RooDataSet variables to numpy arrays
    
        Parameters:
        -----------
        dataset_name : str
        Name of the dataset in workspace ('osData' or 'lsData')
        variables : list of str
            List of variable names to extract (e.g., ['m1', 'm2', 'y1', 'y2', 'phi1', 'phi2'])
    
        Returns:
        --------
        dict : Dictionary with variable names as keys and numpy arrays as values
                Also includes 'weights' if the dataset is weighted
        """
        if not self.workspace:
            raise RuntimeError("No workspace loaded. Make sure PyROOT is installed.")
    
        dataset = self.workspace.data(dataset_name)
        if not dataset:
            raise ValueError(f"Dataset {dataset_name} not found in workspace")
    
        # Initialize dictionaries for data
        data_dict = {var: [] for var in variables}
        weights = []
    
        # Check if dataset is weighted
        is_weighted = dataset.isWeighted()
    
        self._print(f"Extracting {dataset.numEntries()} entries from {dataset_name}...")
    
        # Loop through entries
        for i in range(dataset.numEntries()):
            entry = dataset.get(i)
        
            # Extract each variable value
            for var_name in variables:
                var = self.workspace.var(var_name)
                if var:
                    data_dict[var_name].append(entry[var_name].getVal())
                else:
                    if i == 0:  # Only warn once
                        self._print(f"Warning: Variable {var_name} not found")
        
            # Extract weight if available
            if is_weighted:
                weights.append(dataset.weight())
    
        # Convert to numpy arrays
        for key in data_dict:
            data_dict[key] = np.array(data_dict[key])
    
        if is_weighted:
            data_dict['weights'] = np.array(weights)
            self._print(f"  Dataset is weighted (sum of weights: {np.sum(weights):.1f})")
        else:
            data_dict['weights'] = np.ones(len(data_dict[variables[0]]))
            self._print(f"  Dataset is unweighted")
    
        return data_dict
    
    def select_signal_region(self, nsigma):
        import ROOT
        # Define signal region around the mean
        m1_mean = self.workspace.var("m1_sig_mean").getVal()
        m1_sigma = self.workspace.var("m1_sig_sigma").getVal()
        m2_mean = self.workspace.var("m2_sig_mean").getVal()
        m2_sigma = self.workspace.var("m2_sig_sigma").getVal()

        # Create signal region cuts
        m1_cut = f"m1 > {m1_mean - nsigma*m1_sigma} && m1 < {m1_mean + nsigma*m1_sigma}"
        m2_cut = f"m2 > {m2_mean - nsigma*m2_sigma} && m2 < {m2_mean + nsigma*m2_sigma}"

        # Apply to your dataset
        osData_sigRegion = self.workspace.data("osData").reduce(f"{m1_cut} && {m2_cut}")
        lsData_sigRegion = self.workspace.data("lsData").reduce(f"{m1_cut} && {m2_cut}")

        getattr(self.workspace, 'import')(osData_sigRegion, ROOT.RooFit.Rename("osData_sigRegion"))
        getattr(self.workspace, 'import')(lsData_sigRegion, ROOT.RooFit.Rename("lsData_sigRegion"))


    def plot_mass_fits(self, dataset_name: str = 'osData', save_prefix: str = 'mass_fits'):
        """
        Plot the mass fits using RooFit's plotting capabilities
    
        Parameters:
        -----------
        dataset_name : str
            'osData' or 'lsData'
        save_prefix : str
            Prefix for saved plot files
        """
        import ROOT
    
        if not self.workspace:
            raise RuntimeError("No workspace loaded. Make sure PyROOT is installed.")
    
        dataset = self.workspace.data(dataset_name)
        if not dataset:
            self._print(f"Dataset {dataset_name} not found")
            return
    
        # Load the appropriate parameter snapshot
        snapshot_name = 'osParams' if 'os' in dataset_name.lower() else 'lsParams'
        if self.workspace.loadSnapshot(snapshot_name):
            self._print(f"Loaded snapshot: {snapshot_name}")
    
        # Get variables and model
        m1 = self.workspace.var("m1")
        m2 = self.workspace.var("m2")
        model = self.workspace.pdf("model")
    
        if not all([m1, m2, model]):
            self._print("Error: Could not find required variables/model in workspace")
            return
    
        # Create canvas
        c = ROOT.TCanvas("c", "Mass Fits", 1600, 800)
        c.Divide(2, 1)
    
        # Plot m1
        c.cd(1)
        frame1 = m1.frame(ROOT.RooFit.Title(f"{dataset_name}: m_{{1}} projection"))
        dataset.plotOn(frame1, ROOT.RooFit.Name("data"), ROOT.RooFit.DrawOption("pz"))
    
        # Plot model components
        model.plotOn(frame1, ROOT.RooFit.Name("model"))
        model.plotOn(frame1, ROOT.RooFit.Components("pdfBkgBkg"), 
                    ROOT.RooFit.LineStyle(ROOT.kDashed), 
                    ROOT.RooFit.LineColor(ROOT.kRed), 
                    ROOT.RooFit.Name("BkgBkg"))
        model.plotOn(frame1, ROOT.RooFit.Components("pdfSigSig"), 
                    ROOT.RooFit.LineColor(ROOT.kGreen+2), 
                    ROOT.RooFit.Name("SigSig"))
        frame1.Draw()
    
        # Plot m2
        c.cd(2)
        frame2 = m2.frame(ROOT.RooFit.Title(f"{dataset_name}: m_{{2}} projection"))
        dataset.plotOn(frame2, ROOT.RooFit.Name("data"), ROOT.RooFit.DrawOption("pz"))
    
        model.plotOn(frame2, ROOT.RooFit.Name("model"))
        model.plotOn(frame2, ROOT.RooFit.Components("pdfBkgBkg"), 
                    ROOT.RooFit.LineStyle(ROOT.kDashed), 
                    ROOT.RooFit.LineColor(ROOT.kRed))
        model.plotOn(frame2, ROOT.RooFit.Components("pdfSigSig"), 
                    ROOT.RooFit.LineColor(ROOT.kGreen+2))
        frame2.Draw()
    
        c.SaveAs(f"{save_prefix}_{dataset_name}.pdf")
        c.SaveAs(f"{save_prefix}_{dataset_name}.png")
        self._print(f"✓ Saved mass fit plots: {save_prefix}_{dataset_name}.pdf/png")
    
        return c

    def plot_pair_kinematics(self, region='', nsigma=3, compare: bool = True, normalize: bool = False, figsize: Tuple[float, float] = (14, 10)):
        """Plot pair kinematic distributions, optionally comparing OS and LS
    
        Args:
            region: 'signal' to use signal region data, otherwise use full data
            compare: If True, compare OS and LS distributions
            normalize: If True, normalize each distribution by its total number of events
            figsize: Figure size as (width, height)
        """

        if region == 'signal':
            os_data = self.extract_data_to_numpy(
                'osData_sigRegion', 
                ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
            )

            ls_data = self.extract_data_to_numpy(
                'lsData_sigRegion', 
                ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
            )
        else:
            os_data = self.extract_data_to_numpy(
                'osData', 
                ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
            )

            ls_data = self.extract_data_to_numpy(
                'lsData', 
                ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
            )

        is_weighted = not np.all(os_data['weights'] == 1.0)

        # Prepare weights based on normalization choice
        if normalize:
            n_os_events = len(os_data['weights'])
            n_ls_events = len(ls_data['weights']) if ls_data is not None else 0
        
            os_weights = os_data['weights'] / n_os_events
            ls_weights = ls_data['weights'] / n_ls_events if ls_data is not None else None
        
            # For unweighted case with normalization
            if not is_weighted:
                os_weights = np.ones(n_os_events) / n_os_events
                ls_weights = np.ones(n_ls_events) / n_ls_events if ls_data is not None else None
        else:
            # Use original weights (or no weights)
            os_weights = os_data['weights'] if is_weighted else None
            ls_weights = ls_data['weights'] if (is_weighted and ls_data is not None) else None

        # delta y
        os_data['deltaY'] = (np.abs(os_data['y1'] - os_data['y2']))
        ls_data['deltaY'] = (np.abs(ls_data['y1'] - ls_data['y2']))
        # delta phi
        os_data['deltaPhi'] = np.abs((os_data['phi1'] - os_data['phi2'] + np.pi) % (2 * np.pi) - np.pi)
        ls_data['deltaPhi'] = np.abs((ls_data['phi1'] - ls_data['phi2'] + np.pi) % (2 * np.pi) - np.pi)

        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Define common bin edges for fair comparison
        pt_bins = np.linspace(0.0, 20.0, 39)
        mass_bins = np.linspace(0.0, 10.0, 101)
        # Compute shared bin edges from the combined data range
        dy_min = min(os_data['deltaY'].min(), ls_data['deltaY'].min() if ls_data is not None else np.inf)
        dy_max = max(os_data['deltaY'].max(), ls_data['deltaY'].max() if ls_data is not None else -np.inf)
        deltay_bins = np.linspace(dy_min, dy_max, 21)  # 21 edges = 20 bins
        deltaphi_bins = 20
    
        # Set y-axis label based on normalization
        ylabel = 'Normalized Events' if normalize else 'Events'

        # Pair pT
        if os_data is not None:
            self._plot_with_errors(axes[0, 0], os_data['pairPt'], bins=pt_bins,
                                   weights=os_weights, histtype='step', linewidth=2,
                                   color='#2E86AB', alpha=0.9, label='OS')
        if compare and ls_data is not None:
            self._plot_with_errors(axes[0, 0], ls_data['pairPt'], bins=pt_bins,
                                   weights=ls_weights, histtype='step', linewidth=2,
                                   color='#F18F01', alpha=0.9, label='LS')
        axes[0, 0].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[0, 0].set_ylabel(ylabel)
        axes[0, 0].set_yscale('log')
        axes[0, 0].set_xlim(0.0, 20.0)
        axes[0, 0].set_title('Pair $p_T$', fontweight='bold')
        axes[0, 0].legend(frameon=True, shadow=True)
        axes[0, 0].grid(alpha=0.3)

        # Pair mass
        if os_data is not None:
            self._plot_with_errors(axes[0, 1], os_data['pairMass'], bins=mass_bins,
                                   weights=os_weights, histtype='step', linewidth=2,
                                   color='#C73E1D', alpha=0.9, label='OS')
        if compare and ls_data is not None:
            self._plot_with_errors(axes[0, 1], ls_data['pairMass'], bins=mass_bins,
                                   weights=ls_weights, histtype='step', linewidth=2,
                                   color='#6A4C93', alpha=0.9, label='LS')

        axes[0, 1].set_xlabel(r'$M_{pair}$ [GeV/$c^2$]')
        axes[0, 1].set_xlim(3.5, 7.0)
        axes[0, 1].set_ylabel(ylabel)
        axes[0, 1].set_title('Pair Mass', fontweight='bold')
        axes[0, 1].legend(frameon=True, shadow=True)
        axes[0, 1].grid(alpha=0.3)

        # deltaY
        if os_data is not None and 'deltaY' in os_data:
            self._plot_with_errors(axes[1, 0], os_data['deltaY'], bins=deltay_bins,
                           weights=os_weights, histtype='step', linewidth=2,
                           color='#A23B72', alpha=0.9, label='OS')
        if compare and ls_data is not None and 'deltaY' in ls_data:
            self._plot_with_errors(axes[1, 0], ls_data['deltaY'], bins=deltay_bins,
                           weights=ls_weights, histtype='step', linewidth=2,
                           color='#1982C4', alpha=0.9, label='LS')

        axes[1, 0].set_xlabel(r'$\Delta y$')
        axes[1, 0].set_ylabel(ylabel)
        axes[1, 0].set_title('$\Delta y$', fontweight='bold')
        axes[1, 0].legend(frameon=True, shadow=True)
        axes[1, 0].grid(alpha=0.3)

        # deltaPhi
        if os_data is not None and 'deltaPhi' in os_data:
            self._plot_with_errors(axes[1, 1], os_data['deltaPhi'], bins=deltaphi_bins,
                                   weights=os_weights, histtype='step', linewidth=2,
                                   color='#F18F01', alpha=0.9, label='OS')
        if compare and ls_data is not None and 'deltaPhi' in ls_data:
            self._plot_with_errors(axes[1, 1], ls_data['deltaPhi'], bins=deltaphi_bins,
                                   weights=ls_weights, histtype='step', linewidth=2,
                                   color='#8338EC', alpha=0.9, label='LS')

        axes[1, 1].set_xlabel(r'$\Delta\phi$')
        axes[1, 1].set_ylabel(ylabel)
        axes[1, 1].set_title('$\Delta\phi$', fontweight='bold')
        axes[1, 1].legend(frameon=True, shadow=True)
        axes[1, 1].grid(alpha=0.3)
    
        # Set title based on normalization
        title_suffix = ' (Normalized)' if normalize else ''
        region_suffix = f' in the {nsigma}σ region' if region=='signal' else ''
        plt.suptitle(f"{'OS vs LS Pair Kinematics Comparison' if compare else 'Pair Kinematics'}{region_suffix}{title_suffix}", 
                fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.show()
        return fig, axes
    
    def plot_kinematics(self, data_type: str = 'osData', region='', figsize: Tuple[float, float] = (14, 10)):
        """
        Plot kinematic distributions
        Args:
        - region: 'signal' to use signal region data, otherwise use full data
        """
        if region == 'signal':
            data_type = data_type + '_sigRegion'

        data = self.extract_data_to_numpy(
            data_type, 
            ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
        )
        
        fig, axes = plt.subplots(3, 3, figsize=figsize)

        is_weighted = not np.all(data['weights'] == 1.0)
        
        # pT1
        self._plot_with_errors(axes[0, 0], data['pt1'], bins=120, range=(0.0, 24.0),
                               weights=data['weights'] if is_weighted else None,
                               color='#2E86AB', alpha=0.6, edgecolor='black')

        axes[0, 0].set_xlabel(r'$p_{T,1}$ [GeV/$c$]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_yscale('log')
        axes[0, 0].set_xlim(0.0, 24.0)
        axes[0, 0].set_title(f'{data_type.upper()}: Particle 1 $p_T$', fontweight='bold')
        axes[0, 0].grid(alpha=0.3)
        
        # pT2
        self._plot_with_errors(axes[0, 1], data['pt2'], bins=120, range=(0.0, 24.0),
                               weights=data['weights'] if is_weighted else None,
                               color='#A23B72', alpha=0.6, edgecolor='black')

        axes[0, 1].set_xlabel(r'$p_{T,2}$ [GeV/$c$]')
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_xlim(0.0, 24.0)
        axes[0, 1].set_yscale('log')
        axes[0, 1].set_title(f'{data_type.upper()}: Particle 2 $p_T$', fontweight='bold')
        axes[0, 1].grid(alpha=0.3)
        
        # Pair pT
        self._plot_with_errors(axes[0, 2], data['pairPt'], bins=100, range=(0.0, 35.0),
                               weights=data['weights'] if is_weighted else None,
                               color='#F18F01', alpha=0.6, edgecolor='black')
        
        axes[0, 2].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[0, 2].set_ylabel('Events')
        axes[0, 2].set_xlim(0.0, 35.0)
        axes[0, 2].set_yscale('log')
        axes[0, 2].set_title(f'{data_type.upper()}: Pair $p_T$', fontweight='bold')
        axes[0, 2].grid(alpha=0.3)
        
        # y1
        self._plot_with_errors(axes[1, 0], data['y1'], bins=50,
                               weights=data['weights'] if is_weighted else None,
                               color='#2E86AB', alpha=0.6, edgecolor='black')

        axes[1, 0].set_xlabel(r'$y_{1}$')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_title(f'{data_type.upper()}: $y_{1}$', fontweight='bold')
        axes[1, 0].grid(alpha=0.3)

        # y2
        self._plot_with_errors(axes[1, 1], data['y2'], bins=50,
                               weights=data['weights'] if is_weighted else None,
                               color='#A23B72', alpha=0.6, edgecolor='black')
        axes[1, 1].set_xlabel(r'$y_{2}$')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title(f'{data_type.upper()}: $y_{2}$', fontweight='bold')
        axes[1, 1].grid(alpha=0.3)

        # delta y
        if 'deltaY' not in data:
            data['deltaY'] = (np.abs(data['y1'] - data['y2']))

        if is_weighted:
            self._plot_with_errors(axes[1, 2], data['deltaY'], bins=50,
                                   weights=data['weights'],
                                   color='#F18F01', alpha=0.6, edgecolor='black')
        else:
            self._plot_with_errors(axes[1, 2], data['deltaY'], bins=50,
                                   color='#F18F01', alpha=0.6, edgecolor='black')

        axes[1, 2].set_xlabel(r'$\Delta y$')
        axes[1, 2].set_ylabel('Events')
        axes[1, 2].set_title(f'{data_type.upper()}: $\Delta y$', fontweight='bold')
        axes[1, 2].grid(alpha=0.3)

        # phi1
        self._plot_with_errors(axes[2, 0], data['phi1'], bins=50,
                               weights=data['weights'] if is_weighted else None,
                               color='#2E86AB', alpha=0.6, edgecolor='black')

        axes[2, 0].set_xlabel(r'$\phi_{1}$')
        axes[2, 0].set_ylabel('Events')
        axes[2, 0].set_title(f'{data_type.upper()}: $\phi_{1}$', fontweight='bold')
        axes[2, 0].grid(alpha=0.3)

        # phi2
        self._plot_with_errors(axes[2, 1], data['phi2'], bins=50,
                               weights=data['weights'] if is_weighted else None,
                               color='#A23B72', alpha=0.6, edgecolor='black')
        axes[2, 1].set_xlabel(r'$\phi_{2}$')
        axes[2, 1].set_ylabel('Events')
        axes[2, 1].set_title(f'{data_type.upper()}: $\phi_{2}$', fontweight='bold')
        axes[2, 1].grid(alpha=0.3)

        # delta phi
        if 'deltaPhi' not in data:
            data['deltaPhi'] = np.abs((data['phi1'] - data['phi2'] + np.pi) % (2 * np.pi) - np.pi)

        if is_weighted:
            self._plot_with_errors(axes[2, 2], data['deltaPhi'], bins=50,
                                   weights=data['weights'],
                                   color='#F18F01', alpha=0.6, edgecolor='black')
        else:
            self._plot_with_errors(axes[2, 2], data['deltaPhi'], bins=50,
                                   color='#F18F01', alpha=0.6, edgecolor='black')
        
        axes[2, 2].set_xlabel(r'$\Delta \phi$')
        axes[2, 2].set_ylabel('Events')
        axes[2, 2].set_title(f'{data_type.upper()}: $\Delta \phi$', fontweight='bold')
        axes[2, 2].grid(alpha=0.3)
        
        plt.suptitle(f'{data_type.upper()} Kinematic Distributions', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.show()
        return fig, axes

    def plot_kinematic_distributions(self, dataset_name: str = 'osData', 
                                 save_prefix: str = 'kinematic_dist',
                                 bins: Optional[int] = None):
        """
        Plot kinematic distributions (y, φ, pt, pair variables) with adaptive binning
        and grouped layout.
        """

        self._print(f"Plotting kinematic distributions for {dataset_name}...")
        data = self.extract_data_to_numpy(
            dataset_name, 
            ['y1', 'y2', 'phi1', 'phi2', 'pt1', 'pt2', 'pairPt', 'pairMass']
        )

        is_weighted = not np.all(data['weights'] == 1.0)

        # Define variable metadata: name, label, range, and bin count
        variables = [
            ('y1', r'$y_1$', (-0.8, 0.8), 50),
            ('y2', r'$y_2$', (-0.8, 0.8), 50),
            ('phi1', r'$\phi_1$', (0, 2*np.pi), 50),
            ('phi2', r'$\phi_2$', (0, 2*np.pi), 50),
            ('pt1', r'$p_{T,1}$ [GeV/$c$]', (0, 24.0), 100),
            ('pt2', r'$p_{T,2}$ [GeV/$c$]', (0, 24.0), 100),
            ('pairPt', r'$p_{T}^{pair}$ [GeV/$c$]', (0, 50.0), 100),
            ('pairMass', r'$m_{pair}$ [GeV/$c^2$]', (3.0, 12.0), 18)
        ]

        # Override all bins if user specifies a single value
        if bins is not None:
            variables = [(v, l, r, bins) for (v, l, r, _) in variables]

        # Split variables into groups for visual clarity
        single_vars = variables[:6]  # y1, y2, phi1, phi2, pt1, pt2
        pair_vars = variables[6:]    # pairPt, pairMass

        # Count total rows
        n_rows = (2 if is_weighted else 1) + 2  # single-particle, pair-level (+weight plots)
        fig = plt.figure(figsize=(14, 4 * n_rows))
        gs = GridSpec(n_rows, 3, figure=fig, hspace=0.5, wspace=0.3)

        # ---------------------------
        # Single-particle variables
        # ---------------------------
        for i, (var_name, label, range_vals, nbins) in enumerate(single_vars):
            row, col = divmod(i, 3)
            ax = fig.add_subplot(gs[row, col])
            self._plot_with_errors(ax, data[var_name], bins=nbins, range=range_vals,
                                   weights=data['weights'] if is_weighted else None,
                                   color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.5)
            ax.set_xlabel(label, fontsize=12)
            ax.set_ylabel('Entries', fontsize=11)
            ax.set_title(f'{label}', fontsize=12, pad=8)
            ax.grid(True, alpha=0.3)

        # ---------------------------
        # Pair-level variables
        # ---------------------------
        for i, (var_name, label, range_vals, nbins) in enumerate(pair_vars):
            ax = fig.add_subplot(gs[2, i])
            self._plot_with_errors(ax, data[var_name], bins=nbins, range=range_vals,
                                   weights=data['weights'] if is_weighted else None,
                                   color='darkorange', alpha=0.8, edgecolor='black', linewidth=0.5)
            ax.set_xlabel(label, fontsize=12)
            ax.set_ylabel('Entries', fontsize=11)
            ax.set_title(f'{label}', fontsize=12, pad=8)
            ax.grid(True, alpha=0.3)

        # ---------------------------
        # Weighted plots (optional)
        # ---------------------------
        if is_weighted:
            # Weight distribution
            ax = fig.add_subplot(gs[-2, 0])
            ax.hist(data['weights'], bins=100, alpha=0.7, color='forestgreen',
                    edgecolor='black', linewidth=0.5)
            ax.set_xlabel('Weight', fontsize=12)
            ax.set_ylabel('Entries', fontsize=11)
            ax.set_title('Weight Distribution', fontsize=12)
            ax.set_yscale('log')
            ax.grid(True, alpha=0.3)

            # Weighted pair mass (for example)
            ax = fig.add_subplot(gs[-2, 1])
            self._plot_with_errors(ax, data['pairMass'], bins=60, range=(3.0, 12.0),
                                   weights=data['weights'],
                                   color='crimson', alpha=0.8, edgecolor='black', linewidth=0.5,
                                   histtype='points')
            ax.set_xlabel(r'$m_{pair}$ [GeV/$c^2$]', fontsize=12)
            ax.set_ylabel('Weighted entries', fontsize=11)
            ax.set_title('Weighted pair mass', fontsize=12)
            ax.grid(True, alpha=0.3)

            # Weight stats
            ax = fig.add_subplot(gs[-2, 2])
            stats_text = (
                f"Weight stats\n"
                f"Mean: {np.mean(data['weights']):.3f}\n"
                f"Std:  {np.std(data['weights']):.3f}\n"
                f"Min:  {np.min(data['weights']):.3f}\n"
                f"Max:  {np.max(data['weights']):.3f}\n"
                f"Σw:   {np.sum(data['weights']):.1f}"
            )
            ax.text(0.05, 0.5, stats_text, transform=ax.transAxes,
                    fontsize=11, family='monospace',
                    verticalalignment='center',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.4))
            ax.axis('off')

        plt.suptitle(f"Kinematic Distributions – {dataset_name}",
                     fontsize=16, fontweight='bold', y=0.98)

        plt.savefig(f"{save_prefix}_{dataset_name}.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{save_prefix}_{dataset_name}.png", dpi=300, bbox_inches='tight')
        self._print(f"✓ Saved kinematic plots: {save_prefix}_{dataset_name}.pdf/png")

        return fig


    def plot_2d_correlations(self, dataset_name: str = 'osData', 
                        save_prefix: str = 'correlations',
                        use_weights: bool = True, bins: int = 50):
        """
        Plot 2D correlations between variables
    
        Parameters:
        -----------
        dataset_name : str
            'osData' or 'lsData'
        save_prefix : str
            Prefix for saved plot files
        use_weights : bool
            Whether to use weights for the 2D histogram
        bins : int
            Number of bins for 2D histograms
        """
        # Extract data
        self._print(f"Plotting 2D correlations for {dataset_name} ({'weighted' if use_weights else 'unweighted'})...")
        data = self.extract_data_to_numpy(dataset_name, 
                                        ['m1', 'm2', 'y1', 'y2', 'phi1', 'phi2'])
    
        weights = data['weights'] if use_weights else None
    
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
        correlations = [
            ('m1', 'm2', r'$m_1$ [GeV/c$^2$]', r'$m_2$ [GeV/c$^2$]'),
            ('y1', 'y2', r'$y_1$', r'$y_2$'),
            ('phi1', 'phi2', r'$\phi_1$', r'$\phi_2$'),
            ('m1', 'y1', r'$m_1$ [GeV/c$^2$]', r'$y_1$'),
            ('m2', 'y2', r'$m_2$ [GeV/c$^2$]', r'$y_2$'),
            ('y1', 'phi1', r'$y_1$', r'$\phi_1$'),
        ]
    
        for ax, (var1, var2, label1, label2) in zip(axes.flat, correlations):
            h = ax.hist2d(data[var1], data[var2], bins=bins, 
                        weights=weights, cmap='viridis')
            ax.set_xlabel(label1, fontsize=12)
            ax.set_ylabel(label2, fontsize=12)
            ax.set_title(f'{var1} vs {var2}', fontsize=11)
            plt.colorbar(h[3], ax=ax, label='Entries')
    
        weight_str = 'Weighted' if use_weights else 'Unweighted'
        plt.suptitle(f'2D Correlations: {dataset_name} ({weight_str})', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
    
        plt.savefig(f"{save_prefix}_{dataset_name}_{'weighted' if use_weights else 'unweighted'}.pdf", 
                dpi=300, bbox_inches='tight')
        plt.savefig(f"{save_prefix}_{dataset_name}_{'weighted' if use_weights else 'unweighted'}.png", 
                dpi=300, bbox_inches='tight')
        self._print(f"✓ Saved correlation plots: {save_prefix}_{dataset_name}_"
            f"{'weighted' if use_weights else 'unweighted'}.pdf/png")
    
        return fig

    def plot_all(self, dataset_name: str = 'osData'):
        """
        Generate all plots for a given dataset
        """
        print(f"\n{'=' * 60}")
        print(f"Generating all plots for {dataset_name}")
        print(f"{'=' * 60}\n")
    
        # Mass fits
        self.plot_mass_fits(dataset_name, save_prefix='mass_fits')
    
        # Kinematic distributions
        self.plot_kinematic_distributions(dataset_name, save_prefix='kinematic_dist')
    
        # 2D correlations (weighted and unweighted)
        self.plot_2d_correlations(dataset_name, save_prefix='correlations', 
                             use_weights=False)
        self.plot_2d_correlations(dataset_name, save_prefix='correlations', 
                             use_weights=True)
    
        print(f"\n✓ All plots generated for {dataset_name}!\n")