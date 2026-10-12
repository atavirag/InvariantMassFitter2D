#!/usr/bin/env python3
"""
visualize_fits.py - Visualization functions for fits

This module provides a clean interface for creating publication-quality plots
from RooFit workspaces. Designed for use in Jupyter notebooks.

Usage in Jupyter:
    from visualize_fits import FitVisualizer
    
    viz = FitVisualizer(workspace_file='../../data/fit_workspace.root',
                            os_file='../../data/os_pairs.root',
                            ls_file='../../data/ls_pairs.root')
    
    # Create individual plots
    fig, ax = viz.plot_2d_os(bins=50)
    plt.show()
    
    # Or use convenience methods
    viz.show_correlation()
    viz.show_projections()
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import seaborn as sns
import uproot
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
import warnings

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['font.family'] = 'serif'


class FitVisualizer:
    """
    Visualize 2D invariant mass and kinematic distributions from the pair tree
    
    This class provides methods to create various plots from ROOT files
    containing RooFit workspaces and data. All methods return matplotlib
    figure and axes objects for easy customization in notebooks.
    
    Attributes:
        workspace_file (Path): Path to RooFit workspace ROOT file
        os_file (Path): Path to OS data ROOT file
        ls_file (Path): Path to LS data ROOT file
        os_data (dict): Dictionary containing OS data arrays
        ls_data (dict): Dictionary containing LS data arrays
        workspace_contents (dict): Contents of the workspace
    """
    
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
        print(f"Workspace: {'✓ Loaded' if self.workspace_contents else '✗ Not available'}")
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
        
        f = ROOT.TFile.Open(str(self.workspace_file))
        if not f or f.IsZombie():
            self._print(f"Warning: Could not open {self.workspace_file}")
            return
        
        # Try to get the workspace
        ws = f.Get("w")
        if ws:
            self._print(f"✓ Loaded RooWorkspace 'w' using PyROOT")
            self.workspace_contents['w'] = ws
            
            # List available objects
            if self.verbose:
                self._print("  Available in workspace:")
                for var in ws.allVars():
                    self._print(f"    Variable: {var.GetName()}")
                for pdf in ws.allPdfs():
                    self._print(f"    PDF: {pdf.GetName()}")
                for data in ws.allData():
                    self._print(f"    Dataset: {data.GetName()}")
        else:
            # List all objects in file
            keys = [key.GetName() for key in f.GetListOfKeys()]
            self._print(f"✓ Opened ROOT file with {len(keys)} objects")
            if self.verbose:
                for key in keys[:10]:
                    self._print(f"    - {key}")
        
        self.workspace_contents['root_file'] = f
    
    def list_workspace_contents(self):
        """Print all objects in the workspace"""
        if not self.workspace_contents:
            print("No workspace contents loaded")
            return
        
        print(f"\nWorkspace contents:")
        print("-" * 60)
        
        # If we have PyROOT workspace
        if 'w' in self.workspace_contents:
            ws = self.workspace_contents['w']
            print("\nVariables:")
            for var in ws.allVars():
                print(f"  {var.GetName()}: {var.getVal()}")
            print("\nPDFs:")
            for pdf in ws.allPdfs():
                print(f"  {pdf.GetName()}")
            print("\nDatasets:")
            for data in ws.allData():
                print(f"  {data.GetName()}: {data.numEntries()} entries")
        # If we just have keys from uproot
        elif '_keys' in self.workspace_contents:
            for key in self.workspace_contents['_keys']:
                print(f"  {key}")
        else:
            for key in self.workspace_contents.keys():
                print(f"  {key}")
    
    def get_workspace_object(self, name: str):
        """
        Get an object from the workspace by name
        
        Requires PyROOT. Returns None if not available.
        """
        if 'workspace' not in self.workspace_contents:
            print("Workspace not loaded with PyROOT. Cannot access objects.")
            print("Install PyROOT for full workspace access.")
            return None
        
        ws = self.workspace_contents['workspace']
        
        # Try different object types
        obj = ws.var(name)
        if obj: return obj
        
        obj = ws.pdf(name)
        if obj: return obj
        
        obj = ws.data(name)
        if obj: return obj
        
        obj = ws.function(name)
        if obj: return obj
        
        print(f"Object '{name}' not found in workspace")
        return None
    
    def plot_roofit_frame(self, frame_name: str = None, figsize: Tuple[float, float] = (10, 7)):
        """
        Plot a RooFit frame from the workspace
        
        This requires PyROOT and that the workspace contains saved frames.
        """
        if 'root_file' not in self.workspace_contents:
            print("ROOT file not loaded. Cannot access RooPlot frames.")
            return None, None
        
        import ROOT
        f = self.workspace_contents['root_file']
        
        # If no frame name given, try to find one
        if frame_name is None:
            keys = [key.GetName() for key in f.GetListOfKeys()]
            frame_keys = [k for k in keys if 'frame' in k.lower() or 'plot' in k.lower()]
            if frame_keys:
                frame_name = frame_keys[0]
                print(f"Using frame: {frame_name}")
            else:
                print("No frame specified and none found automatically")
                print(f"Available objects: {keys}")
                return None, None
        
        # Get the frame
        frame = f.Get(frame_name)
        if not frame:
            print(f"Frame '{frame_name}' not found")
            return None, None
        
        # Convert RooPlot to matplotlib
        # This is tricky - we'll extract the histogram data
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get histograms from frame
        n_items = frame.numItems()
        for i in range(n_items):
            obj = frame.getObject(i)
            if obj:
                name = obj.GetName()
                # Try to extract data and plot
                # This is a simplified version - may need adjustment
                ax.set_title(f"RooFit Frame: {frame_name}", fontweight='bold')
                ax.set_xlabel(frame.GetXaxis().GetTitle())
                ax.set_ylabel(frame.GetYaxis().GetTitle())

        plt.tight_layout()
        return fig, ax
    
    # ========== 2D Distribution Plots ==========
    
    def plot_2d_os(self, bins: int = 50, cmap: str = 'viridis', 
                   figsize: Tuple[float, float] = (8, 7)) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create 2D OS mass distribution plot
        
        Parameters:
            bins: Number of bins for histogram
            cmap: Colormap name
            figsize: Figure size (width, height)
            
        Returns:
            fig, ax: Matplotlib figure and axes
        """
        if self.os_data is None:
            raise ValueError("No OS data available")
        
        m1_os, m2_os = self.get_os_masses()
        
        fig, ax = plt.subplots(figsize=figsize)
        
        h, xedges, yedges = np.histogram2d(m1_os, m2_os, bins=bins)
        
        im = ax.imshow(h.T, origin='lower', aspect='auto', 
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                      cmap=cmap, interpolation='nearest')
        
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_title('Opposite Sign (OS)', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Events', rotation=270, labelpad=20)
        
        # Diagonal reference line
        ax.plot([xedges[0], xedges[-1]], [yedges[0], yedges[-1]], 
                'r--', alpha=0.5, linewidth=1.5, label='$m_1 = m_2$')
        ax.legend(loc='upper left')
        
        # Statistics box
        n_events = len(m1_os)
        textstr = f'Events: {n_events:,}\nMean $m_1$: {np.mean(m1_os):.3f}\nMean $m_2$: {np.mean(m2_os):.3f}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        return fig, ax
    
    def plot_2d_ls(self, bins: int = 50, cmap: str = 'plasma',
                   figsize: Tuple[float, float] = (8, 7)) -> Tuple[plt.Figure, plt.Axes]:
        """Create 2D LS mass distribution plot"""
        if self.ls_data is None:
            raise ValueError("No LS data available")
        
        m1_ls, m2_ls = self.get_ls_masses()
        
        fig, ax = plt.subplots(figsize=figsize)
        
        h, xedges, yedges = np.histogram2d(m1_ls, m2_ls, bins=bins)
        
        im = ax.imshow(h.T, origin='lower', aspect='auto',
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                      cmap=cmap, interpolation='nearest')
        
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_title('Like Sign (LS)', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Events', rotation=270, labelpad=20)
        
        ax.plot([xedges[0], xedges[-1]], [yedges[0], yedges[-1]], 
                'r--', alpha=0.5, linewidth=1.5, label='$m_1 = m_2$')
        ax.legend(loc='upper left')
        
        n_events = len(m1_ls)
        textstr = f'Events: {n_events:,}\nMean $m_1$: {np.mean(m1_ls):.3f}\nMean $m_2$: {np.mean(m2_ls):.3f}'
        props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        return fig, ax
    
    def plot_2d_both(self, bins: int = 50, figsize: Tuple[float, float] = (14, 6)):
        """Create side-by-side 2D OS and LS plots"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # OS
        if self.os_data is not None:
            m1_os, m2_os = self.get_os_masses()
            h, xe, ye = np.histogram2d(m1_os, m2_os, bins=bins)
            im1 = ax1.imshow(h.T, origin='lower', aspect='auto',
                            extent=[xe[0], xe[-1], ye[0], ye[-1]],
                            cmap='viridis', interpolation='nearest')
            ax1.set_xlabel(r'$m_1$ [GeV/$c^2$]')
            ax1.set_ylabel(r'$m_2$ [GeV/$c^2$]')
            ax1.set_title('Opposite Sign (OS)', fontweight='bold')
            plt.colorbar(im1, ax=ax1, label='Events')
            ax1.plot([xe[0], xe[-1]], [ye[0], ye[-1]], 'r--', alpha=0.5)
        
        # LS
        if self.ls_data is not None:
            m1_ls, m2_ls = self.get_ls_masses()
            h, xe, ye = np.histogram2d(m1_ls, m2_ls, bins=bins)
            im2 = ax2.imshow(h.T, origin='lower', aspect='auto',
                            extent=[xe[0], xe[-1], ye[0], ye[-1]],
                            cmap='plasma', interpolation='nearest')
            ax2.set_xlabel(r'$m_1$ [GeV/$c^2$]')
            ax2.set_ylabel(r'$m_2$ [GeV/$c^2$]')
            ax2.set_title('Like Sign (LS)', fontweight='bold')
            plt.colorbar(im2, ax=ax2, label='Events')
            ax2.plot([xe[0], xe[-1]], [ye[0], ye[-1]], 'r--', alpha=0.5)
        
        plt.suptitle('2D Invariant Mass Distributions', fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig, (ax1, ax2)
    
    # ========== 1D Projection Plots ==========
    
    def plot_m1_projection(self, bins: int = 60, compare: bool = True,
                          figsize: Tuple[float, float] = (10, 6)) -> Tuple[plt.Figure, plt.Axes]:
        """Plot m1 projection, optionally comparing OS and LS"""
        fig, ax = plt.subplots(figsize=figsize)
        
        if self.os_data is not None:
            m1_os, _ = self.get_os_masses()
            ax.hist(m1_os, bins=bins, histtype='step', linewidth=2, 
                   color='#2E86AB', label=f'OS ({len(m1_os)} events)', alpha=0.9)
        
        if compare and self.ls_data is not None:
            m1_ls, _ = self.get_ls_masses()
            ax.hist(m1_ls, bins=bins, histtype='step', linewidth=2, 
                   color='#F18F01', label=f'LS ({len(m1_ls)} events)', alpha=0.9)
        
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel('Events / bin', fontsize=13)
        ax.set_title('$m_1$ Projection', fontsize=14, fontweight='bold')
        ax.legend(frameon=True, shadow=True)
        ax.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_m2_projection(self, bins: int = 60, compare: bool = True,
                          figsize: Tuple[float, float] = (10, 6)) -> Tuple[plt.Figure, plt.Axes]:
        """Plot m2 projection, optionally comparing OS and LS"""
        fig, ax = plt.subplots(figsize=figsize)
        
        if self.os_data is not None:
            _, m2_os = self.get_os_masses()
            ax.hist(m2_os, bins=bins, histtype='step', linewidth=2,
                   color='#A23B72', label=f'OS ({len(m2_os)} events)', alpha=0.9)
        
        if compare and self.ls_data is not None:
            _, m2_ls = self.get_ls_masses()
            ax.hist(m2_ls, bins=bins, histtype='step', linewidth=2,
                   color='#C73E1D', label=f'LS ({len(m2_ls)} events)', alpha=0.9)
        
        ax.set_xlabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel('Events / bin', fontsize=13)
        ax.set_title('$m_2$ Projection', fontsize=14, fontweight='bold')
        ax.legend(frameon=True, shadow=True)
        ax.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_projections(self, bins: int = 60, figsize: Tuple[float, float] = (14, 6)):
        """Plot both m1 and m2 projections side by side"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # m1
        if self.os_data is not None:
            m1_os, _ = self.get_os_masses()
            ax1.hist(m1_os, bins=bins, histtype='step', linewidth=2,
                    color='#2E86AB', label='OS', alpha=0.9)
        if self.ls_data is not None:
            m1_ls, _ = self.get_ls_masses()
            ax1.hist(m1_ls, bins=bins, histtype='step', linewidth=2,
                    color='#F18F01', label='LS', alpha=0.9)
        ax1.set_xlabel(r'$m_1$ [GeV/$c^2$]')
        ax1.set_ylabel('Events / bin')
        ax1.set_title('$m_1$ Projection', fontweight='bold')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        # m2
        if self.os_data is not None:
            _, m2_os = self.get_os_masses()
            ax2.hist(m2_os, bins=bins, histtype='step', linewidth=2,
                    color='#A23B72', label='OS', alpha=0.9)
        if self.ls_data is not None:
            _, m2_ls = self.get_ls_masses()
            ax2.hist(m2_ls, bins=bins, histtype='step', linewidth=2,
                    color='#C73E1D', label='LS', alpha=0.9)
        ax2.set_xlabel(r'$m_2$ [GeV/$c^2$]')
        ax2.set_ylabel('Events / bin')
        ax2.set_title('$m_2$ Projection', fontweight='bold')
        ax2.legend()
        ax2.grid(alpha=0.3)
        
        plt.suptitle('1D Mass Projections', fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig, (ax1, ax2)
    
    # ========== Correlation Plots ==========
    
    def plot_correlation(self, gridsize: int = 40, figsize: Tuple[float, float] = (10, 9)):
        """Create correlation plot between m1 and m2 for OS data"""
        if self.os_data is None:
            raise ValueError("No OS data available")
        
        m1_os, m2_os = self.get_os_masses()
        
        fig, ax = plt.subplots(figsize=figsize)
        
        hb = ax.hexbin(m1_os, m2_os, gridsize=gridsize, cmap='YlOrRd',
                      mincnt=1, linewidths=0.2)
        
        cb = plt.colorbar(hb, ax=ax)
        cb.set_label('Events per bin', rotation=270, labelpad=20, fontsize=12)
        
        corr = np.corrcoef(m1_os, m2_os)[0, 1]
        
        lims = [min(m1_os.min(), m2_os.min()), max(m1_os.max(), m2_os.max())]
        ax.plot(lims, lims, 'r--', alpha=0.5, linewidth=2, label='$m_1 = m_2$')
        
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_title(f'OS Mass Correlation (ρ = {corr:.3f})', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper left', fontsize=11)
        ax.set_aspect('equal')
        ax.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        return fig, ax
    
    # ========== Kinematic Plots ==========
    
    def plot_kinematics(self, data_type: str = 'os', figsize: Tuple[float, float] = (14, 10)):
        """Plot kinematic distributions"""
        data = self.os_data if data_type.lower() == 'os' else self.ls_data
        if data is None:
            raise ValueError(f"No {data_type.upper()} data available")
        
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        
        # pT1
        axes[0, 0].hist(data['pt1'], bins=200, histtype='stepfilled',
                       color='#2E86AB', alpha=0.6, edgecolor='black')
        axes[0, 0].set_xlabel(r'$p_{T,1}$ [GeV/$c$]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_xlim(0.0, 24.0)
        axes[0, 0].set_title(f'{data_type.upper()}: Particle 1 $p_T$', fontweight='bold')
        axes[0, 0].grid(alpha=0.3)
        
        # pT2
        axes[0, 1].hist(data['pt2'], bins=200, histtype='stepfilled',
                       color='#A23B72', alpha=0.6, edgecolor='black')
        axes[0, 1].set_xlabel(r'$p_{T,2}$ [GeV/$c$]')
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_xlim(0.0, 24.0)
        axes[0, 1].set_title(f'{data_type.upper()}: Particle 2 $p_T$', fontweight='bold')
        axes[0, 1].grid(alpha=0.3)
        
        # Pair pT
        axes[1, 0].hist(data['pairPt'], bins=100, histtype='stepfilled',
                       color='#F18F01', alpha=0.6, edgecolor='black')
        axes[1, 0].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_xlim(0.0, 35.0)
        axes[1, 0].set_title(f'{data_type.upper()}: Pair $p_T$', fontweight='bold')
        axes[1, 0].grid(alpha=0.3)
        
        # Pair mass
        axes[1, 1].hist(data['pairMass'], bins=500, histtype='stepfilled',
                       color='#C73E1D', alpha=0.6, edgecolor='black')
        axes[1, 1].set_xlabel(r'$M_{pair}$ [GeV/$c^2$]')
        axes[1, 1].set_xlim(3.0, 12.0)
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title(f'{data_type.upper()}: Pair Mass', fontweight='bold')
        axes[1, 1].grid(alpha=0.3)
        
        plt.suptitle(f'{data_type.upper()} Kinematic Distributions', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig, axes
    
    # ========== Convenience Methods for Quick Display ==========
    
    def show_2d(self, bins: int = 50):
        """Quick display of 2D distributions"""
        fig, axes = self.plot_2d_both(bins=bins)
        plt.show()
        return fig, axes
    
    def show_projections(self, bins: int = 60):
        """Quick display of 1D projections"""
        fig, axes = self.plot_projections(bins=bins)
        plt.show()
        return fig, axes
    
    def show_correlation(self):
        """Quick display of correlation plot"""
        fig, ax = self.plot_correlation()
        plt.show()
        return fig, ax
    
    def show_kinematics(self, data_type: str = 'os'):
        """Quick display of kinematic distributions"""
        fig, axes = self.plot_kinematics(data_type=data_type)
        plt.show()
        return fig, axes
    
    def show_summary(self, bins: int = 40, figsize: Tuple[float, float] = (16, 12)):
        """Create comprehensive summary plot"""
        if self.os_data is None or self.ls_data is None:
            print("Warning: Need both OS and LS data for full summary")
        
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.35)
        
        m1_os, m2_os = self.get_os_masses()
        
        # 2D OS
        ax1 = fig.add_subplot(gs[0, 0])
        h, xe, ye = np.histogram2d(m1_os, m2_os, bins=bins)
        im = ax1.imshow(h.T, origin='lower', aspect='auto',
                       extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap='viridis')
        ax1.set_title('OS 2D', fontweight='bold')
        ax1.set_xlabel(r'$m_1$')
        ax1.set_ylabel(r'$m_2$')
        plt.colorbar(im, ax=ax1, fraction=0.046)
        
        if self.ls_data is not None:
            m1_ls, m2_ls = self.get_ls_masses()
            
            # 2D LS
            ax2 = fig.add_subplot(gs[0, 1])
            h, xe, ye = np.histogram2d(m1_ls, m2_ls, bins=bins)
            im = ax2.imshow(h.T, origin='lower', aspect='auto',
                           extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap='plasma')
            ax2.set_title('LS 2D', fontweight='bold')
            ax2.set_xlabel(r'$m_1$')
            ax2.set_ylabel(r'$m_2$')
            plt.colorbar(im, ax=ax2, fraction=0.046)
        
        # Correlation
        ax3 = fig.add_subplot(gs[0, 2])
        ax3.hexbin(m1_os, m2_os, gridsize=30, cmap='YlOrRd', mincnt=1)
        corr = np.corrcoef(m1_os, m2_os)[0, 1]
        ax3.set_title(f'Correlation (ρ={corr:.2f})', fontweight='bold')
        ax3.set_xlabel(r'$m_1$')
        ax3.set_ylabel(r'$m_2$')
        
        # Projections
        ax4 = fig.add_subplot(gs[1, 0])
        ax4.hist(m1_os, bins=50, histtype='step', linewidth=2, color='blue', label='OS', alpha=0.8)
        if self.ls_data is not None:
            ax4.hist(m1_ls, bins=50, histtype='step', linewidth=2, color='red', label='LS', alpha=0.8)
        ax4.set_xlabel(r'$m_1$ [GeV/$c^2$]')
        ax4.set_ylabel('Events')
        ax4.set_title('$m_1$ Projection', fontweight='bold')
        ax4.legend()
        ax4.grid(alpha=0.3)
        
        ax5 = fig.add_subplot(gs[1, 1])
        ax5.hist(m2_os, bins=50, histtype='step', linewidth=2, color='blue', label='OS', alpha=0.8)
        if self.ls_data is not None:
            ax5.hist(m2_ls, bins=50, histtype='step', linewidth=2, color='red', label='LS', alpha=0.8)
        ax5.set_xlabel(r'$m_2$ [GeV/$c^2$]')
        ax5.set_ylabel('Events')
        ax5.set_title('$m_2$ Projection', fontweight='bold')
        ax5.legend()
        ax5.grid(alpha=0.3)
        
        # OS-LS difference
        if self.ls_data is not None:
            ax6 = fig.add_subplot(gs[1, 2])
            bins_diff = np.linspace(min(m1_os.min(), m1_ls.min()),
                                   max(m1_os.max(), m1_ls.max()), 50)
            h_os, _ = np.histogram(m1_os, bins=bins_diff)
            h_ls, _ = np.histogram(m1_ls, bins=bins_diff)
            bin_centers = (bins_diff[:-1] + bins_diff[1:]) / 2
            ax6.bar(bin_centers, h_os - h_ls, width=np.diff(bins_diff)[0],
                   color='green', alpha=0.7, edgecolor='black')
            ax6.axhline(0, color='black', linestyle='--', linewidth=1)
            ax6.set_xlabel(r'$m_1$ [GeV/$c^2$]')
            ax6.set_ylabel('OS - LS')
            ax6.set_title('Subtracted', fontweight='bold')
            ax6.grid(alpha=0.3)
        
        # Kinematics
        ax7 = fig.add_subplot(gs[2, 0])
        ax7.hist(self.os_data['pt1'], bins=40, histtype='stepfilled',
                color='steelblue', alpha=0.7, edgecolor='black')
        ax7.set_xlabel(r'$p_{T,1}$ [GeV/$c$]')
        ax7.set_ylabel('Events')
        ax7.set_title('Single $p_T$', fontweight='bold')
        ax7.grid(alpha=0.3)
        
        ax8 = fig.add_subplot(gs[2, 1])
        ax8.hist(self.os_data['pairPt'], bins=40, histtype='stepfilled',
                color='coral', alpha=0.7, edgecolor='black')
        ax8.set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        ax8.set_ylabel('Events')
        ax8.set_title('Pair $p_T$', fontweight='bold')
        ax8.grid(alpha=0.3)
        
        ax9 = fig.add_subplot(gs[2, 2])
        ax9.hist(self.os_data['pairMass'], bins=500, histtype='stepfilled',
                color='orchid', alpha=0.7, edgecolor='black')
        ax9.set_xlabel(r'$M_{pair}$ [GeV/$c^2$]')
        ax9.set_ylabel('Events')
        ax9.set_title('Pair Mass', fontweight='bold')
        ax9.set_xlim(4.0, 12.0)
        
        ax9.grid(alpha=0.3)
        
        plt.suptitle('D Meson Analysis Summary', fontsize=18, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.show()
        return fig
    
    # ========== Save Methods ==========
    
    def save_all_plots(self, output_dir: str = 'plots', dpi: int = 300):
        """Save all plots to a directory"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        print(f"Saving plots to {output_path}/")
        print("-" * 60)
        
        # 2D distributions
        if self.os_data is not None or self.ls_data is not None:
            fig, _ = self.plot_2d_both()
            fig.savefig(output_path / '2d_distributions.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ 2d_distributions.pdf")
        
        # Projections
        if self.os_data is not None or self.ls_data is not None:
            fig, _ = self.plot_projections()
            fig.savefig(output_path / 'projections.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ projections.pdf")
        
        # Correlation
        if self.os_data is not None:
            fig, _ = self.plot_correlation()
            fig.savefig(output_path / 'correlation.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ correlation.pdf")
        
        # OS kinematics
        if self.os_data is not None:
            fig, _ = self.plot_kinematics('os')
            fig.savefig(output_path / 'kinematics_os.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ kinematics_os.pdf")
        
        # LS kinematics
        if self.ls_data is not None:
            fig, _ = self.plot_kinematics('ls')
            fig.savefig(output_path / 'kinematics_ls.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ kinematics_ls.pdf")
        
        # Summary
        if self.os_data is not None:
            fig = self.show_summary()
            fig.savefig(output_path / 'summary.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ summary.pdf")
        
        print("-" * 60)
        print(f"All plots saved to {output_path}/")


# ========== Utility Functions ==========

def quick_look(workspace_file: str):
    """
    Quick visualization of all data
    
    Usage:
        quick_look('workspace.root')
    """
    viz = FitVisualizer(workspace_file)
    viz.show_2d()
    viz.show_projections()
    viz.show_correlation()
    viz.show_kinematics()
    return viz


def main():
    """Command line interface (kept for backward compatibility)"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize 2D invariant mass fits')
    parser.add_argument('--workspace', type=str, 
                       default='fit_workspace.root',
                       help='RooFit workspace file')
    parser.add_argument('--output-dir', type=str,
                       default='plots',
                       help='Output directory for plots')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("2D Invariant Mass Fit Visualizer")
    print("=" * 60)
    
    viz = FitVisualizer(args.workspace)
    viz.save_all_plots(output_dir=args.output_dir)
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()