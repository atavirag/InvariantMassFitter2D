#!/usr/bin/env python3
"""
visualize_raw_data.py - Visualization functions for raw data

This module provides a clean interface for creating publication-quality plots
from RooFit workspaces. Designed for use in Jupyter notebooks.

Usage in Jupyter:
    from visualize_raw_data import RawDataVisualizer
    
    viz = RawDataVisualizer(os_file='../../data/os_pairs.root',
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


class RawDataVisualizer:
    """
    Visualize 2D invariant mass and kinematic distributions from the pair tree
    
    This class provides methods to create various plots from ROOT files
    containing RooFit workspaces and data. All methods return matplotlib
    figure and axes objects for easy customization in notebooks.
    
    Attributes:
        os_file (Path): Path to OS data ROOT file
        ls_file (Path): Path to LS data ROOT file
        os_data (dict): Dictionary containing OS data arrays
        ls_data (dict): Dictionary containing LS data arrays
    """
    
    def __init__(self, os_file: Optional[str] = None, 
                 ls_file: Optional[str] = None, verbose: bool = True):
        """
        Initialize visualizer
        
        Parameters:
            os_file: Path to OS data ROOT file (optional)
            ls_file: Path to LS data ROOT file (optional)
            verbose: Print loading information
        """
        self.os_file = Path(os_file) if os_file else None
        self.ls_file = Path(ls_file) if ls_file else None
        self.verbose = verbose
        
        # Data containers
        self.os_data = None
        self.ls_data = None
        self.has_pyroot = False
        
        # Check for PyROOT
        try:
            import ROOT
            self.has_pyroot = True
        except ImportError:
            pass
        
        # Load everything
        self._load_data()
        
        # Print summary
        if verbose:
            self._print_summary()
    
    def _print_summary(self):
        """Print a summary of what was loaded"""
        print("\n" + "=" * 60)
        print("RawDataVisualizer Summary")
        print("=" * 60)
        print(f"OS Data: {'✓ Loaded' if self.os_data else '✗ Not available'}")
        print(f"LS Data: {'✓ Loaded' if self.ls_data else '✗ Not available'}")
        print(f"PyROOT: {'✓ Available' if self.has_pyroot else '✗ Not available (limited workspace access)'}")
        if not self.has_pyroot:
            print("\nTip: Install PyROOT for full RooFit workspace access:")
            print("     pip install root")
        print("=" * 60 + "\n")
    
    def _print(self, msg: str):
        """Print if verbose mode is on"""
        if self.verbose:
            print(msg)
    
    def _load_data(self):
        """Load data from ROOT files using uproot"""
        # Load OS data
        if self.os_file and self.os_file.exists():
            try:
                with uproot.open(self.os_file) as f:
                    if 'pairs' in f:
                        tree = f['pairs']
                        self.os_data = tree.arrays(
                            ['m1', 'm2', 'mbar1', 'mbar2', 'y1', 'y2', 'phi1', 'phi2', 'pairType', 
                             'pt1', 'pt2', 'pairPt', 'pairMass'], 
                            library='np'
                        )
                        self._print(f"✓ Loaded {len(self.os_data['m1'])} OS events")
                    else:
                        self._print(f"Warning: No 'pairs' tree in {self.os_file}")
            except Exception as e:
                self._print(f"Error loading OS file: {e}")
        
        # Load LS data
        if self.ls_file and self.ls_file.exists():
            try:
                with uproot.open(self.ls_file) as f:
                    if 'pairs' in f:
                        tree = f['pairs']
                        self.ls_data = tree.arrays(
                            ['m1', 'm2', 'mbar1', 'mbar2', 'y1', 'y2', 'phi1', 'phi2', 'pairType',
                             'pt1', 'pt2', 'pairPt', 'pairMass'], 
                            library='np'
                        )
                        self._print(f"✓ Loaded {len(self.ls_data['m1'])} LS events")
                    else:
                        self._print(f"Warning: No 'pairs' tree in {self.ls_file}")
            except Exception as e:
                self._print(f"Error loading LS file: {e}")

    def _prepare_masses(self, data: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract correct mass combinations based on pairType
        
        The pairType is a bitmask:
        - bit 0: DD (both D mesons)
        - bit 1: DbarDbar (both anti-D mesons)
        - bit 2: DDbar (D and anti-D)
        - bit 3: DbarD (anti-D and D)
        """
        pairType = data['pairType']
        m1_vals = np.zeros(len(pairType))
        m2_vals = np.zeros(len(pairType))
        
        # DD (bit 0)
        mask_dd = (pairType & (1 << 0)) != 0
        m1_vals[mask_dd] = data['m1'][mask_dd]
        m2_vals[mask_dd] = data['m2'][mask_dd]
        
        # DbarDbar (bit 1)
        mask_dbdb = (pairType & (1 << 1)) != 0
        m1_vals[mask_dbdb] = data['mbar1'][mask_dbdb]
        m2_vals[mask_dbdb] = data['mbar2'][mask_dbdb]
        
        # DDbar (bit 2)
        mask_ddb = (pairType & (1 << 2)) != 0
        m1_vals[mask_ddb] = data['m1'][mask_ddb]
        m2_vals[mask_ddb] = data['mbar2'][mask_ddb]
        
        # DbarD (bit 3)
        mask_dbd = (pairType & (1 << 3)) != 0
        m1_vals[mask_dbd] = data['mbar1'][mask_dbd]
        m2_vals[mask_dbd] = data['m2'][mask_dbd]
        
        return m1_vals, m2_vals
    
    def get_os_masses(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get OS mass arrays (m1, m2)"""
        if self.os_data is None:
            raise ValueError("No OS data loaded")
        return self._prepare_masses(self.os_data)
    
    def get_ls_masses(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get LS mass arrays (m1, m2)"""
        if self.ls_data is None:
            raise ValueError("No LS data loaded")
        return self._prepare_masses(self.ls_data)
    
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
    def add_deltas(self, data_type: str = 'os'):
        data = self.os_data if data_type.lower() == 'os' else self.ls_data
        if data is None:
            raise ValueError(f"No {data_type.upper()} data available")
        # delta y
        data['deltaY'] = (np.abs(data['y1'] - data['y2']))
        # delta phi
        data['deltaPhi'] = np.abs((data['phi1'] - data['phi2'] + np.pi) % (2 * np.pi) - np.pi)
    
    def plot_kinematics(self, data_type: str = 'os', figsize: Tuple[float, float] = (14, 10)):
        """Plot kinematic distributions"""
        data = self.os_data if data_type.lower() == 'os' else self.ls_data
        if data is None:
            raise ValueError(f"No {data_type.upper()} data available")

        self.add_deltas(data_type)
        fig, axes = plt.subplots(3, 3, figsize=figsize)
        
        # pT1
        axes[0, 0].hist(data['pt1'], bins=120, histtype='stepfilled',
                       color='#2E86AB', alpha=0.6, edgecolor='black')
        axes[0, 0].set_xlabel(r'$p_{T,1}$ [GeV/$c$]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_yscale('log')
        axes[0, 0].set_xlim(0.0, 24.0)
        axes[0, 0].set_title(f'{data_type.upper()}: Particle 1 $p_T$', fontweight='bold')
        axes[0, 0].grid(alpha=0.3)
        
        # pT2
        axes[0, 1].hist(data['pt2'], bins=120, histtype='stepfilled',
                       color='#A23B72', alpha=0.6, edgecolor='black')
        axes[0, 1].set_xlabel(r'$p_{T,2}$ [GeV/$c$]')
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_xlim(0.0, 24.0)
        axes[0, 1].set_yscale('log')
        axes[0, 1].set_title(f'{data_type.upper()}: Particle 2 $p_T$', fontweight='bold')
        axes[0, 1].grid(alpha=0.3)
        
        # Pair pT
        axes[0, 2].hist(data['pairPt'], bins=100, histtype='stepfilled',
                       color='#F18F01', alpha=0.6, edgecolor='black')
        axes[0, 2].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[0, 2].set_ylabel('Events')
        axes[0, 2].set_xlim(0.0, 35.0)
        axes[0, 2].set_yscale('log')
        axes[0, 2].set_title(f'{data_type.upper()}: Pair $p_T$', fontweight='bold')
        axes[0, 2].grid(alpha=0.3)
        
        # y1
        axes[1, 0].hist(data['y1'], bins=50, histtype='stepfilled',
                       color='#2E86AB', alpha=0.6, edgecolor='black')
        axes[1, 0].set_xlabel(r'$y_{1}$')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_title(f'{data_type.upper()}: $y_{1}$', fontweight='bold')
        axes[1, 0].grid(alpha=0.3)

        # y2
        axes[1, 1].hist(data['y2'], bins=50, histtype='stepfilled',
                       color='#A23B72', alpha=0.6, edgecolor='black')
        axes[1, 1].set_xlabel(r'$y_{2}$')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title(f'{data_type.upper()}: $y_{2}$', fontweight='bold')
        axes[1, 1].grid(alpha=0.3)

        # delta y
        axes[1, 2].hist(data['deltaY'], bins=50, histtype='stepfilled',
                       color='#F18F01', alpha=0.6, edgecolor='black')
        axes[1, 2].set_xlabel(r'$\Delta y$')
        axes[1, 2].set_ylabel('Events')
        axes[1, 2].set_title(f'{data_type.upper()}: $\Delta y$', fontweight='bold')
        axes[1, 2].grid(alpha=0.3)

        # phi1
        axes[2, 0].hist(data['phi1'], bins=50, histtype='stepfilled',
                       color='#2E86AB', alpha=0.6, edgecolor='black')
        axes[2, 0].set_xlabel(r'$\phi_{1}$')
        axes[2, 0].set_ylabel('Events')
        axes[2, 0].set_title(f'{data_type.upper()}: $\phi_{1}$', fontweight='bold')
        axes[2, 0].grid(alpha=0.3)

        # phi2
        axes[2, 1].hist(data['phi2'], bins=50, histtype='stepfilled',
                       color='#A23B72', alpha=0.6, edgecolor='black')
        axes[2, 1].set_xlabel(r'$\phi_{2}$')
        axes[2, 1].set_ylabel('Events')
        axes[2, 1].set_title(f'{data_type.upper()}: $\phi_{2}$', fontweight='bold')
        axes[2, 1].grid(alpha=0.3)

        # delta phi
        axes[2, 2].hist(data['deltaPhi'], bins=50, histtype='stepfilled',
                       color='#F18F01', alpha=0.6, edgecolor='black')
        axes[2, 2].set_xlabel(r'$\Delta \phi$')
        axes[2, 2].set_ylabel('Events')
        axes[2, 2].set_title(f'{data_type.upper()}: $\Delta \phi$', fontweight='bold')
        axes[2, 2].grid(alpha=0.3)
        
        plt.suptitle(f'{data_type.upper()} Kinematic Distributions', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig, axes
    
    def plot_pair_kinematics(self, compare: bool = True, figsize: Tuple[float, float] = (14, 10)):
        """Plot pair kinematic distributions, optionally comparing OS and LS"""

        self.add_deltas('os')
        self.add_deltas('ls')

        fig, axes = plt.subplots(2, 2, figsize=figsize)
    
        # Define common bin edges for fair comparison
        pt_bins = np.linspace(0.0, 36.0, 91)
        mass_bins = np.linspace(0.0, 50.0, 501)
        deltay_bins = 50
        deltaphi_bins = 50
    
        # Pair pT
        if self.os_data is not None:
            axes[0, 0].hist(self.os_data['pairPt'], bins=pt_bins, histtype='step', linewidth=2,
                        color='#2E86AB', alpha=0.9, label='OS')
        if compare and self.ls_data is not None:
            axes[0, 0].hist(self.ls_data['pairPt'], bins=pt_bins, histtype='step', linewidth=2,
                       color='#F18F01', alpha=0.9, label='LS')
        axes[0, 0].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_yscale('log')
        axes[0, 0].set_xlim(0.0, 36.0)
        axes[0, 0].set_title('Pair $p_T$', fontweight='bold')
        axes[0, 0].legend(frameon=True, shadow=True)
        axes[0, 0].grid(alpha=0.3)
    
        # Pair mass
        if self.os_data is not None:
            axes[0, 1].hist(self.os_data['pairMass'], bins=mass_bins, histtype='step', linewidth=2,
                       color='#C73E1D', alpha=0.9, label='OS')
        if compare and self.ls_data is not None:
            axes[0, 1].hist(self.ls_data['pairMass'], bins=mass_bins, histtype='step', linewidth=2,
                       color='#6A4C93', alpha=0.9, label='LS')
        axes[0, 1].set_xlabel(r'$M_{pair}$ [GeV/$c^2$]')
        axes[0, 1].set_xlim(3.0, 12.0)
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_title('Pair Mass', fontweight='bold')
        axes[0, 1].legend(frameon=True, shadow=True)
        axes[0, 1].grid(alpha=0.3)
    
        # deltaY
        if self.os_data is not None and 'deltaY' in self.os_data:
            axes[1, 0].hist(self.os_data['deltaY'], bins=deltay_bins, histtype='step', linewidth=2,
                       color='#A23B72', alpha=0.9, label='OS')
        if compare and self.ls_data is not None and 'deltaY' in self.ls_data:
            axes[1, 0].hist(self.ls_data['deltaY'], bins=deltay_bins, histtype='step', linewidth=2,
                       color='#1982C4', alpha=0.9, label='LS')
        axes[1, 0].set_xlabel(r'$\Delta y$')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_title('$\Delta y$', fontweight='bold')
        axes[1, 0].legend(frameon=True, shadow=True)
        axes[1, 0].grid(alpha=0.3)
    
        # deltaPhi
        if self.os_data is not None and 'deltaPhi' in self.os_data:
            axes[1, 1].hist(self.os_data['deltaPhi'], bins=deltaphi_bins, histtype='step', linewidth=2,
                       color='#F18F01', alpha=0.9, label='OS')
        if compare and self.ls_data is not None and 'deltaPhi' in self.ls_data:
            axes[1, 1].hist(self.ls_data['deltaPhi'], bins=deltaphi_bins, histtype='step', linewidth=2,
                       color='#8338EC', alpha=0.9, label='LS')
        axes[1, 1].set_xlabel(r'$\Delta\phi$')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title('$\Delta\phi$', fontweight='bold')
        axes[1, 1].legend(frameon=True, shadow=True)
        axes[1, 1].grid(alpha=0.3)
    
        plt.suptitle('OS vs LS Pair Kinematics Comparison' if compare else 'Pair Kinematics', 
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
    
    def show_pair_kinematics(self):
        """Quick display of kinematic distributions"""
        fig, axes = self.plot_pair_kinematics()
        plt.show()
        return fig, axes
    
    def show_kinematics(self, data_type: str = 'os'):
        """Quick display of kinematic distributions"""
        fig, axes = self.plot_kinematics(data_type=data_type)
        plt.show()
        return fig, axes
    
    def show_pair_kinematics(self):
        """Quick display of kinematic distributions"""
        fig, axes = self.plot_pair_kinematics()
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
        ax7.hist(self.os_data['pt1'], bins=40, histtype='step', linewidth=2, label='OS',
                color='blue', alpha=0.8)
        if self.ls_data is not None:
            ax7.hist(self.ls_data['pt1'], bins=40, histtype='step', linewidth=2, color='red', label='LS', alpha=0.8)
        ax7.set_xlabel(r'$p_{T,1}$ [GeV/$c$]')
        ax7.set_ylabel('Events')
        ax7.set_yscale('log')
        ax7.set_title('Single $p_{T, 1}$', fontweight='bold')
        ax7.grid(alpha=0.3)
        ax7.legend()
        ax7.set_xlim(0.0, 24.0)
        
        ax8 = fig.add_subplot(gs[2, 1])
        ax8.hist(self.os_data['pt2'], bins=40, histtype='step', linewidth=2, label='OS',
                color='blue', alpha=0.8)
        if self.ls_data is not None:
            ax8.hist(self.ls_data['pt2'], bins=40, histtype='step', linewidth=2, color='red', label='LS', alpha=0.8)
        ax8.set_xlabel(r'$p_{T,2}$ [GeV/$c$]')
        ax8.set_ylabel('Events')
        ax8.set_yscale('log')
        ax8.set_title('Single $p_{T, 2}$', fontweight='bold')
        ax8.grid(alpha=0.3)
        ax8.legend()
        ax8.set_xlim(0.0, 24.0)
        
        ax9 = fig.add_subplot(gs[2, 2])
        ax9.hist(self.os_data['pairPt'], bins=40, histtype='step', linewidth=2, label='OS',
                color='blue', alpha=0.8)
        if self.ls_data is not None:
            ax9.hist(self.ls_data['pairPt'], bins=40, histtype='step', linewidth=2, color='red', label='LS', alpha=0.8)
        ax9.set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        ax9.set_ylabel('Events')
        ax9.set_yscale('log')
        ax9.set_title('Pair $p_T$', fontweight='bold')
        ax9.legend()
        ax9.set_xlim(0.0, 36.0)
        
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

        # Pair kinematics
        if self.ls_data is not None and self.os_data is not None:
            fig, _ = self.plot_pair_kinematics()
            fig.savefig(output_path / 'pair_kinematics.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ pair_kinematics.pdf")
        
        # Summary
        if self.os_data is not None:
            fig = self.show_summary()
            fig.savefig(output_path / 'summary.pdf', dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            print("✓ summary.pdf")
        
        print("-" * 60)
        print(f"All plots saved to {output_path}/")


# ========== Utility Functions ==========

def quick_look(os_file: Optional[str] = None, 
               ls_file: Optional[str] = None):
    """
    Quick visualization of all data
    
    Usage:
        quick_look('os_pairs.root', 'ls_pairs.root')
    """
    viz = RawDataVisualizer(os_file, ls_file)
    viz.show_2d()
    viz.show_projections()
    viz.show_correlation()
    viz.show_kinematics()
    viz.show_pair_kinematics()
    return viz


def main():
    """Command line interface (kept for backward compatibility)"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize raw data invariant mass and kinematic distributions')
    parser.add_argument('--os-file', type=str,
                       default='os_pairs.root',
                       help='OS data file')
    parser.add_argument('--ls-file', type=str,
                       default='ls_pairs.root',
                       help='LS data file')
    parser.add_argument('--output-dir', type=str,
                       default='plots',
                       help='Output directory for plots')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Raw data visualizer")
    print("=" * 60)
    
    viz = RawDataVisualizer(args.os_file, args.ls_file)
    viz.save_all_plots(output_dir=args.output_dir)
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()