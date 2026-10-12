#!/usr/bin/env python3
"""
correlation_functions.py - Visualization functions for D meson pair correlations

This module provides a clean interface for creating publication-quality plots
from RooFit workspaces. Designed for use in Jupyter notebooks.

Usage in Jupyter:
    from correlation_functions import CorrelationFunctions

    corr = CorrelationFunctions(file='../../data/os_pairs.root')

    # Create individual plots
    fig, ax = corr.plot_2d_os(bins=50)
    plt.show()

    # Or use convenience methods
    corr.show_kinematics()
    corr.show_pair_kinematics()
    corr.show_2d_delta_eta_phi()
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import seaborn as sns
import uproot
import ROOT
from mpl_toolkits.mplot3d import Axes3D
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


class CorrelationFunctions:
    """
    Visualize 2D invariant mass and kinematic distributions from the pair tree

    This class provides methods to create various plots from ROOT files
    containing RooFit workspaces and data. All methods return matplotlib
    figure and axes objects for easy customization in notebooks.

    Attributes:
        file (Path): Path to data ROOT file
        os_data (dict): Dictionary containing data arrays
    """

    # --- pT-dependent signal/sideband mass windows, from SetPtRanges() ---
    PT_BINS  = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 16., 24.0])
    SIGN_LOW = np.array([1.83632,1.82996,1.82514,1.82106,1.81705,1.81195,1.80820,
                          1.80400,1.80022,1.79685,1.79171,1.79171,1.79171])
    SIGN_UPP = np.array([1.87963,1.88585,1.89095,1.89453,1.89887,1.90420,1.90812,
                          1.91201,1.91714,1.92061,1.92627,1.92627,1.92627])
    LSB_LOW  = np.array([1.78218,1.70601,1.74287,1.72924,1.71477,1.69663,1.68331,
                          1.66899,1.65406,1.64215,1.62350,1.62350,1.62350])
    LSB_UPP  = np.array([1.80384,1.78804,1.77577,1.76597,1.75568,1.74276,1.73327,
                          1.72300,1.71253,1.70403,1.69079,1.69079,1.69079])
    RSB_LOW  = np.array([1.91211,1.92776,1.94032,1.94963,1.96023,1.97339,1.98305,
                          1.99301,2.00483,2.01342,2.02720,2.02720,2.02720])
    RSB_UPP  = np.array([1.93377,1.95570,1.97323,1.98636,2.00114,2.01951,2.03301,
                          2.04701,2.06329,2.07530,2.09448,2.09448,2.09448])

    def __init__(self, file: Optional[str] = None, verbose: bool = True):
        """
        Initialize visualizer

        Parameters:
            file: Path to data ROOT file for D0D0bar correlations (optional)
            verbose: Print loading information
        """
        self.file = Path(file) if file else None
        self.verbose = verbose
        
        # Data containers
        self.data = None
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
        print("CorrelationFunctions Summary")
        print("=" * 60)
        print(f"Data: {'✓ Loaded' if self.data else '✗ Not available'}")
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
        if self.file and self.file.exists():
            try:
                with uproot.open(self.file) as f:
                    if 'pairs' in f:
                        tree = f['pairs']
                        self.data = tree.arrays(
                            ['m1', 'm2', 'mbar1', 'mbar2', 'y1', 'y2', 'phi1', 'phi2', 'pairType',
                             'pt1', 'pt2', 'pairPt', 'pairMass'],
                            library='np'
                        )
                        self._print(f"✓ Loaded {len(self.data['m1'])} events")
                    else:
                        self._print(f"Warning: No 'pairs' tree in {self.file}")
            except Exception as e:
                self._print(f"Error loading file: {e}")

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
        if self.data is None:
            raise ValueError("No data loaded")
        m1_vals, m2_vals = self._prepare_masses(self.data)
        pairType = self.data['pairType']
        os_mask = (pairType & (1 << 2)) != 0 | (pairType & (1 << 3)) != 0
        return m1_vals[os_mask], m2_vals[os_mask]

    def get_ls_masses(self) -> Tuple[np.ndarray, np.ndarray]:
        if self.data is None:
            raise ValueError("No data loaded")
        m1_vals, m2_vals = self._prepare_masses(self.data)
        pairType = self.data['pairType']
        ls_mask = (pairType & (1 << 0)) != 0 | (pairType & (1 << 1)) != 0
        return m1_vals[ls_mask], m2_vals[ls_mask]

    def _get_pt_bin_index(self, pt: np.ndarray) -> np.ndarray:
        """Map each pT value to its bin index in PT_BINS (clipped to valid range)."""
        idx = np.digitize(pt, self.PT_BINS) - 1
        return np.clip(idx, 0, len(self.PT_BINS) - 2)

    def get_signal_sideband_masks(self, data_type: str = 'os', pt_var: str = 'pt1'):
        """
        Build boolean masks (full-length, matching self.data) for the signal
        region, left sideband, and right sideband, using pT-dependent mass
        windows exactly as in SetPtRanges().

        Parameters:
            data_type: 'os' or 'ls'
            pt_var: which pT column defines the bin (default 'pt1', the D
                candidate pT — change to 'pairPt' if that's what the C++
                ranges were actually binned in)
        """
        if self.data is None:
            raise ValueError("No data loaded")

        m1_vals, _ = self._prepare_masses(self.data)
        pt = self.data[pt_var]
        type_mask = self._full_length_type_mask(data_type)
        pt_idx = self._get_pt_bin_index(pt)

        signal_mask = (m1_vals >= self.SIGN_LOW[pt_idx]) & (m1_vals <= self.SIGN_UPP[pt_idx]) & type_mask
        lsb_mask    = (m1_vals >= self.LSB_LOW[pt_idx])  & (m1_vals <= self.LSB_UPP[pt_idx])  & type_mask
        rsb_mask    = (m1_vals >= self.RSB_LOW[pt_idx])  & (m1_vals <= self.RSB_UPP[pt_idx])  & type_mask
        sideband_mask = lsb_mask | rsb_mask

        return signal_mask, sideband_mask, lsb_mask, rsb_mask

    def select_signal_region(self, data_type: str = 'os', pt_var: str = 'pt1'):
        """Constrain self.data to the pT-dependent signal-region mask."""
        signal_mask, _, _, _ = self.get_signal_sideband_masks(data_type, pt_var)
        self._apply_mask(signal_mask)

    def select_sideband_region(self, data_type: str = 'os', pt_var: str = 'pt1', side: str = 'both'):
        """
        Constrain self.data to the pT-dependent sideband mask.
        side: 'left', 'right', or 'both' (LSB ∪ RSB, matching the C++ combined SB).
        """
        _, sideband_mask, lsb_mask, rsb_mask = self.get_signal_sideband_masks(data_type, pt_var)
        mask = {'left': lsb_mask, 'right': rsb_mask, 'both': sideband_mask}[side]
        self._apply_mask(mask)

    def _apply_mask(self, mask: np.ndarray):
        """Shared helper: backs up full data (if not already) and applies a mask."""
        if self._original_data is None:
            self._original_data = {key: val.copy() for key, val in self.data.items()}
        self.data = {key: val[mask] for key, val in self.data.items()}
    
    # ========== 2D Distribution Plots ==========
    
    def plot_2d(self, bins: int = 50, cmap: str = 'viridis', data_type:str='os',
                   figsize: Tuple[float, float] = (8, 7)) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create 2D mass distribution plot
        
        Parameters:
            bins: Number of bins for histogram
            cmap: Colormap name
            data_type: Type of data to plot ('os' or 'ls')
            figsize: Figure size (width, height)
            
        Returns:
            fig, ax: Matplotlib figure and axes
        """
        if self.data is None:
            raise ValueError("No data available")
        
        if data_type.lower() == 'os':
            m1, m2 = self.get_os_masses()
        elif data_type.lower() == 'ls':
            m1, m2 = self.get_ls_masses()
        else:
            raise ValueError(f"Unknown data_type '{data_type}'. Use 'os' or 'ls'.")

        fig, ax = plt.subplots(figsize=figsize)
        
        h, xedges, yedges = np.histogram2d(m1, m2, bins=bins)
        
        im = ax.imshow(h.T, origin='lower', aspect='auto', 
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                      cmap=cmap, interpolation='nearest')
        
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_title(f'{data_type.upper()} Sign', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Events', rotation=270, labelpad=20)
        
        # Diagonal reference line
        ax.plot([xedges[0], xedges[-1]], [yedges[0], yedges[-1]], 
                'r--', alpha=0.5, linewidth=1.5, label='$m_1 = m_2$')
        ax.legend(loc='upper left')
        
        # Statistics box
        n_events = len(m1)
        textstr = f'Events: {n_events:,}\nMean $m_1$: {np.mean(m1):.3f}\nMean $m_2$: {np.mean(m2):.3f}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        return fig, ax

    # ========== 1D Projection Plots ==========
    
    def plot_m1_projection(self, bins: int = 60, data_type: str = 'os', figsize: Tuple[float, float] = (10, 6)) -> Tuple[plt.Figure, plt.Axes]:
        """Plot m1 projection, optionally comparing OS and LS"""
        fig, ax = plt.subplots(figsize=figsize)
        
        if self.data is not None:
            if data_type.lower() == 'os':
                m1_os, _ = self.get_os_masses()
                ax.hist(m1_os, bins=bins, histtype='step', linewidth=2,
                       color='#2E86AB', label=f'OS ({len(m1_os)} events)', alpha=0.9)
            elif data_type.lower() == 'ls':
                m1_ls, _ = self.get_ls_masses()
                ax.hist(m1_ls, bins=bins, histtype='step', linewidth=2,
                       color='#A23B72', label=f'LS ({len(m1_ls)} events)', alpha=0.9)
            else:
                raise ValueError(f"Unknown data_type '{data_type}'. Use 'os' or 'ls'.")
        ax.set_xlabel(r'$m_1$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel('Events / bin', fontsize=13)
        ax.set_title('$m_1$ Projection', fontsize=14, fontweight='bold')
        ax.legend(frameon=True, shadow=True)
        ax.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_m2_projection(self, bins: int = 60, data_type: str = 'os',
                          figsize: Tuple[float, float] = (10, 6)) -> Tuple[plt.Figure, plt.Axes]:
        """Plot m2 projection"""
        fig, ax = plt.subplots(figsize=figsize)

        if self.data is not None:
            if data_type.lower() == 'os':
                _, m2 = self.get_os_masses()
            elif data_type.lower() == 'ls':
                _, m2 = self.get_ls_masses()
            else:
                raise ValueError(f"Unknown data_type '{data_type}'. Use 'os' or 'ls'.")

        ax.set_xlabel(r'$m_2$ [GeV/$c^2$]', fontsize=13)
        ax.set_ylabel('Events / bin', fontsize=13)
        ax.set_title('$m_2$ Projection', fontsize=14, fontweight='bold')
        ax.legend(frameon=True, shadow=True)
        ax.grid(alpha=0.3, linestyle='--')

        plt.tight_layout()
        return fig, ax
    
    # ========== Kinematic Plots ==========

    @staticmethod
    def constrain_angle(angle: np.ndarray, minimum: float = 0.0, harmonic: int = 1) -> np.ndarray:
        """
        Constrain angle(s) to [minimum, minimum + 2π/harmonic).

        Mirrors the behaviour of o2::math_utils::detail::constrainAngle used in O2
        for Δφ calculations (see TVector2::Phi_0_2pi in ROOT).

        Parameters:
            angle:    Angle(s) in radians (scalar or numpy array).
            minimum:  Lower bound of the output range (default: 0 → [0, 2π)).
            harmonic: Harmonic number; period becomes 2π/harmonic (default: 1).

        Returns:
            Angle(s) within [minimum, minimum + 2π/harmonic).
        """
        period = 2.0 * np.pi / harmonic
        angle = np.asarray(angle, dtype=float).copy()
        # Vectorised equivalent of the while-loop pair in the C++ original
        angle -= period * np.floor((angle - minimum) / period)
        return angle

    # We only constrain the mass region for m1
    def constrain_mass_region(self, mass_min, mass_max, data_type: str = 'os'):
        """
        Constrain the data to a specific mass region for both m1 and m2.

        Parameters:
            mass_min: Minimum mass value (inclusive).
            mass_max: Maximum mass value (inclusive).
        """
        if self.data is not None:
            if data_type.lower() == 'os':
                m1_os, _ = self.get_os_masses()
                mask = (m1_os >= mass_min) & (m1_os <= mass_max)
            elif data_type.lower() == 'ls':
                m1_ls, _ = self.get_ls_masses()
                mask = (m1_ls >= mass_min) & (m1_ls <= mass_max)
            else:
                raise ValueError(f"Unknown data_type '{data_type}'. Use 'os' or 'ls'.")
            self.data = {key: val[mask] for key, val in self.data.items()}

    def add_deltas(self, data_type: str = 'os'):
        data = self.data
        if data is None:
            raise ValueError(f"No {data_type.upper()} data available")

        # --- Δy ---
        data['deltaYnoAbs'] = data['y1'] - data['y2']
        data['deltaY'] = np.abs(data['deltaYnoAbs'])

        # --- Δφ ---
        # Raw difference, then constrain to [-π/2, 3π/2) following the O2
        # convention (minimum = -π/2).  This places back-to-back pairs at π
        # and collinear pairs at 0, with no discontinuity in the plotted range.
        raw_dphi = data['phi1'] - data['phi2']
        data['deltaPhinoAbs'] = self.constrain_angle(raw_dphi, minimum=-np.pi / 2.0)
        data['deltaPhi'] = data['deltaPhinoAbs']   # already in [−π/2, 3π/2); take abs if desired
    
    def plot_kinematics(self, data_type: str = 'os', figsize: Tuple[float, float] = (14, 10)):
        """Plot kinematic distributions"""
        data = self.data
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

        fig, axes = plt.subplots(2, 2, figsize=figsize)
    
        # Define common bin edges for fair comparison
        pt_bins = np.linspace(0.0, 36.0, 91)
        mass_bins = np.linspace(0.0, 50.0, 501)
        deltay_bins = 50
        deltaphi_bins = 50
    
        # Pair pT
        if self.data is not None:
            axes[0, 0].hist(self.data['pairPt'], bins=pt_bins, histtype='step', linewidth=2,
                        color='#2E86AB', alpha=0.9, label='OS')
        axes[0, 0].set_xlabel(r'$p_{T,pair}$ [GeV/$c$]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_yscale('log')
        axes[0, 0].set_xlim(0.0, 36.0)
        axes[0, 0].set_title('Pair $p_T$', fontweight='bold')
        axes[0, 0].legend(frameon=True, shadow=True)
        axes[0, 0].grid(alpha=0.3)
    
        # Pair mass
        if self.data is not None:
            axes[0, 1].hist(self.data['pairMass'], bins=mass_bins, histtype='step', linewidth=2,
                       color='#C73E1D', alpha=0.9, label='OS')
        axes[0, 1].set_xlabel(r'$M_{pair}$ [GeV/$c^2$]')
        axes[0, 1].set_xlim(3.0, 12.0)
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_title('Pair Mass', fontweight='bold')
        axes[0, 1].legend(frameon=True, shadow=True)
        axes[0, 1].grid(alpha=0.3)
    
        # deltaY
        if self.data is not None and 'deltaY' in self.data:
            axes[1, 0].hist(self.data['deltaY'], bins=deltay_bins, histtype='step', linewidth=2,
                       color='#A23B72', alpha=0.9, label='OS')
        axes[1, 0].set_xlabel(r'$\Delta y$')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_title('$\Delta y$', fontweight='bold')
        axes[1, 0].legend(frameon=True, shadow=True)
        axes[1, 0].grid(alpha=0.3)
    
        # deltaPhi
        if self.data is not None and 'deltaPhi' in self.data:
            axes[1, 1].hist(self.data['deltaPhi'], bins=deltaphi_bins, histtype='step', linewidth=2,
                       color='#F18F01', alpha=0.9, label='OS')
        axes[1, 1].set_xlabel(r'$\Delta\phi$')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title('$\Delta\phi$', fontweight='bold')
        axes[1, 1].legend(frameon=True, shadow=True)
        axes[1, 1].grid(alpha=0.3)
    
        plt.suptitle('OS vs LS Pair Kinematics Comparison' if compare else 'Pair Kinematics', 
                fontsize=16, fontweight='bold')
        plt.tight_layout()
        return fig, axes

    def plot_2d_delta_eta_phi(self, bins=(50, 50)):
        """
        ROOT LEGO-style 2D histogram of Δη vs Δφ.
        """


        if self.data is None:
            raise ValueError("No data available")

        self.add_deltas("os")

        delta_eta = self.data["deltaYnoAbs"]
        delta_phi = self.data["deltaPhi"]

        # --- determine axis ranges ---
        eta_min, eta_max = np.min(delta_eta), np.max(delta_eta)
        phi_min, phi_max = np.min(delta_phi), np.max(delta_phi)

        # --- ROOT histogram (TH2F) ---
        h = ROOT.TH2F(
            "h_delta_eta_phi",
            "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi",
            bins[0], eta_min, eta_max,
            bins[1], phi_min, phi_max
        )

        # --- fill histogram ---
        for de, dphi in zip(delta_eta, delta_phi):
            h.Fill(de, dphi)

        # --- canvas ---
        c = ROOT.TCanvas("c", "LEGO Plot", 800, 600)

        # --- style (optional but recommended) ---
        ROOT.gStyle.SetOptStat(0)

        # --- LEGO-style draw ---
        h.Draw("LEGO2")   # main LEGO mode

        # Alternative options:
        # h.Draw("LEGO")   # simpler blocks
        # h.Draw("COLZ")   # heatmap style
        # h.Draw("SURF2")  # smooth surface

        # --- axis labels styling (optional fine-tuning) ---
        h.GetXaxis().SetTitle("#Delta#eta")
        h.GetYaxis().SetTitle("#Delta#phi")
        h.GetZaxis().SetTitle("Events")

        c.Update()

        return c, h
    
    # ========== Convenience Methods for Quick Display ==========
    
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


# ========== Utility Functions ==========

def quick_look(file: Optional[str] = None):
    """
    Quick visualization of all data
    
    Usage:
        quick_look('os_pairs.root', 'ls_pairs.root')
    """
    viz = CorrelationFunctions(file)
    viz.show_2d()
    viz.show_projections()
    viz.show_correlation()
    viz.show_kinematics()
    viz.show_pair_kinematics()
    return viz


def main():
    """Command line interface (kept for backward compatibility)"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Correlation functions')
    parser.add_argument('--file', type=str,
                       default='os_pairs.root',
                       help='OS data file')
    parser.add_argument('--output-dir', type=str,
                       default='plots',
                       help='Output directory for plots')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Raw data visualizer")
    print("=" * 60)
    
    corr = CorrelationFunctions(args.os_file, args.ls_file)
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == '__main__':
    main()