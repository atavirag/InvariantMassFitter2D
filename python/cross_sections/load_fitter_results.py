# load_fitter_results.py

import uproot
import glob
import re
from collections import defaultdict
import os

#def load_histograms_multi(hist_names, file_pattern="/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/pt_differential_results/correlations_*.root"):
def load_histograms_multi(hist_names, corrected=True, file_pattern_corrected="/home/andrea/MyMacros/InvariantMassFitter/data/data_2025/fit_results_exponential_corrected_mass_1.73_2.05_eta_0.7_pairPt_*.root", file_pattern_uncorrected="/home/andrea/MyMacros/InvariantMassFitter/data/fit_results_exponential_mass_1.73_2.05_eta_0.7_pairPt_*.root"):
    """
    Load OS/LS histogram pairs from ROOT files with structured names.

    Parameters:
        hist_names (list of str): List of histogram names to load (e.g., ['hYields_os', 'hYields_ls'])
        corrected (bool): If True, load from corrected files; if False, load from uncorrected files
        file_pattern_corrected (str): Glob pattern to find corrected ROOT files
        file_pattern_uncorrected (str): Glob pattern to find uncorrected ROOT files

    Returns:
        dict: Nested dictionary structure:
            {
              '0_2': {
                  'OS': {
                      'hYields': {'values': ..., 'edges': ..., 'errors': ...},
                  },
                  'LS': {
                      'hYields': {'values': ..., 'edges': ..., 'errors': ...},
                  }
              },
              ...
            }
    """
    file_pattern = file_pattern_corrected if corrected else file_pattern_uncorrected
    files = glob.glob(file_pattern)
    
    # Pattern to extract pt bin from filename (e.g., fit_corrected_results_eta_0.7_pairPt_0_2.root)
    pattern = re.compile(r"fit_results_exponential_(?:corrected_)?mass_1\.73_2\.05_eta_0\.7_pairPt_(\d+)_(\d+)\.root")
    
    hist_data = {}

    for fname in files:
        match = pattern.search(fname)
        if not match:
            print(f"⚠️  Filename doesn't match expected pattern: {fname}")
            continue
        
        # Extract pt_bin from filename
        pt_min, pt_max = match.groups()
        pt_bin = f"{pt_min}_{pt_max}"
        print(f"✅ Loading pt bin: {pt_bin} from {os.path.basename(fname)}")

        try:
            with uproot.open(fname) as root_file:
                os_hists = {}
                ls_hists = {}

                # Load OS and LS histograms
                for hist_name in hist_names:
                    os_name = f"{hist_name}_os" if not hist_name.endswith('_os') else hist_name
                    ls_name = f"{hist_name}_ls" if not hist_name.endswith('_ls') else hist_name
                    
                    if os_name not in root_file:
                        raise KeyError(f"Histogram '{os_name}' not found in file for pt bin {pt_bin}")
                    if ls_name not in root_file:
                        raise KeyError(f"Histogram '{ls_name}' not found in file for pt bin {pt_bin}")

                    os_hist = root_file[os_name]
                    ls_hist = root_file[ls_name]

                    # Extract base name without _os/_ls suffix for storage
                    base_name = hist_name.replace('_os', '').replace('_ls', '')
                    
                    os_hists[base_name] = {
                        "values": os_hist.values(),
                        "errors": os_hist.errors(),
                        "edges": os_hist.axis(0).edges()
                    }
                    ls_hists[base_name] = {
                        "values": ls_hist.values(),
                        "errors": ls_hist.errors(),
                        "edges": ls_hist.axis(0).edges()
                    }
                    print(f"  Loaded: {os_name}")

                # Store in the nested structure with pt_bin as key
                hist_data[pt_bin] = {
                    "OS": os_hists,
                    "LS": ls_hists
                }

        except Exception as e:
            print(f"❌ Error loading pt bin {pt_bin}: {e}")

    return hist_data


def load_single_value_histo(file_path, hist_name):
    """
    Loads the value from a single-bin histogram in a ROOT file.

    Parameters:
        file_path (str): Path to the ROOT file
        hist_name (str): Name of the histogram inside the file

    Returns:
        float: The value of the only bin in the histogram
    """
    with uproot.open(file_path) as f:
        if hist_name not in f:
            raise KeyError(f"Histogram '{hist_name}' not found in file: {file_path}")
        
        hist = f[hist_name]
        values = hist.values()

        if len(values) != 1:
            raise ValueError(f"Histogram '{hist_name}' is not a single-bin histogram (has {len(values)} bins).")

        return values[0]


def load_full_histo(file_path, hist_name):
    """
    Loads values and errors from a multiple-bin histogram in a ROOT file.

    Parameters:
        file_path (str): Path to the ROOT file
        hist_name (str): Name of the histogram inside the file

    Returns:
        tuple: (values, errors) arrays from the histogram
    """
    with uproot.open(file_path) as f:
        if hist_name not in f:
            raise KeyError(f"Histogram '{hist_name}' not found in file: {file_path}")
        
        hist = f[hist_name]
        values = hist.values()
        errors = hist.errors()

        return values, errors