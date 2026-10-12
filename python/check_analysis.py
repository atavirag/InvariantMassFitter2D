#!/usr/bin/env python3
"""
check_analysis.py - Diagnostic tool for D meson pair analysis

Quickly check data integrity, fit quality, and plot examples
"""

import sys
from pathlib import Path
import argparse

# Color codes for terminal
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_section(title):
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{title:^60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}\n")

def check_mark(condition):
    return f"{Colors.OKGREEN}✓{Colors.ENDC}" if condition else f"{Colors.FAIL}✗{Colors.ENDC}"

def check_files(os_file, ls_file, workspace):
    """Check if required files exist"""
    print_section("File Check")
    
    files = {
        "OS data": os_file,
        "LS data": ls_file,
        "Workspace": workspace
    }
    
    all_exist = True
    for name, path in files.items():
        exists = Path(path).exists()
        size = Path(path).stat().st_size / (1024**2) if exists else 0
        
        print(f"{check_mark(exists)} {name:15} : {path}")
        if exists:
            print(f"    Size: {size:.2f} MB")
        all_exist = all_exist and exists
    
    return all_exist

def check_dependencies():
    """Check Python dependencies"""
    print_section("Python Dependencies")
    
    deps = {
        "NumPy": "numpy",
        "Matplotlib": "matplotlib",
        "Seaborn": "seaborn",
        "uproot": "uproot",
        "Awkward": "awkward",
        "SciPy": "scipy"
    }
    
    all_installed = True
    for name, module in deps.items():
        try:
            mod = __import__(module)
            version = getattr(mod, '__version__', 'unknown')
            print(f"{check_mark(True)} {name:15} : v{version}")
        except ImportError:
            print(f"{check_mark(False)} {name:15} : NOT INSTALLED")
            all_installed = False
    
    # Check PyROOT separately
    print(f"\n{Colors.BOLD}Optional (for advanced features):{Colors.ENDC}")
    try:
        import ROOT
        print(f"{check_mark(True)} PyROOT         : v{ROOT.__version__}")
    except ImportError:
        print(f"{check_mark(False)} PyROOT         : NOT INSTALLED (optional)")
    
    return all_installed

def analyze_root_file(filename, tree_name="pairs"):
    """Analyze ROOT file contents"""
    try:
        import uproot
    except ImportError:
        print(f"{Colors.WARNING}uproot not available, skipping analysis{Colors.ENDC}")
        return None
    
    try:
        with uproot.open(filename) as f:
            tree = f[tree_name]
            
            print(f"\n{Colors.OKBLUE}Tree: {tree_name}{Colors.ENDC}")
            print(f"  Entries: {tree.num_entries:,}")
            print(f"  Branches: {len(tree.keys())}")
            
            # Show branch info
            print(f"\n{Colors.BOLD}  Branch Information:{Colors.ENDC}")
            for branch in sorted(tree.keys())[:10]:  # Show first 10
                print(f"    - {branch}")
            
            if len(tree.keys()) > 10:
                print(f"    ... and {len(tree.keys()) - 10} more")
            
            # Get some data
            data = tree.arrays(['m1', 'm2', 'pt1', 'pt2', 'pairType'], 
                              library='np', entry_stop=1000)
            
            print(f"\n{Colors.BOLD}  Quick Statistics (first 1000 events):{Colors.ENDC}")
            print(f"    m1:  mean={data['m1'].mean():.4f}, std={data['m1'].std():.4f}")
            print(f"    m2:  mean={data['m2'].mean():.4f}, std={data['m2'].std():.4f}")
            print(f"    pt1: mean={data['pt1'].mean():.2f}, std={data['pt1'].std():.2f}")
            print(f"    pt2: mean={data['pt2'].mean():.2f}, std={data['pt2'].std():.2f}")
            
            # Check pairType distribution
            pair_types = data['pairType']
            dd = sum((pair_types & (1 << 0)) != 0)
            dbdb = sum((pair_types & (1 << 1)) != 0)
            ddb = sum((pair_types & (1 << 2)) != 0)
            dbd = sum((pair_types & (1 << 3)) != 0)
            
            print(f"\n{Colors.BOLD}  Pair Type Distribution:{Colors.ENDC}")
            print(f"    DD:       {dd:6d}")
            print(f"    DbarDbar: {dbdb:6d}")
            print(f"    DDbar:    {ddb:6d}")
            print(f"    DbarD:    {dbd:6d}")
            
            return True
            
    except Exception as e:
        print(f"{Colors.FAIL}Error analyzing file: {e}{Colors.ENDC}")
        return False

def check_workspace(workspace_file):
    """Check RooFit workspace"""
    print_section("Workspace Analysis")
    
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)  # Don't show graphics
        
        f = ROOT.TFile.Open(str(workspace_file))
        if not f or f.IsZombie():
            print(f"{Colors.FAIL}Cannot open workspace file{Colors.ENDC}")
            return False
        
        ws = f.Get("w")
        if not ws:
            print(f"{Colors.FAIL}Workspace 'w' not found{Colors.ENDC}")
            f.Close()
            return False
        
        print(f"{Colors.OKGREEN}Workspace loaded successfully{Colors.ENDC}\n")
        
        # Check for key components
        components = {
            "Variables": ["m1", "m2"],
            "PDFs": ["model"],
            "Datasets": ["osData", "lsData"],
            "Parameters": ["nSigSig", "nSigBkg", "nBkgSig", "nBkgBkg"]
        }
        
        for comp_type, items in components.items():
            print(f"{Colors.BOLD}{comp_type}:{Colors.ENDC}")
            for item in items:
                if comp_type == "Variables":
                    obj = ws.var(item)
                elif comp_type == "PDFs":
                    obj = ws.pdf(item)
                elif comp_type == "Datasets":
                    obj = ws.data(item)
                elif comp_type == "Parameters":
                    obj = ws.var(item)
                
                if obj:
                    if comp_type == "Datasets":
                        print(f"  {check_mark(True)} {item:15} : {obj.numEntries()} entries")
                    elif comp_type == "Parameters":
                        print(f"  {check_mark(True)} {item:15} : {obj.getVal():.2f} ± {obj.getError():.2f}")
                    else:
                        print(f"  {check_mark(True)} {item}")
                else:
                    print(f"  {check_mark(False)} {item}")
        
        # Check for fit results
        print(f"\n{Colors.BOLD}Fit Results:{Colors.ENDC}")
        for result_name in ["osFitResult", "lsFitResult"]:
            result = ws.obj(result_name)
            if result:
                print(f"  {check_mark(True)} {result_name}")
                print(f"      Status: {result.status()}")
                print(f"      Cov Quality: {result.covQual()}")
                print(f"      EDM: {result.edm():.2e}")
            else:
                print(f"  {check_mark(False)} {result_name}")
        
        f.Close()
        return True
        
    except ImportError:
        print(f"{Colors.WARNING}PyROOT not available, cannot check workspace{Colors.ENDC}")
        return False
    except Exception as e:
        print(f"{Colors.FAIL}Error: {e}{Colors.ENDC}")
        return False

def quick_plot(os_file, ls_file, output='quick_check.png'):
    """Create a quick diagnostic plot"""
    print_section("Quick Diagnostic Plot")
    
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        import uproot
    except ImportError as e:
        print(f"{Colors.WARNING}Required library not available: {e}{Colors.ENDC}")
        return False
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Load OS data
        with uproot.open(os_file) as f:
            tree = f['pairs']
            os_data = tree.arrays(['m1', 'm2', 'mbar1', 'mbar2', 'pairType'], 
                                  library='np')
        
        # Load LS data
        with uproot.open(ls_file) as f:
            tree = f['pairs']
            ls_data = tree.arrays(['m1', 'm2', 'mbar1', 'mbar2', 'pairType'], 
                                  library='np')
        
        # Helper function to get correct masses
        def get_masses(data, is_os=True):
            pairType = data['pairType']
            m1 = np.zeros(len(pairType))
            m2 = np.zeros(len(pairType))
            
            if is_os:
                mask_ddb = (pairType & (1 << 2)) != 0
                m1[mask_ddb] = data['m1'][mask_ddb]
                m2[mask_ddb] = data['mbar2'][mask_ddb]
                
                mask_dbd = (pairType & (1 << 3)) != 0
                m1[mask_dbd] = data['mbar1'][mask_dbd]
                m2[mask_dbd] = data['m2'][mask_dbd]
            else:
                mask_dd = (pairType & (1 << 0)) != 0
                m1[mask_dd] = data['m1'][mask_dd]
                m2[mask_dd] = data['m2'][mask_dd]
                
                mask_dbdb = (pairType & (1 << 1)) != 0
                m1[mask_dbdb] = data['mbar1'][mask_dbdb]
                m2[mask_dbdb] = data['mbar2'][mask_dbdb]
            
            return m1[m1 > 0], m2[m2 > 0]  # Filter zeros
        
        m1_os, m2_os = get_masses(os_data, is_os=True)
        m1_ls, m2_ls = get_masses(ls_data, is_os=False)
        
        # Plot 1: OS m1
        axes[0, 0].hist(m1_os, bins=50, histtype='step', linewidth=2, color='blue')
        axes[0, 0].set_xlabel('m1 [GeV/c²]')
        axes[0, 0].set_ylabel('Events')
        axes[0, 0].set_title(f'OS: m1 ({len(m1_os)} events)')
        axes[0, 0].grid(alpha=0.3)
        
        # Plot 2: OS m2
        axes[0, 1].hist(m2_os, bins=50, histtype='step', linewidth=2, color='red')
        axes[0, 1].set_xlabel('m2 [GeV/c²]')
        axes[0, 1].set_ylabel('Events')
        axes[0, 1].set_title(f'OS: m2 ({len(m2_os)} events)')
        axes[0, 1].grid(alpha=0.3)
        
        # Plot 3: LS m1
        axes[1, 0].hist(m1_ls, bins=50, histtype='step', linewidth=2, color='green')
        axes[1, 0].set_xlabel('m1 [GeV/c²]')
        axes[1, 0].set_ylabel('Events')
        axes[1, 0].set_title(f'LS: m1 ({len(m1_ls)} events)')
        axes[1, 0].grid(alpha=0.3)
        
        # Plot 4: LS m2
        axes[1, 1].hist(m2_ls, bins=50, histtype='step', linewidth=2, color='orange')
        axes[1, 1].set_xlabel('m2 [GeV/c²]')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title(f'LS: m2 ({len(m2_ls)} events)')
        axes[1, 1].grid(alpha=0.3)
        
        plt.suptitle('Quick Diagnostic Plot', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output, dpi=150, bbox_inches='tight')
        print(f"{Colors.OKGREEN}✓ Plot saved: {output}{Colors.ENDC}")
        
        return True
        
    except Exception as e:
        print(f"{Colors.FAIL}Error creating plot: {e}{Colors.ENDC}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Diagnostic tool for D meson pair analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python check_analysis.py --all
  python check_analysis.py --files-only
  python check_analysis.py --plot-only
        """)
    
    parser.add_argument('--os-file', type=str, default='os_pairs.root',
                       help='OS data file')
    parser.add_argument('--ls-file', type=str, default='ls_pairs.root',
                       help='LS data file')
    parser.add_argument('--workspace', type=str, default='fit_workspace.root',
                       help='Workspace file')
    parser.add_argument('--all', action='store_true',
                       help='Run all checks (default)')
    parser.add_argument('--files-only', action='store_true',
                       help='Only check files')
    parser.add_argument('--plot-only', action='store_true',
                       help='Only create diagnostic plot')
    parser.add_argument('--output', type=str, default='diagnostic_plot.png',
                       help='Output plot filename')
    
    args = parser.parse_args()
    
    # Default to all checks if no specific flag
    if not (args.files_only or args.plot_only):
        args.all = True
    
    print(f"{Colors.HEADER}{Colors.BOLD}")
    print("╔════════════════════════════════════════════════════════════╗")
    print("║        D Meson Pair Analysis - Diagnostic Tool            ║")
    print("╚════════════════════════════════════════════════════════════╝")
    print(f"{Colors.ENDC}")
    
    results = {}
    
    # Check files
    if args.all or args.files_only:
        results['files'] = check_files(args.os_file, args.ls_file, args.workspace)
    
    # Check dependencies
    if args.all:
        results['deps'] = check_dependencies()
    
    # Analyze data files
    if args.all and results.get('files', False):
        print_section("OS Data Analysis")
        results['os_analysis'] = analyze_root_file(args.os_file)
        
        print_section("LS Data Analysis")
        results['ls_analysis'] = analyze_root_file(args.ls_file)
    
    # Check workspace
    if args.all and Path(args.workspace).exists():
        results['workspace'] = check_workspace(args.workspace)
    
    # Create diagnostic plot
    if (args.all or args.plot_only) and Path(args.os_file).exists() and Path(args.ls_file).exists():
        results['plot'] = quick_plot(args.os_file, args.ls_file, args.output)
    
    # Summary
    print_section("Summary")
    
    total_checks = len(results)
    passed_checks = sum(1 for v in results.values() if v)
    
    print(f"Total checks: {total_checks}")
    print(f"Passed: {Colors.OKGREEN}{passed_checks}{Colors.ENDC}")
    print(f"Failed: {Colors.FAIL}{total_checks - passed_checks}{Colors.ENDC}")
    
    if passed_checks == total_checks:
        print(f"\n{Colors.OKGREEN}{Colors.BOLD}✓ All checks passed! Analysis looks good.{Colors.ENDC}")
        return 0
    else:
        print(f"\n{Colors.WARNING}{Colors.BOLD}⚠ Some checks failed. Review the output above.{Colors.ENDC}")
        return 1

if __name__ == '__main__':
    sys.exit(main())