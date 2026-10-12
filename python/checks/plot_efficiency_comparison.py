import ROOT
import matplotlib.pyplot as plt
import numpy as np

def get_efficiency_values(file_path, histo_name):
    """Extract efficiency values from a TEfficiency histogram"""
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open {file_path}")
        return None, None
    
    # Get the TEfficiency object
    eff_obj = root_file.Get(histo_name)
    if not eff_obj:
        print(f"Error: Cannot find {histo_name} in {file_path}")
        root_file.Close()
        return None, None
    
    # Get the total histogram to find bin edges
    total_hist = eff_obj.GetTotalHistogram()
    n_bins = total_hist.GetNbinsX()
    
    pt_values = []
    eff_values = []
    eff_errors_low = []
    eff_errors_high = []
    pt_errors_low = []
    pt_errors_high = []
    
    for i in range(1, n_bins + 1):
        pt = total_hist.GetBinCenter(i)
        pt_low = total_hist.GetBinWidth(i)/2
        pt_high = total_hist.GetBinWidth(i)/2
        eff = eff_obj.GetEfficiency(i)
        err_low = eff_obj.GetEfficiencyErrorLow(i)
        err_high = eff_obj.GetEfficiencyErrorUp(i)
        
        pt_values.append(pt)
        pt_errors_low.append(pt_low)
        pt_errors_high.append(pt_high)
        eff_values.append(eff)
        eff_errors_low.append(err_low)
        eff_errors_high.append(err_high)
    
    root_file.Close()
    
    return (np.array(pt_values), np.array(eff_values), 
            np.array(eff_errors_low), np.array(eff_errors_high),
            np.array(pt_errors_low), np.array(pt_errors_high))

def plot_efficiency_comparison():
    """Plot efficiency comparison between Prelim and final version"""
    
    # File paths
    prelim_file = "AccEffD0ToKPi_Prelim.root"
    final_file = "AccEffD0ToKPi.root"
    
    # Get efficiency data
    pt_prelim_prompt, eff_prelim_prompt, err_low_prelim_prompt, err_high_prelim_prompt, pt_low_prelim_prompt, pt_high_prelim_prompt = \
        get_efficiency_values(prelim_file, "TEfficiency_AccEffD0ToKPi_PrelimPrompt")
    
    pt_prelim_nonprompt, eff_prelim_nonprompt, err_low_prelim_nonprompt, err_high_prelim_nonprompt, pt_low_prelim_nonprompt, pt_high_prelim_nonprompt = \
        get_efficiency_values(prelim_file, "TEfficiency_AccEffD0ToKPi_PrelimNonprompt")
    
    pt_final_prompt, eff_final_prompt, err_low_final_prompt, err_high_final_prompt, pt_low_final_prompt, pt_high_final_prompt = \
        get_efficiency_values(final_file, "TEfficiency_AccEffD0ToKPiPrompt")
    
    pt_final_nonprompt, eff_final_nonprompt, err_low_final_nonprompt, err_high_final_nonprompt, pt_low_final_nonprompt, pt_high_final_nonprompt = \
        get_efficiency_values(final_file, "TEfficiency_AccEffD0ToKPiNonprompt")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot Prelim data
    prompt_prelim = ax.errorbar(pt_prelim_prompt, eff_prelim_prompt, 
                yerr=[err_low_prelim_prompt, err_high_prelim_prompt],
                xerr=[pt_low_prelim_prompt, pt_high_prelim_prompt],
                fmt='o', color='tomato', label='Prompt Efficiency (Prelim.)', 
                markersize=6, linewidth=2, capsize=0)
    
    nonprompt_prelim = ax.errorbar(pt_prelim_nonprompt, eff_prelim_nonprompt,
                yerr=[err_low_prelim_nonprompt, err_high_prelim_nonprompt],
                xerr=[pt_low_prelim_nonprompt, pt_high_prelim_nonprompt],
                fmt='^', color='dodgerblue', label='Non-Prompt Efficiency (Prelim.)',
                markersize=8, linewidth=2, capsize=0)
    
    # Plot Final data
    ax.errorbar(pt_final_prompt, eff_final_prompt,
                yerr=[err_low_final_prompt, err_high_final_prompt],
                xerr=[pt_low_final_prompt, pt_high_final_prompt],
                fmt='s', color='firebrick', label='Prompt Efficiency (This Thesis)',
                markersize=6, linewidth=2, alpha=0.7, capsize=0)
    
    ax.errorbar(pt_final_nonprompt, eff_final_nonprompt,
                yerr=[err_low_final_nonprompt, err_high_final_nonprompt],
                xerr=[pt_low_final_prompt, pt_high_final_prompt],
                fmt='D', color='darkblue', label='Non-Prompt Efficiency (This Thesis)',
                markersize=6, linewidth=2, alpha=0.7, capsize=0)
    
    #modify error bars for prelim.
    # 🔹 Modify the error bar lines (both x and y) separately:
    for err in prompt_prelim[2]:
        err.set_linestyle('-.')    # dashed error bars

    for err in nonprompt_prelim[2]:
        err.set_linestyle('-.')    # dashed error bars


    # Set log scale
    ax.set_yscale('log')
    
    # Labels and formatting
    ax.set_xlabel(r'$p_{\rm T}$ (GeV/$c$)', fontsize=16)
    ax.set_ylabel('Acceptance × Efficiency', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=16, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    
    # Add text annotations (similar to canvas)
    plt.figtext(0.63, 0.85, 'This Thesis', size=18)
    plt.figtext(0.63, 0.81, r'pp, $\sqrt{s} = 13.6$ TeV', size=16)
    plt.figtext(0.63, 0.77, r'D$^0$ $\rightarrow$ K$^- \pi^+$ and charge conj.', size=16)
    
    plt.tight_layout()
    plt.savefig('efficiency_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('efficiency_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\nPlot saved as:")
    print("  - efficiency_comparison.pdf")
    print("  - efficiency_comparison.png")

if __name__ == "__main__":
    # Set ROOT to batch mode
    ROOT.gROOT.SetBatch(True)
    
    plot_efficiency_comparison()