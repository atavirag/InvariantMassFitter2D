# compute_xsec.py

from .load_fitter_results import load_histograms_multi, load_single_value_histo
from . import xsec_functions as xs
from . import plot_xsecs as plx
import math


# Load histograms
f_prompt = load_single_value_histo('/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/raw_prompt_frac_pol0_swapped.root', "hFitParam")
lumi = load_single_value_histo('/home/andrea/MyMacros/InvariantMassFitter/test/Integrated_Lumi.root', "hCorrLumi")
data = load_histograms_multi(["hYields"])

# BR of D⁰ --> K⁺pi⁻
# Source: PDG
BR, dBR = 0.03950, 0.0156

delta_y = 1.4

# Single D0 cross section from ALICE paper at 13 TeV
# D0 cross section in 1.0 < pt < 24.0 = 700.95 +- 24.49 (stat) + 21.54 - 21.94 (sys)
# D0 cross section in 2.0 < pt < 24.0 = 447.95 +- 7.98 (stat) + 13.66 - 14.30 (sys)
single_D0_xsec_1_24 = [700.95, 24.49, 21.54, 21.94]
single_D0_xsec_2_24 = [447.95, 7.98, 13.66, 14.30]

single_D0_xsec_1_24_eta_07 = [x * 1.4 for x in single_D0_xsec_1_24]
single_D0_xsec_2_24_eta_07 = [x * 1.4 for x in single_D0_xsec_2_24]

# Systematic uncertainties
# relative systematic sources: 
#       [f_prompt, lumi, trigger_eff, BDT_sel, yield, ITS_TPC, single_track, track_tuner, BR]

# Access data
pt_bins = ["0_2", "2_4", "4_6", "6_12", "12_24"]
int_pt_bin = "0_100"
xsecs_os = []
xsecs_ls = []

xsec_os_sum = [0, 0, 0]
xsec_ls_sum = [0, 0, 0]

os_over_ls_all = []

int_rel_sys = [0.02, 0.10, 0.05, 0.06, 0.05, 0.08, 0.05, 0.02, dBR]
rel_sys = {
    'pt_bin': pt_bins,
    'uncert': [[0.02, 0.10, 0.05, 0.08, 0.05, 0.08, 0.05, 0.02, dBR],
               [0.02, 0.10, 0.05, 0.08, 0.06, 0.08, 0.05, 0.02, dBR],
               [0.02, 0.10, 0.05, 0.06, 0.05, 0.08, 0.05, 0.02, dBR],
               [0.02, 0.10, 0.05, 0.08, 0.07, 0.08, 0.05, 0.02, dBR],
               [0.02, 0.10, 0.05, 0.15, 0.11, 0.08, 0.05, 0.02, dBR]
            ]
    }

bdt_sel_uncert = [0.08, 0.08, 0.06, 0.08, 0.15]
yield_uncert = [0.05, 0.06, 0.05, 0.07, 0.11]


# for pt_bin in pt_bins:
#     # TODO: add efficiencies
#     """ # load yields
#     os_yield = data[pt_bin]["OS"]["hYields"]["values"]
#     ls_yield = data[pt_bin]["LS"]["hYields"]["values"]

#     os_yield_stat = data[pt_bin]["OS"]["hYields"]["errors"]
#     ls_yield_stat = data[pt_bin]["LS"]["hYields"]["errors"]

#     # Compute cross section with corrected yields
#     xsec_os_raw = xs.cross_section_raw(os_yield[0], os_yield_stat[0], f_prompt, BR, lumi, delta_y)
#     print("Bin ", pt_bin, "xsec_os:", xsec_os_raw) """

#     # load corrected yields
#     os_yield_corr = data[pt_bin]["OS"]["hYieldsCorr"]["values"][0]
#     ls_yield_corr = data[pt_bin]["LS"]["hYieldsCorr"]["values"][0]

#     os_yield_corr_stat = data[pt_bin]["OS"]["hYieldsCorr"]["errors"][0]
#     ls_yield_corr_stat = data[pt_bin]["LS"]["hYieldsCorr"]["errors"][0]

#     pt_bin_parts = pt_bin.split('_')  # splits into ['0', '100']
#     pt = [int(part) for part in pt_bin_parts]  # converts each to int
#     pt_bin_index = pt_bins.index(pt_bin)

#     xsec_os = xs.diff_xsec_with_uncerts(os_yield_corr, os_yield_corr_stat, rel_sys['uncert'][pt_bin_index], f_prompt, BR, lumi, delta_y, pt)
#     xsec_ls = xs.diff_xsec_with_uncerts(ls_yield_corr, ls_yield_corr_stat, rel_sys['uncert'][pt_bin_index], f_prompt, BR, lumi, delta_y, pt)

#     xsec_os_times_pt = xs.xsec_with_uncerts(os_yield_corr, os_yield_corr_stat, rel_sys['uncert'][pt_bin_index], f_prompt, BR, lumi, delta_y)
#     xsec_ls_times_pt = xs.xsec_with_uncerts(ls_yield_corr, ls_yield_corr_stat, rel_sys['uncert'][pt_bin_index], f_prompt, BR, lumi, delta_y)

#     os_over_ls = xs.os_over_ls(xsec_os, xsec_ls)

#     xsecs_os.append(xsec_os)
#     xsecs_ls.append(xsec_ls)

#     os_over_ls_all.append(os_over_ls)

#     xsec_os_sum[0] = xsec_os_sum[0] + xsec_os_times_pt[0]
#     xsec_ls_sum[0] = xsec_ls_sum[0] + xsec_ls_times_pt[0]

#     xsec_os_sum[1] = xsec_os_sum[1] + xsec_os_times_pt[1]
#     xsec_ls_sum[1] = xsec_ls_sum[1] + xsec_ls_times_pt[1]

#     xsec_os_sum[2] = (yield_uncert[pt_bin_index]**2 + bdt_sel_uncert[pt_bin_index]**2)*xsec_os_times_pt[0]**2 + xsec_os_sum[2]
#     xsec_ls_sum[2] = (yield_uncert[pt_bin_index]**2 + bdt_sel_uncert[pt_bin_index]**2)*xsec_ls_times_pt[0]**2 + xsec_ls_sum[2]

# # Add the uncertainties which are common accross pT ranges
# xsec_os_sum[2] += (0.02**2 + 0.1**2 + 0.05**2 + 0.08**2 + 0.05**2 + 0.02**2 + 0.01**2) * xsec_os_sum[0]**2
# xsec_ls_sum[2] += (0.02**2 + 0.1**2 + 0.05**2 + 0.08**2 + 0.05**2 + 0.02**2 + 0.01**2) * xsec_ls_sum[0]**2

# plx.plot_differential_xsecs(pt_bins, xsecs_os, "OS")
# plx.plot_differential_xsecs(pt_bins, xsecs_ls, "LS")
# plx.plot_differential_os_over_ls(pt_bins, os_over_ls_all)
# print("Plot saved")

# integrated cross section
# load corrected yields
int_os_yield_corr = data["OS"]["hYields"]["values"][0]
int_ls_yield_corr = data["LS"]["hYields"]["values"][0]

int_os_yield_corr_stat = data["OS"]["hYields"]["errors"][0]
int_ls_yield_corr_stat = data["LS"]["hYields"]["errors"][0]

int_xsec_os = xs.xsec_with_uncerts(int_os_yield_corr, int_os_yield_corr_stat, int_rel_sys, f_prompt, BR, lumi, delta_y)
int_xsec_ls = xs.xsec_with_uncerts(int_ls_yield_corr, int_ls_yield_corr_stat, int_rel_sys, f_prompt, BR, lumi, delta_y)

os_over_ls = xs.os_over_ls(int_xsec_os, int_xsec_ls)

print(f"Bin: {int_pt_bin} | xsec OS/LS: {os_over_ls[0]} +- {os_over_ls[1]} (stat.) +- {os_over_ls[2]} (sys.)")
print(f"       | xsec OS: {int_xsec_os[0]} +- {int_xsec_os[1]} (stat.) +- {int_xsec_os[2]} (sys.)")
print(f"       | xsec LS: {int_xsec_ls[0]} +- {int_xsec_ls[1]} (stat.) +- {int_xsec_ls[2]} (sys.)")

# print("")
# print(f"       | sum of differential xsec OS: {xsec_os_sum[0]} +- {xsec_os_sum[1]} (stat.) +- {math.sqrt(xsec_os_sum[2])} (sys.)" )
# print(f"       | sum of differential xsec LS: {xsec_ls_sum[0]} +- {xsec_ls_sum[1]} (stat.) +- {math.sqrt(xsec_ls_sum[2])} (sys.)")
plx.plot_int_cross_xsecs(int_xsec_os, int_xsec_ls)

# Compute the ratio of cross sections
r_os = xs.not_sigma_eff(single_D0_xsec_1_24_eta_07, int_xsec_os)
r_ls = xs.not_sigma_eff(single_D0_xsec_1_24_eta_07, int_xsec_ls)
# from mu b to mb

r_os = [x * 0.001 for x in r_os]
r_ls = [x * 0.001 for x in r_ls]

print("  ")
print(f"Bin: {int_pt_bin} | ratio OS: {r_os[0]} +- {r_os[1]} (stat.) +- {r_os[2]} (sys.)")
print(f"Bin: {int_pt_bin} | ratio LS: {r_ls[0]} +- {r_ls[1]} (stat.) +- {r_ls[2]} (sys.)")

plx.plot_sigma_eff(r_os, r_ls)
