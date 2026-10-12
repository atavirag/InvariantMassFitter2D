# plot_xsecs.py

import matplotlib.pyplot as plt
import numpy as np

def plot_differential_xsecs(pt_bins, xsecs, pair_type):
    """
    Input:
        pt_bins: list of pt bins written as "low_high"
        xsecs: list of cross sections written as [xsec, stat_err, sys_err]
    Output:
        Plot of the xsecs put together with stat errs as lines and sys errs as boxes
    """

    # Extract edges as list of tuples of floats
    bin_edges = [tuple(map(float, b.split("_"))) for b in pt_bins]

    # Calculate bin centers
    bin_centers = [(low + high) / 2 for low, high in bin_edges]
    bin_errors = [(high - low) / 2 for low, high in bin_edges]

    # Extract values and uncertainties
    xsec_vals = [x[0] for x in xsecs]
    xsec_stat = [x[1] for x in xsecs]
    xsec_sys = [x[2] for x in xsecs]

    bin_centers = np.array(bin_centers)
    xsec_vals = np.array(xsec_vals)
    xsec_stat = np.array(xsec_stat)
    xsec_sys = np.array(xsec_sys)

    # Plot
    plt.figure(figsize=(8, 6))
    plt.yscale('log')
    #plt.xscale('log')

    # Systematic uncertainty band (boxes)
    plt.bar(
        bin_centers,
        2 * xsec_sys,          # height of the box: ±sys
        width=0.9,             # width of the box (can be adjusted)
        bottom=xsec_vals - xsec_sys,
        color="lightblue",
        alpha=0.7,
        label="Systematic uncertainty",
        align='center',
        edgecolor='none'
    )

    # Central values with statistical error bars
    plt.errorbar(
        bin_centers,
        xsec_vals,
        yerr=xsec_stat,
        xerr=bin_errors,
        fmt='o',
        color='blue',
        label='Cross section (stat. error)',
        capsize=0
    )

    plt.tick_params(axis='both', which='major', labelsize=14)  # for major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # optional, for minor ticks

    plt.xlabel(r"$p_{\mathrm{T}}^{\mathrm{DD}}$ (GeV/$c$)", fontsize=16)
    if pair_type == "LS":
        plt.ylabel(r"d$\sigma_{\mathrm{D}\mathrm{D}}$/d$p_{\mathrm{T}}^{\mathrm{DD}}$ ($\mu$b)", fontsize=18)
    elif pair_type == "OS":
        plt.ylabel(r"d$\sigma_{\mathrm{D}\bar{\mathrm{D}}}$/d$p_{\mathrm{T}}^{\mathrm{DD}}$ ($\mu$b)", fontsize=18)
    plt.title("")
    add_text(pair_type)
    plt.legend(frameon=False,fontsize=14, loc='lower left')
    plt.tick_params(axis='y', right=True)
    plt.tick_params(axis='x', top=True)
    plt.tick_params(direction='in')
    plt.tight_layout()

    plt.savefig(f"test/cross_sections/differential_xsec_{pair_type}.pdf")
    plt.close()

    return

def plot_differential_xsecs_both(pt_bins, xsecs_os, xsecs_ls):
    """
    Input:
        pt_bins: list of pt bins written as "low_high"
        xsecs_os: list of cross sections for OS written as [xsec, stat_err, sys_err]
        xsecs_ls: list of cross sections for LS written as [xsec, stat_err, sys_err]
    Output:
        Plot of both OS and LS xsecs with stat errs as lines and sys errs as boxes
    """
    bin_edges = [tuple(map(float, b.split("_"))) for b in pt_bins]
    bin_centers = np.array([(low + high) / 2 for low, high in bin_edges])
    bin_errors = np.array([(high - low) / 2 for low, high in bin_edges])

    def extract(xsecs):
        return (np.array([x[0] for x in xsecs]),
                np.array([x[1] for x in xsecs]),
                np.array([x[2] for x in xsecs]))

    os_vals, os_stat, os_sys = extract(xsecs_os)
    ls_vals, ls_stat, ls_sys = extract(xsecs_ls)

    # Slight horizontal offset so boxes don't overlap
    # offset = bin_errors * 0.05
    offset = 0

    plt.figure(figsize=(8, 6))
    plt.yscale('log')

    # OS systematic boxes
    #plt.bar(bin_centers - offset, 2 * os_sys, width=bin_errors * 0.9,
    #        bottom=os_vals - os_sys, color='lightcoral', alpha=0.7,
    #        align='center', edgecolor='none')
    # LS systematic boxes
    #plt.bar(bin_centers + offset, 2 * ls_sys, width=bin_errors * 0.9,
    #        bottom=ls_vals - ls_sys, color='lightblue', alpha=0.7,
    #        align='center', edgecolor='none')

    # OS central values + stat errors
    plt.errorbar(bin_centers - offset, os_vals, yerr=os_stat, xerr=bin_errors,
                 fmt='o', color='red', label=r'$\mathrm{D}^{0}\bar{\mathrm{D}}^{0}$', capsize=0)
    # LS central values + stat errors
    plt.errorbar(bin_centers + offset, ls_vals, yerr=ls_stat, xerr=bin_errors,
                 fmt='s', color='blue', label=r'$\mathrm{D}^{0}\mathrm{D}^{0} + \bar{\mathrm{D}}^{0}\bar{\mathrm{D}}^{0}$', capsize=0)

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=12)
    
    plt.xlabel(r"$p_{\mathrm{T}}^{\mathrm{DD}}$ (GeV/$c$)", fontsize=16)
    plt.ylabel(r"d$\sigma$/d$p_{\mathrm{T}}^{\mathrm{DD}}$ ($\mu$b)", fontsize=18)
    plt.title("")
    plt.legend(frameon=False, fontsize=14, loc='center right')
    plt.tick_params(axis='y', right=True)
    plt.minorticks_on()
    plt.tick_params(axis='x', which='minor', direction='in', top=True)
    plt.tick_params(direction='in')
    plt.tight_layout()

    ax = plt.gca()  # get current axes
    ax.text(0.55, 0.9, 'ALICE Work in progress', fontsize=16, color='black', transform=ax.transAxes)
    ax.text(0.55, 0.84, r'$\mathrm{D}^0 \rightarrow K^- \pi^+$ and charge conj.', fontsize=14, color='black', transform=ax.transAxes)
    ax.text(0.60, 0.78, r'pp, $\sqrt{\mathrm{s}}$ = 13.6 TeV', fontsize=14, color='black', transform=ax.transAxes)

    plt.savefig("test/cross_sections/differential_xsec_both.pdf")
    plt.close()

def plot_differential_os_over_ls(pt_bins, os_over_ls):
    """
    Input:
        pt_bins: list of pt bins written as "low_high"
        os_over_ls: list of ratios written as [os_over_ls, stat_err, sys_err]
    Output:
        Plot of the os_over_ls put together with stat errs as lines and sys errs as boxes
    """

    # Extract edges as list of tuples of floats
    bin_edges = [tuple(map(float, b.split("_"))) for b in pt_bins]

    # Calculate bin centers
    bin_centers = [(low + high) / 2 for low, high in bin_edges]
    bin_errors = [(high - low) / 2 for low, high in bin_edges]

    # Extract values and uncertainties
    xsec_vals = [x[0] for x in os_over_ls]
    xsec_stat = [x[1] for x in os_over_ls]
    xsec_sys = [x[2] for x in os_over_ls]

    bin_centers = np.array(bin_centers)
    xsec_vals = np.array(xsec_vals)
    xsec_stat = np.array(xsec_stat)
    xsec_sys = np.array(xsec_sys)

    # Plot
    plt.figure(figsize=(8, 6))

    # Systematic uncertainty band (boxes)
    plt.bar(
        bin_centers,
        2 * xsec_sys,          # height of the box: ±sys
        width=0.9,             # width of the box (can be adjusted)
        bottom=xsec_vals - xsec_sys,
        color="lightblue",
        alpha=0.7,
        label="Systematic uncertainty",
        align='center',
        edgecolor='none'
    )

    # Central values with statistical error bars
    plt.errorbar(
        bin_centers,
        xsec_vals,
        yerr=xsec_stat,
        xerr=bin_errors,
        fmt='o',
        color='blue',
        label='Cross section (stat. error)',
        capsize=0
    )

    plt.ylim(0.0, 10.0)
    plt.tick_params(axis='both', which='major', labelsize=14)  # for major ticks
    plt.tick_params(axis='both', which='minor', labelsize=12)  # optional, for minor ticks

    plt.xlabel(r"$p_{\mathrm{T}}^{\mathrm{DD}}$ (GeV/$c$)", fontsize=16)
    plt.ylabel(r"(d$\sigma_{\mathrm{D}\mathrm{D}}$/d$p_{\mathrm{T}}^{\mathrm{DD}}$) / (d$\sigma_{\mathrm{D}\bar{\mathrm{D}}}$/d$p_{\mathrm{T}}^{\mathrm{DD}}$)", fontsize=18)
    plt.title("")
    plt.legend(frameon=False,fontsize=14, loc='center left')
    plt.tick_params(axis='y', right=True)
    plt.tick_params(axis='x', top=True)
    plt.tick_params(direction='in')
    plt.tight_layout()

    ax = plt.gca()  # get current axes
    ax.text(0.05, 0.9, 'ALICE Work in progress', fontsize=16, color='black', transform=ax.transAxes)
    ax.text(0.05, 0.84, r'$\mathrm{D}^0 \rightarrow K^- \pi^+$ and charge conj.', fontsize=14, color='black', transform=ax.transAxes)
    ax.text(0.05, 0.78, r'pp, $\sqrt{\mathrm{s}}$ = 13.6 TeV', fontsize=14, color='black', transform=ax.transAxes)

    plt.savefig(f"test/cross_sections/differential_os_over_ls.pdf")
    plt.close()

    return

def add_text(pair_type):
    # Add text to the plot
    ax = plt.gca()  # get current axes

    # Use transform=ax.transAxes for relative positioning
    ax.text(0.6, 0.9, 'ALICE Work in progress', fontsize=16, color='black', transform=ax.transAxes)
    ax.text(0.57, 0.84, r'$\mathrm{D}^0 \rightarrow K^- \pi^+$ and charge conj.', fontsize=14, color='black', transform=ax.transAxes)
    ax.text(0.7, 0.78, r'pp, $\sqrt{\mathrm{s}}$ = 13.6 TeV', fontsize=14, color='black', transform=ax.transAxes)

    if pair_type == "LS":
        ax.text(0.8, 0.72, r'$\mathrm{D}^{0}\mathrm{D}^{0} + \bar{\mathrm{D}}^{0}\bar{\mathrm{D}}^{0}$', fontsize=14, color='black', transform=ax.transAxes)
    elif pair_type == "OS":
        ax.text(0.8, 0.72, r'$\mathrm{D}^{0}\bar{\mathrm{D}}^{0}$', fontsize=14, color='black', transform=ax.transAxes)

    return

def plot_int_cross_xsecs(xsec_os, xsec_ls):
    """
    Input:
        xsecs: cross section written as [xsec, stat_err, sys_err]
    Output:
        Plot of the xsecs put together with stat errs as lines and sys errs as boxes
    """

    # Extract cross section components
    xsec_os_val = xsec_os[0]
    xsec_os_stat = xsec_os[1]
    xsec_os_sys = xsec_os[2]

    xsec_ls_val = xsec_ls[0]
    xsec_ls_stat = xsec_ls[1]
    xsec_ls_sys = xsec_ls[2]

    # Plot cross sections
    plt.figure(figsize=(4, 4))

    plt.errorbar(xsec_os_val, 1, xerr=xsec_os_stat, marker='o', color='blue', label='OS')
    plt.errorbar(xsec_ls_val, 2, xerr=xsec_ls_stat, marker='o', color='blue', label='LS')

    plt.bar(
        x=xsec_os_val,                          # X position of the bar
        height=0.5,                               # Height of the box (arbitrary)
        width=2 * xsec_os_sys,                  # Total width (left + right uncertainty)
        bottom=0.75,                               # Bottom of the bar (optional)
        color='lightblue',
        alpha=0.7,
        label='Systematic uncertainty',
        align='center',
        edgecolor='none'
    )

    plt.bar(
        x=xsec_ls_val,                          # X position of the bar
        height=0.5,                               # Height of the box (arbitrary)
        width=2 * xsec_ls_sys,                  # Total width (left + right uncertainty)
        bottom=1.75,                               # Bottom of the bar (optional)
        color='lightblue',
        alpha=0.7,
        label='Systematic uncertainty',
        align='center',
        edgecolor='none'
    )

    plt.xlabel(r"$\sigma^{\mathrm{DD}}$ $(\mu$b)", fontsize=13)
    plt.xscale('log')
    plt.xlim(30, 350)
    plt.ylim(0, 4)
    plt.tick_params(axis='y', labelleft=False, left=False)
    plt.tight_layout()

    ax = plt.gca()  # get current axes

    # Use transform=ax.transAxes for relative positioning
    ax.text(0.4, 0.92, 'ALICE Work in progress', fontsize=13, color='black', transform=ax.transAxes)
    ax.text(0.33, 0.85, r'$\mathrm{D}^0 \rightarrow K^- \pi^+$ and charge conj.', fontsize=12, color='black', transform=ax.transAxes)
    ax.text(0.57, 0.78, r'pp, $\sqrt{\mathrm{s}}$ = 13.6 TeV', fontsize=12, color='black', transform=ax.transAxes)
    ax.text(0.57, 0.71, r'$1 < p_{\mathrm{T}} < 24$ GeV/$c$', fontsize=12, color='black', transform=ax.transAxes)

    ax.text(0.57, 0.48, r'$\mathrm{D}^{0}\mathrm{D}^{0} + \bar{\mathrm{D}}^{0}\bar{\mathrm{D}}^{0}$', fontsize=13, color='black', transform=ax.transAxes)
    ax.text(0.75, 0.22, r'$\mathrm{D}^{0}\bar{\mathrm{D}}^{0}$', fontsize=13, color='black', transform=ax.transAxes)


    plt.savefig("test/cross_sections/integrated_cross_sections.pdf")
    plt.close()

    return

def plot_sigma_eff(r_os, r_ls):
    """
    Input:
        xsecs: cross section written as [xsec, stat_err, sys_err]
    Output:
        Plot of the xsecs put together with stat errs as lines and sys errs as boxes
    """

    # Extract cross section components
    r_os_val = r_os[0]
    r_os_stat = r_os[1]
    r_os_sys = r_os[2]

    r_ls_val = r_ls[0]
    r_ls_stat = r_ls[1]
    r_ls_sys = r_ls[2]

    # Plot cross sections
    plt.figure(figsize=(4, 4))

    plt.errorbar(r_os_val, 1, xerr=r_os_stat, marker='o', color='blue', label='OS')
    plt.errorbar(r_ls_val, 2, xerr=r_ls_stat, marker='o', color='blue', label='LS')

    plt.bar(
        x=r_os_val,                          # X position of the bar
        height=0.5,                               # Height of the box (arbitrary)
        width=2 * r_os_sys,                  # Total width (left + right uncertainty)
        bottom=0.75,                               # Bottom of the bar (optional)
        color='lightblue',
        alpha=0.7,
        label='Systematic uncertainty',
        align='center',
        edgecolor='none'
    )

    plt.bar(
        x=r_ls_val,                          # X position of the bar
        height=0.5,                               # Height of the box (arbitrary)
        width=2 * r_ls_sys,                  # Total width (left + right uncertainty)
        bottom=1.75,                               # Bottom of the bar (optional)
        color='lightblue',
        alpha=0.7,
        label='Systematic uncertainty',
        align='center',
        edgecolor='none'
    )

    plt.xlabel(r"$\mathcal{R}^{\mathrm{DD}}$ (mb)", fontsize=13)
    plt.xscale('log')
    plt.xlim(0.5, 100)
    plt.ylim(0, 4)
    plt.tick_params(axis='y', labelleft=False, left=False)
    plt.tight_layout()

    ax = plt.gca()  # get current axes

    # Use transform=ax.transAxes for relative positioning
    ax.text(0.4, 0.92, 'ALICE Work in progress', fontsize=13, color='black', transform=ax.transAxes)
    ax.text(0.32, 0.85, r'$\mathrm{D}^0 \rightarrow K^- \pi^+$ and charge conj.', fontsize=12, color='black', transform=ax.transAxes)
    ax.text(0.55, 0.78, r'pp, $\sqrt{\mathrm{s}}$ = 13.6 TeV', fontsize=12, color='black', transform=ax.transAxes)
    ax.text(0.57, 0.71, r'$1 < p_{\mathrm{T}} < 24$ GeV/$c$', fontsize=12, color='black', transform=ax.transAxes)

    ax.text(0.57, 0.48, r'$\mathrm{D}^{0}\mathrm{D}^{0} + \bar{\mathrm{D}}^{0}\bar{\mathrm{D}}^{0}$', fontsize=13, color='black', transform=ax.transAxes)
    ax.text(0.57, 0.22, r'$\mathrm{D}^{0}\bar{\mathrm{D}}^{0}$', fontsize=13, color='black', transform=ax.transAxes)

    plt.savefig("test/cross_sections/ratio_cross_sections.pdf")
    print("Figure saved as test/cross_sections/ratio_cross_sections.pdf")
    plt.close()

    return