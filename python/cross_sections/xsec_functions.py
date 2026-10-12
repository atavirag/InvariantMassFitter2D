import ROOT
import math

def cross_section(N_corr, dN_corr, f_prompt, BR, lumi, delta_y):
    """
    Input parameters:
        N_corr: number of events corrected by efficiencies
        f_prompt: prompt fraction
        BR: branching ratio of D⁰ --> K⁺pi⁻
        Lumi: integrated luminosity corrected by BC cuts
    Output:
        cross section, xsec, with statistical error, xsec_stat_err
    """

    if lumi == 0:
        raise ValueError("Luminosity values must be non-zero to avoid division by zero.")

    xsec = N_corr * f_prompt * f_prompt / (BR * BR * lumi * delta_y)
    xsec = xsec / 10**6 # Change to micro barn^(-1)
    xsec_stat_err = (dN_corr / N_corr) * xsec
    return [xsec, xsec_stat_err]

def os_over_ls(xsec_os, xsec_ls):
    """
    Input parameters:
        xsec_os: opposite-sign cross section [xsec, xsec_stat, xsec_sys]
        xsec_ls: like-sign cross section [xsec, xsec_stat, xsec_sys]
    Output:
        [os_over_ls, os_over_ls_stat, os_over_ls_sys]
    """
    os_over_ls = xsec_os[0] / xsec_ls[0]
    os_over_ls_stat = os_over_ls * math.sqrt((xsec_os[1]/xsec_os[0])**2 + (xsec_ls[1]/xsec_ls[0])**2)
    os_over_ls_sys = os_over_ls * math.sqrt((xsec_os[2]/xsec_os[0])**2 + (xsec_ls[2]/xsec_ls[0])**2)

    return [os_over_ls, os_over_ls_stat, os_over_ls_sys]

def cross_section_raw(N_raw, dN_raw, eff, deff, f_prompt, BR, lumi, delta_y):
    """
    Input parameters:
        N_raw: number of raw events
        eff: efficiency in the pT interval
        f_prompt: prompt fraction
        BR: branching ratio of D⁰ --> K⁺pi⁻
        Lumi: integrated luminosity corrected by BC cuts
    Output: 
        cross section, xsec, with statistical error, xsec_stat_err
    """

    if lumi == 0 or eff == 0:
        raise ValueError("Luminosity and efficiency values must be non-zero to avoid division by zero.")

    # Central value
    xsec = N_raw * f_prompt * f_prompt / (BR * BR * lumi * eff * eff * delta_y)
    xsec = xsec / 10**6 # Change to micro barn^(-1)

    # Statistical error
    xsec_stat_err = (dN_raw / N_raw + deff / eff) * xsec

    return [xsec, xsec_stat_err]


def sys_uncert(sys_list):
    """
    Calculates the total systematic uncertainty by summing in quadrature.

    Parameters:
        sys_list (list of float): individual systematic uncertainties

    Returns:
        float: total systematic uncertainty
    """
    total_sys_squared = sum(x**2 for x in sys_list)
    total_sys = math.sqrt(total_sys_squared)
    return total_sys

def not_sigma_eff(xsec_singleD, xsec_doubleD):
    """
    Calculate the 'not' effective cross section (sigma_eff) using:
        sigma_eff = (1 / 4) * sigma_single² / sigma_double

    Parameters:
        xsec_singleD: [value, stat_err, sys_err_up, sys_err_down]
        xsec_doubleD: [value, stat_err, sys_err]

    Returns:
        [not_sigma_eff, stat_uncertainty, sys_uncertainty_up, sys_uncertainty_down]
    """
    x_s, stat_s, sys_up_s, sys_down_s = xsec_singleD
    x_d, stat_d, sys_d = xsec_doubleD

    if x_s == 0 or x_d == 0:
        raise ValueError("Cross section values must be non-zero to avoid division by zero.")

    # Central value
    eff = 0.25 * x_s ** 2 / x_d

    # Relative errors
    rel_stat = (2 * stat_s / x_s) + (stat_d / x_d)
    rel_sys_up = (2 * sys_up_s / x_s) + (sys_d / x_d)
    rel_sys_down = (2 * sys_down_s / x_s) + (sys_d / x_d)

    # Uncertainties
    eff_stat = rel_stat * eff
    eff_sys_up = rel_sys_up * eff
    eff_sys_down = rel_sys_down * eff

    return [eff, eff_stat, eff_sys_up, eff_sys_down]

def xsec_with_uncerts(yield_value, yield_error, rel_sys, f_prompt, BR, lumi, delta_y):
    """
    Compute cross section and total systematic uncertainty.

    Returns:
        list: [xsec, stat_err, sys_err]
    """
    xsec = cross_section(yield_value, yield_error, f_prompt, BR, lumi, delta_y)
    abs_sys = [xsec[0] * rel for rel in rel_sys]
    sys_total = sys_uncert(abs_sys)
    xsec.append(sys_total)
    return xsec

def diff_xsec_with_uncerts(yield_value, yield_error, rel_sys, f_prompt, BR, lumi, delta_y, pt_range):
    """
    Compute differential cross section accounting for bin width and total systematic uncertainty.

    Returns:
        list: [xsec, stat_err, sys_err]
    """

    xsec_list = xsec_with_uncerts(yield_value, yield_error, rel_sys, f_prompt, BR, lumi, delta_y)
    pt_width = pt_range[1] - pt_range[0]
    xsec = [x / pt_width for x in xsec_list]

    return xsec