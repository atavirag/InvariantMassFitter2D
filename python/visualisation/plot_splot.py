#!/usr/bin/env python3
"""
plot_splot.py - visualise the raw sPlot distributions written by SPlotter.

For every variable it produces:
  1. <var>_components_os.pdf   : the four OS components overlaid
  2. <var>_components_ls.pdf   : the four LS components overlaid
  3. <var>_ss_os_vs_ls.pdf     : SgnSgn OS vs LS (+ OS/LS ratio pad)

Usage:
  python plot_splot.py ../../data/splot/splot_raw_exponential_mass_1.73_2.05_eta_0.7_pairPt_0_100.root \
        --outdir figures/splot [--variables pt1 deltaPhi] [--no-norm] [--logy]
"""
import argparse
import os
import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

SPECIES = ["SgnSgn", "SgnBkg", "BkgSgn", "BkgBkg"]
LABELS = {
    "SgnSgn": "Signal - Signal",
    "SgnBkg": "Signal - Background",
    "BkgSgn": "Background - Signal",
    "BkgBkg": "Background - Background",
}
COLORS = {
    "SgnSgn": ROOT.kRed + 1,
    "SgnBkg": ROOT.kAzure - 3,
    "BkgSgn": ROOT.kSpring + 3,
    "BkgBkg": ROOT.kOrange + 7,
}
MARKERS = {
    "SgnSgn": ROOT.kFullCircle,
    "SgnBkg": ROOT.kFullDiamond,
    "BkgSgn": ROOT.kFullCross,
    "BkgBkg": ROOT.kFullTriangleUp,
}

TAG_COLORS = {"os": ROOT.kAzure - 3, "ls": ROOT.kRed + 1}
TAG_LABELS = {"os": "OS (D^{0}#bar{D}^{0})", "ls": "LS (D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0})"}

# Axis titles; variables not listed fall back to their name
AXIS_TITLES = {
    "pt1": "#it{p}_{T}^{1} (GeV/#it{c})",
    "pt2": "#it{p}_{T}^{2} (GeV/#it{c})",
    "pairPt": "#it{p}_{T}^{pair} (GeV/#it{c})",
    "pairMass": "#it{m}_{pair} (GeV/#it{c}^{2})",
    "deltaY": "#Delta#it{y}",
    "deltaPhi": "#Delta#varphi (rad)",
}


# ---------------------------------------------------------------- helpers
def discover_variables(f):
    """Variables found as h_<tag>_<var>_<species> in the file."""
    variables = []
    for key in f.GetListOfKeys():
        name = key.GetName()
        for tag in ("os", "ls"):
            for sp in SPECIES:
                pre, suf = f"h_{tag}_", f"_{sp}"
                if name.startswith(pre) and name.endswith(suf):
                    var = name[len(pre):-len(suf)]
                    if var not in variables:
                        variables.append(var)
    return variables


def get_hist(f, tag, var, sp):
    h = f.Get(f"h_{tag}_{var}_{sp}")
    if not h:
        return None
    h = h.Clone(f"{h.GetName()}_plot")
    h.SetDirectory(0)
    return h


def prepare(h, color, norm, marker):
    """Unit-area normalisation (if requested) and styling. Returns the raw integral."""
    raw = h.Integral()
    if norm and raw > 0:
        h.Scale(1.0 / raw)
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(marker)
    h.SetMarkerSize(0.8)
    h.SetLineWidth(2)
    return raw


def set_y_range(hists, logy):
    ymax = max(h.GetMaximum() for h in hists)
    ymin = min(h.GetMinimum() for h in hists)
    if logy:
        lo, hi = 1e-4 * ymax, 5.0 * ymax
    else:
        lo, hi = min(0.0, 1.2 * ymin), 1.4 * ymax   # sWeights can give negative bins
    for h in hists:
        h.SetMinimum(lo)
        h.SetMaximum(hi)


def draw_label(text, x=0.18, y=0.88):
    lt = ROOT.TLatex()
    lt.SetNDC()
    lt.SetTextSize(0.045)
    lt.DrawLatex(x, y, text)


# ---------------------------------------------------------------- plots
def plot_components(f, tag, var, outdir, norm, logy):
    """Four components of one sample (OS or LS) overlaid."""
    hists, entries = [], []
    for sp in SPECIES:
        h = get_hist(f, tag, var, sp)
        if h is None:
            continue
        raw = prepare(h, COLORS[sp], norm, MARKERS[sp])
        hists.append(h)
        entries.append((h, f"{LABELS[sp]}"))
    if not hists:
        print(f"  [{tag}] {var}: no histograms found, skipped")
        return

    c = ROOT.TCanvas(f"c_{tag}_{var}", "", 800, 700)
    c.SetLogy(logy)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    set_y_range(hists, logy)

    first = hists[0]
    first.SetTitle("")
    first.GetXaxis().SetTitle(AXIS_TITLES.get(var, var))
    first.GetYaxis().SetTitle("Normalised sWeighted entries" if norm else "sWeighted entries")
    first.GetYaxis().SetTitleOffset(1.5)
    first.Draw("E1")
    for h in hists[1:]:
        h.Draw("E1 SAME")

    leg = ROOT.TLegend(0.6, 0.73, 0.89, 0.93)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    for h, label in entries:
        leg.AddEntry(h, label, "lep")
    leg.Draw()
    draw_label(TAG_LABELS[tag])

    c.SaveAs(os.path.join(outdir, f"{var}_components_{tag}.pdf"))


def plot_ss_os_vs_ls(f, var, outdir, norm, logy):
    """SgnSgn OS vs LS, with the OS/LS ratio below."""
    h_os = get_hist(f, "os", var, "SgnSgn")
    h_ls = get_hist(f, "ls", var, "SgnSgn")
    if h_os is None or h_ls is None:
        print(f"  {var}: SgnSgn OS or LS missing, skipped")
        return

    raw_os = prepare(h_os, TAG_COLORS["os"], norm, marker=20)
    raw_ls = prepare(h_ls, TAG_COLORS["ls"], norm, marker=24)

    c = ROOT.TCanvas(f"c_ss_{var}", "", 800, 900)
    p1 = ROOT.TPad("p1", "", 0, 0.30, 1, 1)
    p2 = ROOT.TPad("p2", "", 0, 0.00, 1, 0.30)
    p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.14)
    p2.SetTopMargin(0.03);    p2.SetBottomMargin(0.35); p2.SetLeftMargin(0.14)
    p1.Draw(); p2.Draw()

    # --- upper pad
    p1.cd()
    p1.SetLogy(logy)
    set_y_range([h_os, h_ls], logy)
    h_os.SetTitle("")
    h_os.GetYaxis().SetTitle("Normalised sWeighted entries" if norm else "sWeighted entries")
    h_os.GetYaxis().SetTitleOffset(1.5)
    h_os.GetXaxis().SetLabelSize(0)
    h_os.Draw("E1")
    h_ls.Draw("E1 SAME")

    leg = ROOT.TLegend(0.55, 0.72, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_os, f"OS sig-sig (N = {raw_os:.3g})", "lep")
    leg.AddEntry(h_ls, f"LS sig-sig (N = {raw_ls:.3g})", "lep")
    leg.Draw()

    # --- ratio pad
    p2.cd()
    ratio = h_os.Clone(f"ratio_{var}")
    ratio.SetDirectory(0)
    ratio.Divide(h_ls)
    ratio.SetLineColor(ROOT.kAzure-3)
    ratio.SetMarkerColor(ROOT.kAzure-3)
    ratio.SetMinimum(0.0)
    ratio.SetMaximum(3.0)
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("OS / LS")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(0.10); ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetYaxis().SetTitleOffset(0.55)
    ratio.GetXaxis().SetTitle(AXIS_TITLES.get(var, var))
    ratio.GetXaxis().SetTitleSize(0.12); ratio.GetXaxis().SetLabelSize(0.10)
    ratio.Draw("E1")
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0, ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw()

    c.SaveAs(os.path.join(outdir, f"{var}_ss_os_vs_ls.pdf"))

def get_component_hists(f, tag, var, norm):
    """Styled (species, hist) pairs for one sample and variable."""
    out = []
    for sp in SPECIES:
        h = get_hist(f, tag, var, sp)
        if h is None:
            continue
        prepare(h, COLORS[sp], norm, MARKERS[sp])
        out.append((sp, h))
    return out


def get_ss_hists(f, var, norm):
    """Styled (OS, LS) SgnSgn hists, or None if one is missing."""
    h_os = get_hist(f, "os", var, "SgnSgn")
    h_ls = get_hist(f, "ls", var, "SgnSgn")
    if h_os is None or h_ls is None:
        return None
    prepare(h_os, TAG_COLORS["os"], norm, marker=20)
    prepare(h_ls, TAG_COLORS["ls"], norm, marker=24)
    return h_os, h_ls


def draw_overview(outpath, panels, legend_entries, norm, logy, tag_label=None):
    """
    One pad per variable.
    panels         : list of (variable name, [styled hists])
    legend_entries : list of (hist, label) taken from one panel (same style in all pads)
    """
    n = len(panels)
    if n == 0:
        return
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    npads = nrows * ncols

    c = ROOT.TCanvas("c_" + os.path.basename(outpath), "", 600 * ncols, 500 * nrows)
    c.Divide(ncols, nrows, 0.001, 0.001)

    for i, (var, hists) in enumerate(panels, start=1):
        pad = c.cd(i)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.04)
        pad.SetTopMargin(0.06)
        pad.SetBottomMargin(0.14)
        pad.SetLogy(logy)

        set_y_range(hists, logy)
        first = hists[0]
        first.SetTitle("")
        first.GetXaxis().SetTitle(AXIS_TITLES.get(var, var))
        first.GetYaxis().SetTitle("Normalised sWeighted entries" if norm else "sWeighted entries")
        first.GetYaxis().SetTitleOffset(1.6)
        first.Draw("E1")
        for h in hists[1:]:
            h.Draw("E1 SAME")
        if tag_label:
            draw_label(tag_label, x=0.20, y=0.86)

    # Legend only once: in the first spare pad if there is one, otherwise in pad 1
    if npads > n:
        c.cd(n + 1)
        leg = ROOT.TLegend(0.05, 0.25, 0.95, 0.75)
        leg.SetTextSize(0.05)
    else:
        c.cd(1)
        leg = ROOT.TLegend(0.52, 0.66, 0.93, 0.92)
        leg.SetTextSize(0.035)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for h, label in legend_entries:
        leg.AddEntry(h, label, "lep")
    leg.Draw()

    c.SaveAs(outpath)
    print(f"  overview saved: {outpath}")

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="splot_raw_*.root file(s) from SPlotter")
    ap.add_argument("--outdir", default="figures/splot")
    ap.add_argument("--variables", nargs="*", default=None, help="subset; default: all in the file")
    ap.add_argument("--no-norm", action="store_true", help="plot sWeighted counts instead of unit area")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--overview-only", action="store_true",
                    help="skip the per-variable files, only write the overview canvases")
    args = ap.parse_args()

    for path in args.inputs:
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            print(f"Cannot open {path}")
            continue

        outdir = os.path.join(args.outdir, os.path.splitext(os.path.basename(path))[0])
        os.makedirs(outdir, exist_ok=True)
        variables = args.variables or discover_variables(f)
        norm = not args.no_norm
        print(f"{path}: {len(variables)} variable(s) -> {outdir}")

        comp_panels = {"os": [], "ls": []}
        comp_legend = {"os": None, "ls": None}
        ss_panels, ss_legend = [], None

        for var in variables:
            print(f"  {var}")
            if not args.overview_only:
                plot_components(f, "os", var, outdir, norm, args.logy)
                plot_components(f, "ls", var, outdir, norm, args.logy)
                plot_ss_os_vs_ls(f, var, outdir, norm, args.logy)

            # collect for the overview canvases
            for tag in ("os", "ls"):
                hs = get_component_hists(f, tag, var, norm)
                if hs:
                    comp_panels[tag].append((var, [h for _, h in hs]))
                    if comp_legend[tag] is None:
                        comp_legend[tag] = [(h, LABELS[sp]) for sp, h in hs]
            ss = get_ss_hists(f, var, norm)
            if ss:
                ss_panels.append((var, list(ss)))
                if ss_legend is None:
                    ss_legend = [(ss[0], "OS sig-sig"), (ss[1], "LS sig-sig")]

        for tag in ("os", "ls"):
            draw_overview(os.path.join(outdir, f"overview_components_{tag}.pdf"),
                          comp_panels[tag], comp_legend[tag], norm, args.logy,
                          tag_label=TAG_LABELS[tag])
        draw_overview(os.path.join(outdir, "overview_ss_os_vs_ls.pdf"),
                      ss_panels, ss_legend, norm, args.logy)

        f.Close()


if __name__ == "__main__":
    main()