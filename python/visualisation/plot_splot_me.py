#!/usr/bin/env python3
"""
plot_splot_me.py - visualise the ME sPlot distributions written by MESPlotter.

For every variable it produces:
  1. <var>_components_me.pdf : the four components (SgnSgn, SgnBkg, BkgSgn, BkgBkg) overlaid
For every 2D pair <x>_<y>:
  2. <x>_<y>_components_me_2D.pdf : four LEGO pads
  3. <x>_<y>_ss_me_2D.pdf         : SgnSgn LEGO
  4. <x>_<y>_ss_norm_me_2D.pdf    : SgnSgn normalized to the plateau (if h2_*_SgnSgn_norm exists)
plus overview canvases with all variables.

Usage:
  python plot_splot_me.py splot_me_mass_1.73_2.05_eta_0.7.root \
        --outdir figures/splot_me [--variables pt1 deltaPhi] [--no-norm] [--logy]
"""
import argparse
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

TAG = "me"
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
SAMPLE_LABEL = "Mixed events"

AXIS_TITLES = {
    "pt1": "#it{p}_{T}^{1} (GeV/#it{c})",
    "pt2": "#it{p}_{T}^{2} (GeV/#it{c})",
    "pairPt": "#it{p}_{T}^{pair} (GeV/#it{c})",
    "pairMass": "#it{m}_{pair} (GeV/#it{c}^{2})",
    "deltaY": "#Delta#it{y}",
    "deltaPhi": "#Delta#varphi (rad)",
    "deltaEta": "#Delta#eta",
    "deltaPhiCorrelations": "#Delta#varphi (rad)",
    "xJ": "#it{x}_{J}",
}


# ---------------------------------------------------------------- helpers
def discover(f, prefix):
    """Names X found as <prefix>X_<species> in the file."""
    found = []
    for key in f.GetListOfKeys():
        name = key.GetName()
        for sp in SPECIES:
            suf = f"_{sp}"
            if name.startswith(prefix) and name.endswith(suf):
                x = name[len(prefix):-len(suf)]
                if x and x not in found:
                    found.append(x)
    return found


def get_hist(f, name):
    h = f.Get(name)
    if not h:
        return None
    h = h.Clone(f"{name}_plot")
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


def get_component_hists(f, var, norm):
    """Styled (species, hist) pairs for one variable."""
    out = []
    for sp in SPECIES:
        h = get_hist(f, f"h_{TAG}_{var}_{sp}")
        if h is None:
            continue
        prepare(h, COLORS[sp], norm, MARKERS[sp])
        out.append((sp, h))
    return out


# ---------------------------------------------------------------- 1D plots
def plot_components(f, var, outdir, norm, logy):
    pairs = get_component_hists(f, var, norm)
    if not pairs:
        print(f"  {var}: no histograms found, skipped")
        return
    hists = [h for _, h in pairs]

    c = ROOT.TCanvas(f"c_{TAG}_{var}", "", 800, 700)
    c.SetLogy(logy)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetTicks()
    c.SetGrid()
    set_y_range(hists, logy)

    first = hists[0]
    first.SetTitle("")
    first.GetXaxis().SetTitle(AXIS_TITLES.get(var, var))
    first.GetYaxis().SetTitle("Normalised sWeighted entries" if norm else "sWeighted entries")
    first.GetYaxis().SetTitleOffset(1.5)
    first.Draw("E1")
    for h in hists[1:]:
        h.Draw("E1 SAME")

    leg = ROOT.TLegend(0.6, 0.73, 0.86, 0.93)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    for sp, h in pairs:
        leg.AddEntry(h, LABELS[sp], "lep")
    leg.Draw()
    draw_label(SAMPLE_LABEL)

    c.SaveAs(os.path.join(outdir, f"{var}_components_{TAG}.pdf"))


def draw_overview(outpath, panels, legend_entries, norm, logy):
    """One pad per variable. panels: list of (variable, [styled hists])."""
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
        draw_label(SAMPLE_LABEL, x=0.20, y=0.86)

    # Legend once: in the first spare pad if there is one, otherwise in pad 1
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


# ---------------------------------------------------------------- 2D helpers
def axis_titles_2d(key, h):
    """(x title, y title). Assumes variable names contain no underscore; otherwise
    falls back to the titles stored in the histogram."""
    parts = key.split("_")
    if len(parts) == 2:
        return AXIS_TITLES.get(parts[0], parts[0]), AXIS_TITLES.get(parts[1], parts[1])
    return h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()


def prepare2d(h, color, norm):
    """Unit-volume normalisation (if requested) and LEGO styling. Returns the raw integral."""
    raw = h.Integral()
    if norm and raw > 0:
        h.Scale(1.0 / raw)
    h.SetFillColor(color)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(1)
    return raw


def style_lego_pad(pad, logz):
    pad.SetLeftMargin(0.14)
    pad.SetRightMargin(0.06)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.14)
    pad.SetTheta(25)
    pad.SetPhi(-125)
    pad.SetLogz(logz)


def draw_lego(h, key, ztitle):
    xt, yt = axis_titles_2d(key, h)
    h.SetTitle("")
    h.GetXaxis().SetTitle(xt)
    h.GetYaxis().SetTitle(yt)
    h.GetZaxis().SetTitle(ztitle)
    for ax, off in ((h.GetXaxis(), 1.8), (h.GetYaxis(), 1.8), (h.GetZaxis(), 1.4)):
        ax.SetTitleOffset(off)
        ax.SetTitleSize(0.04)
        ax.SetLabelSize(0.035)
    h.Draw("LEGO2")


# ---------------------------------------------------------------- 2D plots
def plot_components_2d(f, key, outdir, norm, logz):
    pads = []
    for sp in SPECIES:
        h = get_hist(f, f"h2_{TAG}_{key}_{sp}")
        if h is None:
            continue
        prepare2d(h, COLORS[sp], norm)
        pads.append((sp, h))
    if not pads:
        print(f"  {key}: no 2D histograms found, skipped")
        return

    ztitle = "Normalised sWeighted entries" if norm else "sWeighted entries"
    c = ROOT.TCanvas(f"c2_{TAG}_{key}", "", 1200, 1000)
    c.Divide(2, 2, 0.001, 0.001)
    for i, (sp, h) in enumerate(pads, start=1):
        pad = c.cd(i)
        style_lego_pad(pad, logz)
        draw_lego(h, key, ztitle)
        draw_label(f"{SAMPLE_LABEL}: {LABELS[sp]}", x=0.14, y=0.94)

    c.SaveAs(os.path.join(outdir, f"{key}_components_{TAG}_2D.pdf"))


def plot_ss_2d(f, key, outdir, norm, logz):
    """SgnSgn only, single large LEGO canvas."""
    h = get_hist(f, f"h2_{TAG}_{key}_SgnSgn")
    if h is None:
        print(f"  {key}: SgnSgn missing (2D), skipped")
        return
    raw = prepare2d(h, COLORS["SgnSgn"], norm)

    c = ROOT.TCanvas(f"c2_ss_{TAG}_{key}", "", 900, 750)
    style_lego_pad(c, logz)
    draw_lego(h, key, "Normalised sWeighted entries" if norm else "sWeighted entries")
    draw_label(f"{SAMPLE_LABEL}: sig-sig (N = {raw:.3g})", x=0.14, y=0.94)
    c.SaveAs(os.path.join(outdir, f"{key}_ss_{TAG}_2D.pdf"))


def plot_ss_norm_2d(f, key, outdir, logz):
    """SgnSgn already scaled to 1 on the plateau by MESPlotter (no extra normalisation)."""
    h = get_hist(f, f"h2_{TAG}_{key}_SgnSgn_norm")
    if h is None:
        return
    prepare2d(h, COLORS["SgnSgn"], norm=False)

    c = ROOT.TCanvas(f"c2_ssn_{TAG}_{key}", "", 900, 750)
    style_lego_pad(c, logz)
    draw_lego(h, key, "Pairs (plateau = 1)")
    draw_label(f"{SAMPLE_LABEL}: sig-sig, normalized", x=0.14, y=0.94)
    c.SaveAs(os.path.join(outdir, f"{key}_ss_norm_{TAG}_2D.pdf"))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="splot_me_*.root file(s) from MESPlotter")
    ap.add_argument("--outdir", default="figures/splot_me")
    ap.add_argument("--variables", nargs="*", default=None, help="subset; default: all in the file")
    ap.add_argument("--variables2D", nargs="*", default=None,
                    help="2D subset as <x>_<y> (e.g. deltaEta_deltaPhi); default: all in the file")
    ap.add_argument("--no-norm", action="store_true", help="plot sWeighted counts instead of unit area")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--no-2d", action="store_true", help="skip the 2D LEGO plots")
    ap.add_argument("--overview-only", action="store_true",
                    help="skip the per-variable files, only write the overview canvas")
    args = ap.parse_args()

    for path in args.inputs:
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            print(f"Cannot open {path}")
            continue

        outdir = os.path.join(args.outdir, os.path.splitext(os.path.basename(path))[0])
        os.makedirs(outdir, exist_ok=True)
        variables = args.variables or discover(f, f"h_{TAG}_")
        norm = not args.no_norm
        print(f"{path}: {len(variables)} variable(s) -> {outdir}")

        panels, legend = [], None
        for var in variables:
            print(f"  {var}")
            if not args.overview_only:
                plot_components(f, var, outdir, norm, args.logy)
            pairs = get_component_hists(f, var, norm)
            if pairs:
                panels.append((var, [h for _, h in pairs]))
                if legend is None:
                    legend = [(h, LABELS[sp]) for sp, h in pairs]

        if not args.no_2d:
            keys2d = args.variables2D or discover(f, f"h2_{TAG}_")
            print(f"  {len(keys2d)} 2D variable pair(s)")
            for key in keys2d:
                print(f"  {key} (2D)")
                plot_components_2d(f, key, outdir, norm, args.logy)
                plot_ss_2d(f, key, outdir, norm, args.logy)
                plot_ss_norm_2d(f, key, outdir, args.logy)

        draw_overview(os.path.join(outdir, f"overview_components_{TAG}.pdf"),
                      panels, legend, norm, args.logy)
        f.Close()


if __name__ == "__main__":
    main()