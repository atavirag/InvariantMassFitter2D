#!/usr/bin/env python3
"""
divide_se_me.py - SE/ME correlation: divides the 2D sWeighted SE histograms by the 2D ME histogram.

SE file : splot_raw_*.root from SPlotter   (h2_<os|ls>_<x>_<y>_SgnSgn)
ME file : splot_me_*.root  from MESPlotter (h2_me_<x>_<y>_SgnSgn_norm, or _SgnSgn with --me-raw)

Output: <outdir>/<se-file-name>_over_me.root with h2_<tag>_<x>_<y>_SEoverME, plus one LEGO pdf per ratio.

Usage:
  python divide_se_me.py --se splot_raw_....root --me splot_me_....root \
        --outdir figures/se_over_me [--tags os ls] [--variables2D deltaEta_deltaPhi] [--me-raw]
"""
import argparse
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

SPECIES = ["SgnSgn", "SgnBkg", "BkgSgn", "BkgBkg"]
AXIS_TITLES = {
    "deltaY": "#Delta#it{y}",
    "deltaPhi": "#Delta#varphi (rad)",
    "deltaEta": "#Delta#eta",
    "deltaPhiCorrelations": "#Delta#varphi (rad)",
}
TAG_LABELS = {"os": "OS (D^{0}#bar{D}^{0})", "ls": "LS (D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0})"}
TAG_COLORS = {"os": ROOT.kAzure - 3, "ls": ROOT.kRed + 1}


def discover(f, tag):
    """2D keys '<x>_<y>' found as h2_<tag>_<x>_<y>_SgnSgn."""
    pre, suf = f"h2_{tag}_", "_SgnSgn"
    return [k.GetName()[len(pre):-len(suf)] for k in f.GetListOfKeys()
            if k.GetName().startswith(pre) and k.GetName().endswith(suf)]


def clone(f, name):
    h = f.Get(name)
    if not h:
        return None
    h = h.Clone(f"{name}_clone")
    h.SetDirectory(0)
    if h.GetSumw2N() == 0:
        h.Sumw2()
    return h


def same_binning(a, b):
    ax, bx, ay, by = a.GetXaxis(), b.GetXaxis(), a.GetYaxis(), b.GetYaxis()
    return (ax.GetNbins() == bx.GetNbins() and ay.GetNbins() == by.GetNbins()
            and abs(ax.GetXmin() - bx.GetXmin()) < 1e-9 and abs(ax.GetXmax() - bx.GetXmax()) < 1e-9
            and abs(ay.GetXmin() - by.GetXmin()) < 1e-9 and abs(ay.GetXmax() - by.GetXmax()) < 1e-9)


def axis_titles(key, h):
    parts = key.split("_")
    if len(parts) == 2:
        return AXIS_TITLES.get(parts[0], parts[0]), AXIS_TITLES.get(parts[1], parts[1])
    return h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()


def plot_ratio(h, tag, key, outdir):
    xt, yt = axis_titles(key, h)
    h.SetTitle("")
    h.GetXaxis().SetTitle(xt)
    h.GetYaxis().SetTitle(yt)
    h.GetZaxis().SetTitle("SE / ME")
    for ax, off in ((h.GetXaxis(), 1.8), (h.GetYaxis(), 1.8), (h.GetZaxis(), 1.4)):
        ax.SetTitleOffset(off)
        ax.SetTitleSize(0.04)
        ax.SetLabelSize(0.035)
    h.SetFillColor(TAG_COLORS[tag])
    h.SetLineColor(ROOT.kBlack)

    c = ROOT.TCanvas(f"c_{tag}_{key}", "", 900, 750)
    c.SetLeftMargin(0.14); c.SetRightMargin(0.06); c.SetTopMargin(0.08); c.SetBottomMargin(0.14)
    c.SetTheta(25); c.SetPhi(-125)
    h.Draw("LEGO2")
    lt = ROOT.TLatex()
    lt.SetNDC(); lt.SetTextSize(0.045)
    lt.DrawLatex(0.14, 0.94, f"{TAG_LABELS[tag]}: sig-sig SE / ME")
    c.SaveAs(os.path.join(outdir, f"{key}_SEoverME_{tag}.pdf"))

def plot_projections(hse, hme, tag, key, outdir, fout):
    """SE/ME projected on each axis (project SE and ME separately, then divide)."""
    xt, yt = axis_titles(key, hse)
    c = ROOT.TCanvas(f"c_proj_{tag}_{key}", "", 1400, 600)
    c.Divide(2, 1)
    keep = []
    for i, (axis, title) in enumerate((("X", xt), ("Y", yt)), start=1):
        pad = c.cd(i)
        pad.SetLeftMargin(0.14); pad.SetRightMargin(0.05); pad.SetTicks(); pad.SetGrid()
        pse = getattr(hse, f"Projection{axis}")(f"p{axis}_se_{tag}_{key}")
        pme = getattr(hme, f"Projection{axis}")(f"p{axis}_me_{tag}_{key}")
        pr = pse.Clone(f"p{axis}_{tag}_{key}_SEoverME")
        pr.Divide(pse, pme)
        pr.SetTitle("")
        pr.GetXaxis().SetTitle(title)
        pr.GetYaxis().SetTitle("SE / ME (projection)")
        pr.GetYaxis().SetTitleOffset(1.5)
        pr.SetLineColor(TAG_COLORS[tag]); pr.SetMarkerColor(TAG_COLORS[tag])
        pr.SetMarkerStyle(20); pr.SetLineWidth(2)
        pr.Draw("E1")
        lt = ROOT.TLatex()
        lt.SetNDC(); lt.SetTextSize(0.045)
        lt.DrawLatex(0.18, 0.88, f"{TAG_LABELS[tag]}: sig-sig SE / ME")
        fout.cd()
        pr.Write()
        keep += [pse, pme, pr, lt]
    c.SaveAs(os.path.join(outdir, f"{key}_SEoverME_proj_{tag}.pdf"))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--se", required=True, help="SE sPlot file (SPlotter output)")
    ap.add_argument("--me", required=True, help="ME sPlot file (MESPlotter output)")
    ap.add_argument("--outdir", default="figures/se_over_me")
    ap.add_argument("--tags", nargs="+", default=["os", "ls"], choices=["os", "ls"])
    ap.add_argument("--variables2D", nargs="*", default=None,
                    help="subset as <x>_<y>; default: all pairs present in both files")
    ap.add_argument("--me-raw", action="store_true",
                    help="divide by h2_me_*_SgnSgn instead of the plateau-normalized _SgnSgn_norm")
    ap.add_argument("--me-min-frac", type=float, default=0.1,
                    help="mask bins where ME content < frac * max(ME) (default 0.1)")
    ap.add_argument("--max-relerr", type=float, default=0.3,
                    help="mask bins where the ME relative error exceeds this (default 0.3; <=0 disables)")
    ap.add_argument("--xrange", type=float, nargs=2, default=None, metavar=("LO", "HI"),
                    help="restrict the displayed first-axis range (e.g. --xrange -1.2 1.2)")
    args = ap.parse_args()

    fse = ROOT.TFile.Open(args.se, "READ")
    fme = ROOT.TFile.Open(args.me, "READ")
    if not fse or fse.IsZombie() or not fme or fme.IsZombie():
        print("Cannot open input files")
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    out_name = os.path.splitext(os.path.basename(args.se))[0] + "_over_me.root"
    fout = ROOT.TFile.Open(os.path.join(args.outdir, out_name), "RECREATE")
    me_suffix = "SgnSgn" if args.me_raw else "SgnSgn_norm"

    for tag in args.tags:
        keys = args.variables2D or discover(fse, tag)
        for key in keys:
            hse = clone(fse, f"h2_{tag}_{key}_SgnSgn")
            hme = clone(fme, f"h2_me_{key}_{me_suffix}")
            if hse is None or hme is None:
                print(f"  [{tag}] {key}: missing SE or ME histogram, skipped")
                continue
            if not same_binning(hse, hme):
                print(f"  [{tag}] {key}: SE and ME binning differ, skipped")
                continue

            ratio = hse.Clone(f"h2_{tag}_{key}_SEoverME")
            ratio.SetDirectory(0)
            ratio.Divide(hse, hme)
            nm = mask_low_me([ratio, hse, hme], hme, args.me_min_frac, args.max_relerr)
            print(f"  [{tag}] {key}: masked {nm}/{ratio.GetNbinsX() * ratio.GetNbinsY()} bins")
            if args.xrange:
                for h in (ratio, hse, hme):
                    h.GetXaxis().SetRangeUser(*args.xrange)
            fout.cd()
            ratio.Write()
            plot_ratio(ratio, tag, key, args.outdir)
            plot_projections(hse, hme, tag, key, args.outdir, fout)
            print(f"  [{tag}] {key}: done")

    fout.Close()
    print(f"Written {os.path.join(args.outdir, out_name)}")
    return 0

def mask_low_me(hists, hme, min_frac, max_relerr):
    """Zero every histogram in `hists` where the ME is too small or too uncertain.
    Returns the number of masked bins."""
    ref = hme.GetMaximum()
    nmask = 0
    for ix in range(1, hme.GetNbinsX() + 1):
        for iy in range(1, hme.GetNbinsY() + 1):
            c, e = hme.GetBinContent(ix, iy), hme.GetBinError(ix, iy)
            bad = c <= min_frac * ref
            if max_relerr > 0 and c > 0 and e / c > max_relerr:
                bad = True
            if bad:
                for h in hists:
                    h.SetBinContent(ix, iy, 0.)
                    h.SetBinError(ix, iy, 0.)
                nmask += 1
    return nmask

if __name__ == "__main__":
    raise SystemExit(main())