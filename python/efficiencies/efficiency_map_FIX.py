#!/usr/bin/env python3
import ROOT
import json
import os

import sys
sys.path.append("../utils")
import plot_style_utils as plutils

# Use modern ROOT style
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.ForceStyle()

# -----------------------------------------------------------------------------
# Global flags
verbose = False
debug = False
# -----------------------------------------------------------------------------


def efficiency_map():
    global verbose, debug

    # Load JSON config
    json_path = "/home/andrea/MyMacros/Double-Dzero-correlations/.config/config_efficiency_map.json"
    if not os.path.exists(json_path):
        print(f"❌ Error: JSON file not found at {json_path}")
        return

    with open(json_path, "r") as f:
        jsonData = json.load(f)

    verbose = jsonData["verbose"]
    debug = jsonData["debug"]

    # Extract parameters
    fileNameGen = jsonData["files"]["filenameGen"]
    fileNameReco = jsonData["files"]["filenameReco"]
    fileNameData = jsonData["files"]["filenameData"]
    outputDirName = jsonData["outputDir"]
    outputFileName = jsonData["outputFile"]

    hNameReco = jsonData["histos"]["hNameReco"]
    hNameGen = jsonData["histos"]["hNameGen"]
    hNameData = jsonData["histos"]["hNameData"]
    hNameRecoPrompt = jsonData["histos"]["hNameRecoPrompt"]
    hNameGenPrompt = jsonData["histos"]["hNameGenPrompt"]
    hNameRecoNonPrompt = jsonData["histos"]["hNameRecoNonPrompt"]
    hNameGenNonPrompt = jsonData["histos"]["hNameGenNonPrompt"]

    # Load data
    fGen = ROOT.TFile.Open(fileNameGen)
    fReco = ROOT.TFile.Open(fileNameReco)
    fData = ROOT.TFile.Open(fileNameData)

    if not fGen or fGen.IsZombie() or not fReco or fReco.IsZombie() or not fData or fData.IsZombie():
        print("❌ ERROR: One or more input files could not be opened")
        return

    # Prepare output
    os.makedirs(outputDirName, exist_ok=True)
    output_path = os.path.join(outputDirName, outputFileName)
    outputFile = ROOT.TFile(output_path, "RECREATE")

    dirname = "hf-correlator-d-meson-pairs"
    dirReco = fReco.Get(dirname)
    dirGen = fGen.Get(dirname)
    dirData = fData.Get(dirname)

    # Load histograms
    hPtVsYReco = dirReco.Get(hNameReco)
    hPtVsYGen = dirGen.Get(hNameGen)
    hPtVsYRecoPrompt = dirReco.Get(hNameRecoPrompt)
    hPtVsYGenPrompt = dirGen.Get(hNameGenPrompt)
    hPtVsYRecoNonPrompt = dirReco.Get(hNameRecoNonPrompt)
    hPtVsYGenNonPrompt = dirGen.Get(hNameGenNonPrompt)
    hNContrib = dirData.Get(hNameData)

    # Integrals
    dataIntegral = hNContrib.Integral()
    recoIntegral = hPtVsYReco.ProjectionZ().Integral()
    genIntegral = hPtVsYGen.ProjectionZ().Integral()

    hNContribReco = hPtVsYReco.ProjectionZ()
    hNContribGen = hPtVsYGen.ProjectionZ()

    if verbose:
        print("✅ Data loaded successfully")

    # Create 2D maps
    def make_2d_clone(hist3d, name, title):
        h2 = hist3d.Project3D("yx")
        h2.Reset()
        h2.SetNameTitle(name, title)
        return h2

    hEfficiencyMap = make_2d_clone(hPtVsYReco, "hEfficiencyMap", "Efficiency map")
    hPtVsYRecoWeighted = make_2d_clone(hPtVsYReco, "hPtVsYRecoWeighted", "Reco weighted")
    hPtVsYGenWeighted = make_2d_clone(hPtVsYGen, "hPtVsYGenWeighted", "Gen weighted")

    hEfficiencyMapPrompt = make_2d_clone(hPtVsYRecoPrompt, "hEfficiencyMapPrompt", "Efficiency map prompt")
    hPtVsYRecoWeightedPrompt = make_2d_clone(hPtVsYRecoPrompt, "hPtVsYRecoWeightedPrompt", "Reco weighted prompt")
    hPtVsYGenWeightedPrompt = make_2d_clone(hPtVsYGenPrompt, "hPtVsYGenWeightedPrompt", "Gen weighted prompt")

    hEfficiencyMapNonPrompt = make_2d_clone(hPtVsYRecoNonPrompt, "hEfficiencyMapNonPrompt", "Efficiency map nonprompt")
    hPtVsYRecoWeightedNonPrompt = make_2d_clone(hPtVsYRecoNonPrompt, "hPtVsYRecoWeightedNonPrompt", "Reco weighted nonprompt")
    hPtVsYGenWeightedNonPrompt = make_2d_clone(hPtVsYGenNonPrompt, "hPtVsYGenWeightedNonPrompt", "Gen weighted nonprompt")

    # Get bin counts
    nbinsPt = hPtVsYGen.GetXaxis().GetNbins()
    nbinsY = hPtVsYGen.GetYaxis().GetNbins()
    nbinsNContrib = hPtVsYGen.GetZaxis().GetNbins()

    y_min = hPtVsYGen.GetYaxis().FindBin(-0.5)
    y_max = hPtVsYGen.GetYaxis().FindBin(0.5)

    if debug:
        print(f"Nbins Pt: {nbinsPt}, Y: {nbinsY}, NContrib: {nbinsNContrib}")

    # Core computation
    for ipt in range(1, nbinsPt + 1):
        for iy in range(y_min, y_max + 1):
            weightReco = weightGen = 0.0
            weightRecoPrompt = weightGenPrompt = 0.0
            weightRecoNonPrompt = weightGenNonPrompt = 0.0

            for inb in range(1, nbinsNContrib + 1):
                binReco = hPtVsYReco.GetBinContent(ipt, iy, inb)
                binGen = hPtVsYGen.GetBinContent(ipt, iy, inb)
                binRecoPrompt = hPtVsYRecoPrompt.GetBinContent(ipt, iy, inb)
                binGenPrompt = hPtVsYGenPrompt.GetBinContent(ipt, iy, inb)
                binRecoNonPrompt = hPtVsYRecoNonPrompt.GetBinContent(ipt, iy, inb)
                binGenNonPrompt = hPtVsYGenNonPrompt.GetBinContent(ipt, iy, inb)

                nData = hNContrib.GetBinContent(inb)
                nReco = hNContribReco.GetBinContent(inb)
                nGen = hNContribGen.GetBinContent(inb)

                dataRecoWeight = dataGenWeight = 1.0
                if all(v > 0 for v in [nReco, nGen, nData, recoIntegral, genIntegral]):
                    dataRecoWeight = (nData / dataIntegral) / (nReco / recoIntegral)
                    dataGenWeight = (nData / dataIntegral) / (nGen / genIntegral)

                weightReco += dataRecoWeight * binReco
                weightGen += dataGenWeight * binGen

                weightRecoPrompt += dataRecoWeight * binRecoPrompt
                weightGenPrompt += dataGenWeight * binGenPrompt

                weightRecoNonPrompt += dataRecoWeight * binRecoNonPrompt
                weightGenNonPrompt += dataGenWeight * binGenNonPrompt

            # Fill weighted maps
            if weightReco > 0:
                hPtVsYRecoWeighted.SetBinContent(ipt, iy, weightReco)
            if weightGen > 0:
                hPtVsYGenWeighted.SetBinContent(ipt, iy, weightGen)

            if weightRecoPrompt > 0:
                hPtVsYRecoWeightedPrompt.SetBinContent(ipt, iy, weightRecoPrompt)
            if weightGenPrompt > 0:
                hPtVsYGenWeightedPrompt.SetBinContent(ipt, iy, weightGenPrompt)

            if weightRecoNonPrompt > 0:
                hPtVsYRecoWeightedNonPrompt.SetBinContent(ipt, iy, weightRecoNonPrompt)
            if weightGenNonPrompt > 0:
                hPtVsYGenWeightedNonPrompt.SetBinContent(ipt, iy, weightGenNonPrompt)

    # Rebin (optional, adjust factor)
    for h in [
        hPtVsYRecoWeighted,
        hPtVsYGenWeighted,
        hEfficiencyMap,
        hPtVsYRecoWeightedPrompt,
        hPtVsYGenWeightedPrompt,
        hEfficiencyMapPrompt,
        hPtVsYRecoWeightedNonPrompt,
        hPtVsYGenWeightedNonPrompt,
        hEfficiencyMapNonPrompt,
    ]:
        h.RebinX(5)

    # Compute efficiency maps
    for ipt in range(1, hPtVsYRecoWeighted.GetXaxis().GetNbins() + 1):
        for iy in range(y_min, y_max):
            def safe_ratio(num, den):
                return num / den if den != 0 else 0

            eff = safe_ratio(hPtVsYRecoWeighted.GetBinContent(ipt, iy), hPtVsYGenWeighted.GetBinContent(ipt, iy))
            effPrompt = safe_ratio(hPtVsYRecoWeightedPrompt.GetBinContent(ipt, iy), hPtVsYGenWeightedPrompt.GetBinContent(ipt, iy))
            effNonPrompt = safe_ratio(hPtVsYRecoWeightedNonPrompt.GetBinContent(ipt, iy), hPtVsYGenWeightedNonPrompt.GetBinContent(ipt, iy))

            hEfficiencyMap.SetBinContent(ipt, iy, eff)
            hEfficiencyMapPrompt.SetBinContent(ipt, iy, effPrompt)
            hEfficiencyMapNonPrompt.SetBinContent(ipt, iy, effNonPrompt)

    # Style and save
    hEfficiencyMap.GetZaxis().SetTitle("(Acc x #epsilon)")
    hEfficiencyMap.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")

    c2d = ROOT.TCanvas("cEfficiencyMap", "Efficiency map", 1000, 700)
    hEfficiencyMap.Draw("colz")

    plutils.draw_all_info(0.12, 0.4, 0.8, 0.94, 'This Thesis', 'pp', '13.6', 'D0')

    """  ptTitle = ROOT.TPaveText(0.12, 0.8, 0.4, 0.94, "brNDC")
    ptTitle.SetTextFont(42)
    ptTitle.AddText("#scale[1.35]{This Thesis}")
    ptTitle.Draw("same")

    pt = ROOT.TPaveText(0.12, 0.7, 0.5, 0.9, "brNDC")
    pt.SetTextFont(42)
    pt.AddText("#scale[1.2]{pp, #sqrt{#it{s}} = 13.6 TeV}")
    pt.AddText("D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.")
    pt.Draw("same") """

    # Save histograms and canvases
    outputFile.cd()
    hEfficiencyMap.Write()
    hEfficiencyMapPrompt.Write()
    hEfficiencyMapNonPrompt.Write()
    c2d.Write()
    outputFile.Close()

    print(f"✅ Plots and histograms saved to {output_path}")
    ROOT.gApplication.Terminate(0)


# Run if executed directly
if __name__ == "__main__":
    efficiency_map()
