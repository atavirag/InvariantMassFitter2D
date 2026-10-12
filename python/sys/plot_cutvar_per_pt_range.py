import ROOT

# Open the ROOT file in read-only mode
f = ROOT.TFile.Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/CutVarD0_2024_triggered.root", "READ")

# List of canvas names
canvases = [
    "cRawYieldVsCut_pt1_2",
    "cRawYieldVsCut_pt2_3",
    "cRawYieldVsCut_pt3_4",
    "cRawYieldVsCut_pt4_5",
    "cRawYieldVsCut_pt5_6",
    "cRawYieldVsCut_pt6_7",
    "cRawYieldVsCut_pt7_8",
    "cRawYieldVsCut_pt8_10",
    "cRawYieldVsCut_pt10_12",
    "cRawYieldVsCut_pt12_16",
    "cRawYieldVsCut_pt16_24",
]

ROOT.TGaxis.SetMaxDigits(3)

for cname in canvases:
    canvas = f.Get(cname)
    if not canvas:
        print(f"Canvas {cname} not found!")
        continue

    new_canvas = ROOT.TCanvas(f"{cname}_new", "", 800, 800)

    new_canvas.cd()
    canvas.DrawClonePad()

    pad_clone = new_canvas.GetPad(0)  # get the pad that was just cloned/drawn
    if pad_clone:
        pad_clone.SetFillStyle(0)
        pad_clone.SetBorderSize(0)

        pad_clone.SetTopMargin(0.05)
        pad_clone.SetLeftMargin(0.17)
        pad_clone.SetRightMargin(0.02)
        pad_clone.SetBottomMargin(0.12)

        pad_clone.Modified()

    # Now add your text overlays
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)

    pt_range = cname.split("_pt")[-1]
    pt_low, pt_high = pt_range.split("_")

    latex.DrawLatex(0.6, 0.6, f"{pt_low} < p_{{T}} < {pt_high} GeV/c")
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.21, 0.88, "This Thesis")

    new_canvas.Modified()
    new_canvas.Update()

    new_canvas.SaveAs(f"/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/cRawYieldVsCut_pt{pt_low}_{pt_high}.pdf")
    new_canvas.Close()

f.Close()
