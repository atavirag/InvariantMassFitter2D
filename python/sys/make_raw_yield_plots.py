import ROOT
import re
from array import array

def draw_filled_function(tf1, color=ROOT.kOrange-3, alpha=0.4, npoints=300):
    xmin, xmax = tf1.GetXmin(), tf1.GetXmax()
    x_vals = [xmin + i * (xmax - xmin) / npoints for i in range(npoints + 1)]
    y_vals = [tf1.Eval(x) for x in x_vals]

    # Add baseline to form the filled shape
    x_vals += [xmax, xmin]
    y_vals += [0, 0]

    graph = ROOT.TGraph(len(x_vals), array('d', x_vals), array('d', y_vals))
    graph.SetFillStyle(1001)  # 1001 = solid fill
    graph.SetFillColorAlpha(color, alpha)
    graph.SetLineColor(color)
    graph.SetLineWidth(3)
    graph.Draw("F SAME")
    graph.Draw("L SAME")
    return graph  # return it so it doesn't get garbage collected

# Open the ROOT file in read-only mode
f = ROOT.TFile.Open("./rawYield_swapped_1_5.root", "READ")

# Get the histograms with fit results
hSigma = f.Get("hRawYieldsSigma")
hMean = f.Get("hRawYieldsMean")
hSignificance = f.Get("hRawYieldsSignificance")
hSignal = f.Get("hRawYieldsSignal")
hReflOverSgn = f.Get("hReflectionOverSignal")

# Check if histograms exist
if not all([hSigma, hMean, hSignificance, hSignal]):
    print("Error: Could not find required histograms in the ROOT file")
    f.Close()
    exit(1)

# Get all histogram names that match the MassForFit pattern
keys = f.GetListOfKeys()
mass_hists = [k.GetName() for k in keys if re.match(r"MassForFit\d*", k.GetName())]

ROOT.TGaxis.SetMaxDigits(3)

print("Found mass histograms:", mass_hists)

canvas = ROOT.gROOT.FindObject("canvasMass0")

# Print basic canvas information
print(f"Canvas name: {canvas.GetName()}")
print(f"Canvas title: {canvas.GetTitle()}")
print(f"Canvas width: {canvas.GetWw()}")
print(f"Canvas height: {canvas.GetWh()}")

# Get list of all objects drawn on the canvas
canvas.cd()  # Make sure we're in the right canvas
objects = canvas.GetListOfPrimitives()

print(f"Number of objects in canvas: {objects.GetEntries()}")

# Initialize the functions list (start empty, we'll fill it)
functions = []
fTot, fBkg, fRefl = None, None, None
hist = f.Get("MassForFit0")

# Get objects from canvas
objects = canvas.GetListOfPrimitives()  # Assuming 'canvas' is your original canvas

# Print each object and extract the functions you want
for i in range(objects.GetEntries()):
    obj = objects.At(i)
    print(f"Object {i}: {obj.GetName()} (Type: {obj.ClassName()})")
    
    # Check if it's a function (TF1) and extract by name pattern
    if obj.InheritsFrom("RooCurve"):
        print(obj.GetName())
        if "Tot" in obj.GetName():
            fTot = obj.Clone(f"{obj.GetName()}_copy")
            functions.append(fTot)
        elif "Bkg" in obj.GetName():
            fBkg = obj.Clone(f"{obj.GetName()}_copy")
            functions.append(fBkg)
        elif "Refl" in obj.GetName():
            fRefl = obj.Clone(f"{obj.GetName()}_copy")
            functions.append(fRefl)

# Create new canvas
canvas = ROOT.TCanvas(f"c0", "", 800, 800)
canvas.SetTopMargin(0.04)
canvas.SetLeftMargin(0.15)
canvas.SetRightMargin(0.03)
canvas.SetBottomMargin(0.12)
ROOT.gStyle.SetOptStat(0)


# Add ticks on all sides
canvas.SetTickx()
canvas.SetTicky()

# Add a grid
canvas.SetGridx()
canvas.SetGridy()

ROOT.gStyle.SetGridColor(ROOT.kGray + 2)
ROOT.gStyle.SetGridStyle(3)  # Dotted
ROOT.gStyle.SetGridWidth(1)

# Draw histogram
hist.SetMarkerStyle(20)
hist.SetMarkerSize(1.5)
hist.SetMarkerColor(ROOT.kBlack)

hist.SetTitle("")
hist.GetXaxis().SetTitle("m_{K#pi} (GeV/#it{c}^{2})")
hist.GetYaxis().SetTitle("Counts per 5 MeV/#it{c}^{2}")
hist.GetXaxis().SetTitleSize(0.045)
hist.GetYaxis().SetTitleSize(0.045)
hist.GetXaxis().SetLabelSize(0.04)
hist.GetYaxis().SetLabelSize(0.04)

# Style and draw functions
fTot.SetLineColor(ROOT.kBlue - 4)
fBkg.SetLineColor(ROOT.kRed)
fBkg.SetLineStyle(ROOT.kDashed)
if fRefl:
    fRefl.SetLineColor(ROOT.kGreen + 2)
    fRefl.SetLineStyle(9)
    fRefl.SetLineWidth(4)
    
fTot.SetLineWidth(5)
fBkg.SetLineWidth(4)

    
# Determine max y from fTot by sampling points
n_points = 500

xmin = fTot.GetXaxis().GetXmin()
xmax = fTot.GetXaxis().GetXmax()

# Set x-axis range to the fit function range
hist.GetXaxis().SetRangeUser(1.74, 2.05)

ymax = max([fTot.Eval(xmin + i * (xmax - xmin) / n_points) for i in range(n_points + 1)])

# Set Y range a bit above the max
hist.GetYaxis().SetRangeUser(0, ymax * 1.3)

hist.Draw("E")
    
# Draw other functions
fBkg.Draw("SAME")
if fRefl:
    fRefl.Draw("SAME")
fTot.Draw("SAME")

# Legend
legend = ROOT.TLegend(0.16, 0.23, 0.45, 0.46)
legend.SetTextFont(42)
legend.SetBorderSize(0)     # Removes the outer box border
legend.SetFillStyle(0)      # Makes the background transparent
legend.SetTextSize(0.037)
legend.AddEntry(hist, "Data", "lep")
legend.AddEntry(fTot, "Total fit function", "l")
if fRefl:
    legend.AddEntry(fRefl, r"Reflected K#minus#pi", "l")
legend.AddEntry(fBkg, "Background", "l")
legend.Draw()

sigma_val = hSigma.GetBinContent(1)*1000
sigma_err = hSigma.GetBinError(1)*1000

mean_val = hMean.GetBinContent(1)
mean_err = hMean.GetBinError(1)

significance_val = hMean.GetBinContent(1)
significance_err = hSignificance.GetBinError(1)

soverb_val = hSignal.GetBinContent(1)
soverb_err = hSignal.GetBinError(1)

reflOverSgn = hReflOverSgn.GetBinContent(1)

# Add text: pT range and "This Thesis"
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.045)
latex.DrawLatex(0.17, 0.88, "This Thesis")
latex.SetTextSize(0.04)
latex.DrawLatex(0.57, 0.26, f"{5.0} < #it{{p}}_{{T}} < {24.0} GeV/#it{{c}}")
latex.DrawLatex(0.61, 0.2, r"pp, #sqrt{s} = 13.6 TeV")
latex.DrawLatex(0.52, 0.88, r"D^{0}#rightarrow K^{#minus}#pi^{+} and charge conj.")
latex.SetTextSize(0.035)
    
if sigma_val > 0:
    latex.DrawLatex(0.57, 0.83, f"#sigma = ({sigma_val:.2f} #pm {sigma_err:.2f}) MeV/#it{{c}}^{{2}}")
if mean_val > 0:
    latex.DrawLatex(0.57, 0.78, f"#mu = ({mean_val:.2f} #pm {mean_err:.2f}) MeV/#it{{c}}^{{2}}")
if significance_val > 0:
    latex.DrawLatex(0.57, 0.73, f"Signif.(3#sigma) = {significance_val:.1f}")
if soverb_val > 0:
    latex.DrawLatex(0.57, 0.68, f"S = {soverb_val:.0f} #pm {soverb_err:.0f}")
if reflOverSgn > 0:
    latex.DrawLatex(0.57, 0.63, f"Refl./Signal = {reflOverSgn:.3f}")

canvas.Update()
canvas.SaveAs(f"rawYield_swapped_1_5_bonito.pdf")
print(f"Canvas saved as rawYield_swapped_5_24_wihtRefl_bonito.pdf")
canvas.Close()

# Close the ROOT file safely
f.Close()
print("Done processing all histograms")
