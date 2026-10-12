import ROOT

def fit_histogram():
    # Open the ROOT file
    file = ROOT.TFile.Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/CutVarD0_swapped_models.root", "READ")
    if not file or file.IsZombie():
        print("Error: Could not open file.")
        return
    
    # Get the histogram
    hist = file.Get("hRawFracPrompt")
    if not hist:
        print("Error: Histogram 'hRawFracPrompt' not found in file.")
        file.Close()
        return
    
    # Define the fitting function (pol0: constant function)
    # Customize histogram appearance
    hist.SetTitle("")  # Set title
    hist.SetMarkerSize(1.2)
    hist.SetLineWidth(3)
    hist.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")  # X-axis label
    hist.GetYaxis().SetTitle("Prompt fraction")  # Y-axis label
    hist.GetXaxis().SetTitleSize(0.045)  # X-axis label
    hist.GetYaxis().SetTitleSize(0.045)  # Y-axis label
    hist.GetXaxis().SetLabelSize(0.04)  # X-axis label
    hist.GetYaxis().SetLabelSize(0.04)  # Y-axis label
    hist.GetXaxis().SetRangeUser(1, 24)  # Set X-axis range
    hist.GetYaxis().SetRangeUser(0.95, 1.01)   # Set Y-axis range
    
    # Define the fitting function (pol0: constant function)
    fit_func = ROOT.TF1("fit_func", "pol0", 1.0, 24.0)  # Fit within specified range
    fit_func.SetLineColor(ROOT.kBlue-4)  # Set fit line color to red
    fit_func.SetLineWidth(3)  # Increase line width
    fit_func.SetLineStyle(9)  # Dashed line

    # Fit the histogram
    fit_result = hist.Fit(fit_func, "S")  # "S" option returns a TFitResultPtr
    
    # Print the fit results
    if fit_result:
        fit_result.Print("V")  # Print verbose fit result details
    
    # Draw the histogram and the fit
    canvas = ROOT.TCanvas("canvas", "Fit Result", 800, 400)
    # Create new canvas
    canvas.SetTopMargin(0.04)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.02)
    canvas.SetBottomMargin(0.12)

    # Add ticks on all sides
    canvas.SetTickx()
    canvas.SetTicky()

    # Add a grid
    canvas.SetGridx()
    canvas.SetGridy()

    ROOT.gStyle.SetGridColor(ROOT.kGray + 2)
    ROOT.gStyle.SetGridStyle(3)  # Dotted
    ROOT.gStyle.SetGridWidth(1)

    hist.Draw("PE")
    fit_func.Draw("same")

    # Add fit results as text on the plot
    fit_param = fit_func.GetParameter(0)  # Get the parameter of pol0
    fit_error = fit_func.GetParError(0)  # Get the error of the parameter

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.055)
    latex.DrawLatex(0.17, 0.82, "This Thesis")
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.17, 0.75, r"pp, #sqrt{s} = 13.6 TeV")
    latex.DrawLatex(0.17, 0.67, r"D^{0}#rightarrow K^{+}#pi^{#minus} and charge conj.")

    legend = ROOT.TLegend(0.6, 0.65, 0.85, 0.850)  # (x1, y1, x2, y2) in NDC
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)  # Transparent background
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)

    legend.AddEntry(hist, "Data", "pe")       # 'f' for filled band
    legend.AddEntry(fit_func, f"Pol. 0 Fit: {fit_param:.4f} #pm {fit_error:.4f}", "l")         # 'pe' for point+error

    legend.Draw()

    # Save the plot as an image
    img_filename = "/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/raw_prompt_frac_pol0_swapped.pdf"
    root_filename = "/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/cut_variation/raw_prompt_frac_pol0_swapped.root"
    canvas.SaveAs(img_filename)
    outfile = ROOT.TFile.Open(root_filename, "UPDATE")
    outfile.cd()
    canvas.Write()
    # Create a histogram to store fit_param and its error
    hFitParam = ROOT.TH1F("hFitParam", "Fit Parameter;Parameter Value;Counts", 1, 0, 1)
    hFitParam.SetBinContent(1, fit_param)
    hFitParam.SetBinError(1, fit_error)  # Store the uncertainty

    # Write to the ROOT file
    hFitParam.Write()
    outfile.Close()  # Close the file after writing

    # Close the file
    file.Close()

# Run the function in the notebook
test = fit_histogram()
