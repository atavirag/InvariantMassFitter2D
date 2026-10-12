import ROOT
import os


def set_histo_style(h, y_title):
    h.SetMarkerStyle(21)
    h.SetMarkerSize(1.2)
    h.SetMarkerColor(ROOT.kBlue-4)
    h.SetLineColor(ROOT.kBlue-4)
    h.SetLineWidth(2)

    h.SetTitle("")
    h.GetYaxis().SetTitle(y_title)
    h.GetXaxis().SetTitle("BDT cut")
    h.GetYaxis().SetTitleSize(0.055)
    h.GetYaxis().SetTitleOffset(1.5)
    h.GetXaxis().SetTitleSize(0.055)
    h.GetYaxis().SetLabelSize(0.06)
    h.GetXaxis().SetLabelSize(0.07)

def set_canvas_style(canvas):
    ROOT.gStyle.SetOptStat(0)
    canvas.SetTopMargin(0.06)
    canvas.SetLeftMargin(0.17)
    canvas.SetRightMargin(0.03)
    canvas.SetBottomMargin(0.12)

    # Add ticks on all sides
    canvas.SetTickx()
    canvas.SetTicky()

    # Add a grid
    #canvas.SetGridx()
    #canvas.SetGridy()

    ROOT.gStyle.SetGridColor(ROOT.kGray + 3)
    ROOT.gStyle.SetGridStyle(3)  # Dotted
    ROOT.gStyle.SetGridWidth(1)

def plot_bdt_systematics():
    # Set ROOT style for better looking plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetTitleOffset(1.2, "Y")
    ROOT.gROOT.SetBatch(True)  # Run in batch mode for cleaner output
    
    # Open the ROOT file
    file = ROOT.TFile.Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/bdt_systematics_results.root", "READ")
    if not file or file.IsZombie():
        print("Error: Cannot open file bdtSystematicsResults.root")
        return
    
    # Define pt bins
    pt_min = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
    pt_max = [2.0, 3.0, 5.0, 7.0, 10.0, 24.0]
    
    # Get the graph (same for all plots)
    graph = file.Get("Graph")
    if not graph:
        print("Warning: Graph not found in file")
    else:
        graph.SetLineColor(ROOT.kBlack)
        graph.SetLineWidth(2)
        graph.SetMarkerColor(ROOT.kBlack)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.0)
    
    # Loop through pt bins and create one canvas per pt range
    for i, (pmin, pmax) in enumerate(zip(pt_min, pt_max)):
        # Create canvas for this pt bin with 2x2 division
        canvas_name = f"c_pt_{pmin}_{pmax}"
        canvas = ROOT.TCanvas(canvas_name, f"BDT Systematics: {pmin} < p_T < {pmax} GeV/c", 1200, 600)
        canvas.Divide(3, 2, 0.01, 0.01)  # 2x2 grid with small margins
        for primitive in canvas.GetListOfPrimitives():
            if isinstance(primitive, ROOT.TPad):
                set_canvas_style(primitive)

        # Construct names
        raw_yield_canvas_name = f"c_{pmin}_{pmax}"
        effs_canvas_name = f"c_effs_{pmin}_{pmax}"
        ratio_canvas_name = f"c_ratio_{pmin}_{pmax}"
        divide_canvas_name = f"c_divide_{pmin}_{pmax}"
        width_name = f"hWidthSys_{pmin}_{pmax}"
        chi2_name = f"hChi2Sys_{pmin}_{pmax}"
        
        # Get canvas and histograms
        c_effs = file.Get(effs_canvas_name)
        c_ratio = file.Get(ratio_canvas_name)
        c_divide = file.Get(divide_canvas_name)
        c_raw_yield = file.Get(raw_yield_canvas_name)
        h_width = file.Get(width_name)
        h_chi2 = file.Get(chi2_name)
        
        # Pad 1: Raw Yield Systematics (from canvas)
        canvas.cd(1)
        
        if c_raw_yield:
            # Clone all objects from the original canvas
            c_raw_yield.cd()
            primitives = c_raw_yield.GetListOfPrimitives()
            
            # Switch back to our canvas pad
            canvas.cd(1)
            
            # Find the first histogram to set up the axes properly
            first_hist = None
            for primitive in primitives:
                if primitive.ClassName() == "TH1D":
                    first_hist = primitive
                    break
            
            if first_hist:
                # Draw the first histogram to establish the frame
                first_hist_clone = first_hist.Clone()
                set_histo_style(first_hist_clone, "Raw Yield")
                first_hist_clone.Draw("PE")
                
                # Now draw all other objects
                for primitive in primitives:
                    if primitive == first_hist:
                        continue  # Skip the first histogram as it's already drawn
                    
                    primitive_clone = primitive.Clone()
                    
                    if primitive.ClassName() == "TFrame":
                        continue  # Skip TFrame
                    elif primitive.ClassName() == "TH1D":
                        set_histo_style(primitive_clone, "Raw Yield")
                        primitive_clone.Draw("SAME PE")
                    elif primitive.ClassName() == "TLine":
                        primitive_clone.SetLineWidth(2)
                        primitive_clone.Draw("SAME")
                    elif primitive.ClassName() == "TPaveText":
                        # Skip the original TPaveText, we'll add our own
                        continue
                    else:
                        primitive_clone.Draw("SAME")
            
            # Add pt range text
            text1 = ROOT.TLatex()
            text1.SetNDC()
            text1.SetTextSize(0.075)
            text1.SetTextFont(42)
            text1.DrawLatex(0.2, 0.87, f"This Thesis")
            text1.DrawLatex(0.2, 0.8, r'pp, #sqrt{#it{s}} = 13.6 TeV')

            
            # Update the pad
            ROOT.gPad.Update()
        else:
            # Create empty plot if canvas not found
            dummy1 = ROOT.TH1F("dummy1", "Raw Yield Systematics", 10, 0, 10)
            dummy1.GetXaxis().SetTitle("Bin")
            dummy1.GetYaxis().SetTitle("Systematic Uncertainty")
            dummy1.Draw()
            text1 = ROOT.TLatex()
            text1.SetNDC()
            text1.SetTextSize(0.05)
            text1.SetTextFont(42)
            text1.DrawLatex(0.3, 0.5, "No Data Available")
        
        # Pad 2: Width Systematics
        canvas.cd(2)

        if h_width:
            set_histo_style(h_width, 'Width (GeV/#it{c}^{2})')
            h_width.GetYaxis().SetRangeUser(h_width.GetMinimum()-0.003, h_width.GetMaximum()+0.003)
            h_width.GetXaxis().SetLabelSize(0.07)
            h_width.Draw("PE")

            text1 = ROOT.TLatex()
            text1.SetNDC()
            text1.SetTextSize(0.075)
            text1.SetTextFont(42)
            text1.DrawLatex(0.2, 0.85, f"{pmin} < p_{{T}} < {pmax} GeV/c")
            text1.DrawLatex(0.2, 0.75, r"D^{0}#rightarrow K^{#minus}#pi^{+} and charge conj.")

        else:
            # Create empty histogram for consistent axes
            dummy2 = ROOT.TH1F("dummy2", "Width Systematics", 10, 0, 10)
            dummy2.GetXaxis().SetTitle("Bin")
            dummy2.GetYaxis().SetTitle("Systematic Uncertainty")
            dummy2.Draw()
            text2 = ROOT.TLatex()
            text2.SetNDC()
            text2.SetTextSize(0.05)
            text2.SetTextFont(42)
            text2.DrawLatex(0.3, 0.5, "No Data Available")
        
        # Pad 3: Chi2 Systematics
        canvas.cd(3)
        
        if h_chi2:
            set_histo_style(h_chi2, '#chi^{2}/ndf')
            h_chi2.GetYaxis().SetRangeUser(0., 2.5)
            h_chi2.GetXaxis().SetLabelSize(0.07)
            h_chi2.Draw("PE")

        else:
            # Create empty histogram for consistent axes
            dummy3 = ROOT.TH1F("dummy3", "#chi^{2} Systematics", 10, 0, 10)
            dummy3.GetXaxis().SetTitle("Bin")
            dummy3.GetYaxis().SetTitle("Systematic Uncertainty")
            dummy3.Draw()
            text3 = ROOT.TLatex()
            text3.SetNDC()
            text3.SetTextSize(0.05)
            text3.SetTextFont(42)
            text3.DrawLatex(0.3, 0.5, "No Data Available")
        
        # Pad 5: Graph
        canvas.cd(5)
        
        if c_divide:
            c_divide.cd()
            primitives_divide = c_divide.GetListOfPrimitives()
            
            # Switch back to our canvas pad
            canvas.cd(5)
            
            # Find the first histogram to set up the axes properly
            first_hist_divide = None
            for primitive in primitives_divide:
                if primitive.ClassName() == "TH1D":
                    first_hist_divide = primitive
                    break
            
            if first_hist_divide:
                # Draw the first histogram to establish the frame
                first_hist_clone_divide = first_hist_divide.Clone()
                set_histo_style(first_hist_clone_divide, "#it{N}_{raw}/(Acc. #times #varepsilon)")
                first_hist_clone_divide.Draw("PE")
                
                # Now draw all other objects
                for primitive in primitives_divide:
                    if primitive == first_hist_divide:
                        continue  # Skip the first histogram as it's already drawn
                    
                    primitive_clone = primitive.Clone()
                    
                    if primitive.ClassName() == "TFrame":
                        continue  # Skip TFrame
                    elif primitive.ClassName() == "TH1D":
                        set_histo_style(primitive_clone, "(Acc. #times #varepsilon)")
                        primitive_clone.Draw("SAME PE")
                    elif primitive.ClassName() == "TLine":
                        primitive_clone.SetLineWidth(2)
                        primitive_clone.Draw("SAME")
                    elif primitive.ClassName() == "TPaveText":
                        # Skip the original TPaveText, we'll add our own
                        continue
                    else:
                        primitive_clone.Draw("SAME")

            
            # Update the pad
            ROOT.gPad.Update()
        
            """
            if graph:
            # Clone the graph to avoid modifying the original
            graph_clone = graph.Clone(f"graph_pt_{pmin}_{pmax}")
            graph_clone.SetTitle("")
            set_histo_style(graph_clone, "Systematic uncertainty")
            graph_clone.SetMarkerColor(ROOT.kRed)
            graph_clone.SetLineColor(ROOT.kRed)
            graph_clone.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
            graph_clone.GetXaxis().SetRangeUser(1.0, 24.0)
            graph_clone.GetXaxis().SetLabelSize(0.055)
            graph_clone.Draw("AP")
            """
        else:
            # Create empty graph for consistent display
            dummy_graph = ROOT.TGraph(2)
            dummy_graph.SetPoint(0, 0, 0)
            dummy_graph.SetPoint(1, 1, 1)
            dummy_graph.SetTitle("Graph")
            dummy_graph.GetXaxis().SetTitle("X")
            dummy_graph.GetYaxis().SetTitle("Y")
            dummy_graph.Draw("ALP")
            dummy_graph.SetMarkerColor(0)  # Make invisible
            dummy_graph.SetLineColor(0)
            
            text4 = ROOT.TLatex()
            text4.SetNDC()
            text4.SetTextSize(0.05)
            text4.SetTextFont(42)
            text4.DrawLatex(0.3, 0.5, "No Data Available")
        
        # Pad 1: Raw Yield Systematics (from canvas)
        canvas.cd(4)
        
        if c_effs:
            # Clone all objects from the original canvas
            c_effs.cd()
            primitives_effs = c_effs.GetListOfPrimitives()
            
            # Switch back to our canvas pad
            canvas.cd(4)
            
            # Find the first histogram to set up the axes properly
            first_hist_effs = None
            for primitive in primitives_effs:
                if primitive.ClassName() == "TH1D":
                    first_hist_effs = primitive
                    break
            
            if first_hist_effs:
                # Draw the first histogram to establish the frame
                first_hist_clone_effs = first_hist_effs.Clone()
                set_histo_style(first_hist_clone_effs, "(Acc. #times #varepsilon)")
                first_hist_clone_effs.Draw("PE")
                
                # Now draw all other objects
                for primitive in primitives_effs:
                    if primitive == first_hist_effs:
                        continue  # Skip the first histogram as it's already drawn
                    
                    primitive_clone = primitive.Clone()
                    
                    if primitive.ClassName() == "TFrame":
                        continue  # Skip TFrame
                    elif primitive.ClassName() == "TH1D":
                        set_histo_style(primitive_clone, "(Acc. #times #varepsilon)")
                        primitive_clone.Draw("SAME PE")
                    elif primitive.ClassName() == "TLine":
                        primitive_clone.SetLineWidth(2)
                        primitive_clone.Draw("SAME")
                    elif primitive.ClassName() == "TPaveText":
                        # Skip the original TPaveText, we'll add our own
                        continue
                    else:
                        primitive_clone.Draw("SAME")

            
            # Update the pad
            ROOT.gPad.Update()
        else:
            # Create empty plot if canvas not found
            dummy4 = ROOT.TH1F("dummy4", "Efficiencies", 10, 0, 10)
            dummy4.GetXaxis().SetTitle("Bin")
            dummy4.GetYaxis().SetTitle("Systematic Uncertainty")
            dummy4.Draw()
            text4 = ROOT.TLatex()
            text4.SetNDC()
            text4.SetTextSize(0.05)
            text4.SetTextFont(42)
            text4.DrawLatex(0.3, 0.5, "No Data Available")

        # Pad 1: Raw Yield Systematics (from canvas)
        canvas.cd(6)
        
        if c_ratio:
            # Clone all objects from the original canvas
            c_ratio.cd()
            primitives_ratio = c_ratio.GetListOfPrimitives()
            
            # Switch back to our canvas pad
            canvas.cd(6)
            
            # Find the first histogram to set up the axes properly
            first_hist_ratio = None
            for primitive in primitives_ratio:
                if primitive.ClassName() == "TH1D":
                    first_hist_ratio = primitive
                    break
            
            if first_hist_ratio:
                # Draw the first histogram to establish the frame
                first_hist_clone_ratio = first_hist_ratio.Clone()
                set_histo_style(first_hist_clone_ratio, "Trial/Default")
                first_hist_clone_ratio.Draw("PE")
                
                # Now draw all other objects
                for primitive in primitives_ratio:
                    if primitive == first_hist_ratio:
                        continue  # Skip the first histogram as it's already drawn
                    
                    primitive_clone = primitive.Clone()
                    
                    if primitive.ClassName() == "TFrame":
                        continue  # Skip TFrame
                    elif primitive.ClassName() == "TH1D":
                        set_histo_style(primitive_clone, "Trial/Default")
                        primitive_clone.GetXaxis().SetLabelSize(0.055)
                        primitive_clone.Draw("SAME PE")
                    elif primitive.ClassName() == "TLine":
                        primitive_clone.SetLineWidth(2)
                        primitive_clone.Draw("SAME")
                    elif primitive.ClassName() == "TPaveText":
                        # Skip the original TPaveText, we'll add our own
                        continue
                    else:
                        primitive_clone.Draw("SAME")

            
            # Update the pad
            ROOT.gPad.Update()
        else:
            # Create empty plot if canvas not found
            dummy5 = ROOT.TH1F("dummy5", "Efficiencies", 10, 0, 10)
            dummy5.GetXaxis().SetTitle("Bin")
            dummy5.GetYaxis().SetTitle("Systematic Uncertainty")
            dummy5.Draw()
            text5 = ROOT.TLatex()
            text5.SetNDC()
            text5.SetTextSize(0.05)
            text5.SetTextFont(42)
            text5.DrawLatex(0.3, 0.5, "No Data Available")
        
        # Update all pads
        for primitive in canvas.GetListOfPrimitives():
            if isinstance(primitive, ROOT.TPad):
                set_canvas_style(primitive)
        canvas.Update()
        
        # Save canvas
        output_name = f"BDT_Systematics_4pads_pt_{pmin}_{pmax}"
        canvas.SaveAs(f"/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/{output_name}.pdf")
        #canvas.SaveAs(f"{output_name}.png")
        
        print(f"Saved 4-pad plot for pt range {pmin}-{pmax} GeV/c")
        
        # Keep canvas in memory (optional)
        canvas.Draw()
    
    file.Close()
    print("All 4-pad plots created successfully!")

if __name__ == "__main__":
    plot_bdt_systematics()