import ROOT
import os

def modify_canvas(input_file, output_file, canvas_name="cRatioMass"):
    """
    Modify a ROOT canvas by removing stats, title, adding custom text, and adjusting layout.
    
    Args:
        input_file (str): Path to input ROOT file
        output_file (str): Path to output ROOT file
        canvas_name (str): Name of the canvas to modify
    """
    
    # Open the input file
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist!")
        return False
        
    file_in = ROOT.TFile(input_file, "READ")
    if not file_in or file_in.IsZombie():
        print(f"Error: Could not open input file {input_file}")
        return False
    
    # Get the canvas
    canvas = file_in.Get(canvas_name)
    if not canvas:
        print(f"Error: Canvas {canvas_name} not found in {input_file}")
        file_in.Close()
        return False
    
    # Create a new canvas with more square dimensions and tighter margins
    new_canvas = ROOT.TCanvas(f"modified_{canvas_name}", "", 800, 800)
    new_canvas.SetLeftMargin(0.12)
    new_canvas.SetRightMargin(0.15)  # Slightly larger for color palette
    new_canvas.SetTopMargin(0.12)
    new_canvas.SetBottomMargin(0.09)
    
    # Get all primitives from original canvas
    primitives = canvas.GetListOfPrimitives()
    
    # Find the 2D histogram
    hist2d = None
    for primitive in primitives:
        if primitive.InheritsFrom("TH2"):
            hist2d = primitive
            break
    
    if not hist2d:
        print(f"Error: No 2D histogram found in canvas {canvas_name}")
        file_in.Close()
        return False
    
    # Clone the histogram to avoid issues when file is closed
    hist2d_copy = hist2d.Clone(f"modified_{hist2d.GetName()}")
    
    # Remove title and stats
    hist2d_copy.SetTitle("")
    hist2d_copy.GetXaxis().SetTitleOffset(1.2)
    hist2d_copy.GetYaxis().SetTitleOffset(1.78)
    hist2d_copy.SetStats(0)
    
    # Draw the histogram
    new_canvas.cd()
    hist2d_copy.Draw("COLZ")
    
    # Add red dashed lines at masses 1.84 and 1.88
    # Get the axis ranges to draw lines across the full plot
    x_min = hist2d_copy.GetXaxis().GetXmin()
    x_max = hist2d_copy.GetXaxis().GetXmax()
    y_min = hist2d_copy.GetYaxis().GetXmin()
    y_max = hist2d_copy.GetYaxis().GetXmax()
    
    # Create lines
    lines = []
    
    # Vertical lines at x = 1.84 and x = 1.88
    line_x1 = ROOT.TLine(1.84, y_min, 1.84, y_max)
    line_x2 = ROOT.TLine(1.88, y_min, 1.88, y_max)
    
    # Horizontal lines at y = 1.84 and y = 1.88
    line_y1 = ROOT.TLine(x_min, 1.84, x_max, 1.84)
    line_y2 = ROOT.TLine(x_min, 1.88, x_max, 1.88)
    
    # Style the lines (red, dashed)
    for line in [line_x1, line_x2, line_y1, line_y2]:
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)  # Dashed line
        line.SetLineWidth(3)
        line.Draw("same")
        lines.append(line)  # Keep reference to prevent deletion
    
    # Update the canvas
    new_canvas.Update()
    
    # Save to new file
    file_out = ROOT.TFile(output_file, "RECREATE")
    new_canvas.Write()
    file_out.Close()
    
    # Clean up
    file_in.Close()
    
    print(f"Modified canvas saved to {output_file}")
    return True

def process_correlation_files():
    """
    Process both correlation files and create modified versions.
    """
    
    # Input files
    ls_file = "/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/integrated_results/correlations_LS_expo_1.74_2.05_pt_1_24_ptPair_0_100.root"
    os_file = "/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/integrated_results/correlations_OS_expo_1.74_2.05_pt_1_24_ptPair_0_100.root"
    
    # Output files
    ls_output = "/home/andrea/MyMacros/Double-Dzero-correlations/test/comparisons/modified_correlations_LS_expo_1.72_2.05_pt_1_24_ptPair_0_100.root"
    os_output = "/home/andrea/MyMacros/Double-Dzero-correlations/test/comparisons/modified_correlations_OS_expo_1.72_2.05_pt_1_24_ptPair_0_100.root"
    
    # Process LS file
    print("Processing LS correlation file...")
    success_ls = modify_canvas(ls_file, ls_output)
    
    # Process OS file
    print("Processing OS correlation file...")
    success_os = modify_canvas(os_file, os_output)
    
    if success_ls and success_os:
        print("\nBoth files processed successfully!")
        print(f"Modified files created:")
        print(f"  - {ls_output}")
        print(f"  - {os_output}")
    else:
        print("\nSome files could not be processed. Check error messages above.")

def create_comparison_canvas():
    """
    Optional: Create a comparison canvas showing both LS and OS plots side by side.
    """
    
    ls_file = "/home/andrea/MyMacros/Double-Dzero-correlations/test/comparisons/modified_correlations_LS_expo_1.72_2.05_pt_1_24_ptPair_0_100.root"
    os_file = "/home/andrea/MyMacros/Double-Dzero-correlations/test/comparisons/modified_correlations_OS_expo_1.72_2.05_pt_1_24_ptPair_0_100.root"
    
    if not os.path.exists(ls_file) or not os.path.exists(os_file):
        print("Modified files not found. Run process_correlation_files() first.")
        return
    
    # Open files
    file_ls = ROOT.TFile(ls_file, "READ")
    file_os = ROOT.TFile(os_file, "READ")
    
    canvas_ls = file_ls.Get("modified_cRatioMass")
    canvas_os = file_os.Get("modified_cRatioMass")
    
    if not canvas_ls or not canvas_os:
        print("Could not find modified canvases")
        return
    
    # Create comparison canvas
    comparison = ROOT.TCanvas("comparison", "", 1600, 800)
    comparison.Divide(2, 1)
    
    # Draw LS
    comparison.cd(1)
    canvas_ls.DrawClonePad()
    
    # Add LS label
    latex_ls = ROOT.TLatex()
    latex_ls.SetNDC()
    latex_ls.SetTextSize(0.04)
    latex_ls.SetTextFont(62)
    latex_ls.DrawLatex(0.65, 0.91, "Like Sign (LS)")

    # Add custom LaTeX text
    # "This Thesis" text
    latex1 = ROOT.TLatex()
    latex1.SetNDC()
    latex1.SetTextSize(0.045)
    latex1.SetTextFont(42)
    latex1.SetTextAlign(11)  # Left aligned
    latex1.DrawLatex(0.13, 0.96, "This Thesis")
    
    # Physics info text
    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.035)
    latex2.SetTextFont(42)
    latex2.SetTextAlign(11)  # Left aligned
    latex2.DrawLatex(0.13, 0.91, "pp collisions, #sqrt{s} = 13.6 TeV")
    
    # Draw OS
    comparison.cd(2)
    canvas_os.DrawClonePad()
    
    # Add OS label
    latex_os = ROOT.TLatex()
    latex_os.SetNDC()
    latex_os.SetTextSize(0.04)
    latex_os.SetTextFont(62)
    latex_os.DrawLatex(0.56, 0.91, "Opposite Sign (OS)")
    
    # Save comparison
    comparison.SaveAs("/home/andrea/MyMacros/Double-Dzero-correlations/samples/comparison_LS_OS.pdf")
    comparison.SaveAs("/home/andrea/MyMacros/Double-Dzero-correlations/samples/comparison_LS_OS.png")
    
    print("Comparison canvas saved as comparison_LS_OS.pdf and comparison_LS_OS.png")
    
    # Clean up
    file_ls.Close()
    file_os.Close()

if __name__ == "__main__":
    # Set ROOT to batch mode to avoid opening graphics windows
    ROOT.gROOT.SetBatch(True)
    
    # Process the files
    process_correlation_files()
    
    # Optionally create comparison
    create_comparison_canvas()
    
    print("\nScript completed!")
