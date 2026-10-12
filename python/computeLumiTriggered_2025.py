import ROOT
from ROOT import TFile, TH1F
import math

verbose = True  # Enable verbose to see what's happening

def compute_lumi_triggered():

    file = TFile.Open("~/MyMacros/InvariantMassFitter/data/AnalysisResults/AnalysisResults_2025_full.root", "READ")
    if not file or file.IsZombie():
        print(">> ERROR File not well readout")
        return
    
    print("Available directories and objects:")
    for key in file.GetListOfKeys():
        obj = key.ReadObj()
        print(f"  {key.GetName()} - Type: {obj.ClassName()}")

    dir_root = file.Get("eventselection-run3/luminosity")
    if not dir_root:
        print("ERROR: Directory 'eventselection-run3/luminosity' not found!")
        return
    
    hLumiTVXafterBCcuts = dir_root.Get("hLumiTVXafterBCcuts")

    # Get lumi after BC cuts
    #bin_idx = hLumiTVXafterBCcuts.GetXaxis().FindBin("fHfDoubleCharm2P")
    lumi = hLumiTVXafterBCcuts.GetEntries()*10**(-6)
    print(f"      Lumi after BC cuts: {lumi} ")

    # Save corrected luminosity
    hCorrLumi = TH1F("hCorrLumi", "Corrected Integrated Luminosity;Luminosity (pb^{-1});Counts", 1, 0, 1)
    hCorrLumi.SetBinContent(1, lumi)
    hCorrLumi.SetBinError(1, 0)

    outFile = TFile("/home/andrea/MyMacros/InvariantMassFitter/test/Integrated_Lumi_2025.root", "UPDATE")
    hCorrLumi.Write()
    outFile.Close()


if __name__ == "__main__":
    compute_lumi_triggered()