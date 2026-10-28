import ROOT
from ROOT import TFile, TH1F
import math

verbose = True  # Enable verbose to see what's happening

def compute_lumi_triggered():

    file = TFile.Open("~/analysis/double-D/InvariantMassFitter2D/data/AnalysisResults/AnalysisResults_data_2026.root", "READ")
    if not file or file.IsZombie():
        print(">> ERROR File not well readout")
        return
    
    print("Available directories and objects:")
    for key in file.GetListOfKeys():
        obj = key.ReadObj()
        print(f"  {key.GetName()} - Type: {obj.ClassName()}")

    dir_root = file.Get("hf-candidate-creator-2prong/Zorro")
    if not dir_root:
        print("ERROR: Directory 'hf-candidate-creator-2prong/Zorro' not found!")
        return
    
    print(f"\nContents of Zorro directory:")
    for key in dir_root.GetListOfKeys():
        obj = key.ReadObj()
        print(f"  {key.GetName()} - Type: {obj.ClassName()}")

    directoryValues = {}

    def process_dir(dir, path=""):
        values = {"analysedTriggers": -1, "inspectedTVX": -1, "scaler": -1}
        
        print(f"\nProcessing directory: '{path}'")
        for key in dir.GetListOfKeys():
            obj = key.ReadObj()
            print(f"  Found: {key.GetName()} - Type: {obj.ClassName()}")
            
            if obj.InheritsFrom("TDirectory"):
                new_path = f"{path}/{obj.GetName()}" if path else obj.GetName()
                process_dir(obj, new_path)
            elif obj.InheritsFrom("TH1"):
                hist = obj
                name = hist.GetName()  # This has the full path
                key_name = key.GetName()  # This has just the object name
                print(f"    Histogram name: {name}, Key name: {key_name}")
                
                # Use key.GetName() instead of hist.GetName()
                if key_name == "AnalysedTriggers":
                    bin_idx = hist.GetXaxis().FindBin("fHfDoubleCharm2P")
                    values["analysedTriggers"] = hist.GetBinContent(bin_idx)
                    print(f"      Found AnalysedTriggers: {values['analysedTriggers']} (bin {bin_idx})")
                elif key_name == "InspectedTVX":
                    values["inspectedTVX"] = hist.Integral()
                    print(f"      Found InspectedTVX: {values['inspectedTVX']}")
                elif key_name == "Selections":
                    bin_idx = hist.GetXaxis().FindBin("fHfDoubleCharm2P")
                    values["scaler"] = hist.GetBinContent(bin_idx)
                    print(f"      Found Selections: {values['scaler']} (bin {bin_idx})")

        if values["analysedTriggers"] != -1 or values["inspectedTVX"] != -1 or values["scaler"] != -1:
            print(f"  Saving values for path: '{path}'")
            directoryValues[path] = values
        else:
            print(f"  No values found for path: '{path}'")

    process_dir(dir_root, "")
    
    print(f"\n\nTotal directories processed: {len(directoryValues)}")
    print(f"Directory values: {directoryValues}")

    if not directoryValues:
        print("ERROR: No data found in any directory!")
        return

    totAnalysedTriggers = 0.
    totScalers = 0.
    totInspectedTVX = 0.

    for dirPath, values in directoryValues.items():
        nEvents = (values["inspectedTVX"] * values["analysedTriggers"]) / values["scaler"]
        totAnalysedTriggers += values["analysedTriggers"]
        totScalers += values["scaler"]
        totInspectedTVX += values["inspectedTVX"]

        if verbose:
            print(f"Directory: {dirPath}")
            print(f"  AnalysedTriggersOfInterest: {values['analysedTriggers']}")
            print(f"  InspectedTVX: {values['inspectedTVX']}")
            print(f"  Selections (DoubleCharm): {values['scaler']}")
            print(f"  Number of events: {nEvents}")

    totNev = totInspectedTVX * totAnalysedTriggers / totScalers

    # old 2024, before VdM scan: xSecTVX = 59.4  # mb
    xSecTVX = 52.8  # mb, value for 2024 after VdM
    lumi = totNev / (xSecTVX * 1e9)  # pb^-1

    if verbose:
        print(f"  Total AnalysedTriggersOfInterest: {totAnalysedTriggers}")
        print(f"  Total InspectedTVX: {totInspectedTVX}")
        print(f"  Total Selections (DoubleCharm): {totScalers}")
        print(f"  Integrated Luminosity TVX BEFORE BC CORRECTION: {lumi}")

    # Apply BC correction
    dirBC = file.Get("eventselection-run3")
    hCounterTVX = dirBC.Get("luminosity/hCounterTVX")
    hCounterTVXafterBCcuts = dirBC.Get("luminosity/hCounterTVXafterBCcuts")
    hLumiTVXafterBCcuts = dirBC.Get("luminosity/hLumiTVXafterBCcutsRCT")

    countsTVX = hCounterTVX.GetEntries()
    countsTVXafterBCcuts = hCounterTVXafterBCcuts.GetEntries()

    corrLumi = lumi * (countsTVXafterBCcuts / countsTVX)  # pb^-1
    #corrLumi = hLumiTVXafterBCcuts.GetEntries(1) * 1e6 # pb^-1
    print(f"  Integrated Luminosity: {corrLumi}")

    # Save corrected luminosity
    hCorrLumi = TH1F("hCorrLumi", "Corrected Integrated Luminosity;Luminosity (pb^{-1});Counts", 1, 0, 1)
    hCorrLumi.SetBinContent(1, corrLumi)
    hCorrLumi.SetBinError(1, 0)

    outFile = TFile("~/analysis/double-D/InvariantMassFitter2D/data/Integrated_Lumi_2026.root", "RECREATE")
    hCorrLumi.Write()
    outFile.Close()


if __name__ == "__main__":
    compute_lumi_triggered()
