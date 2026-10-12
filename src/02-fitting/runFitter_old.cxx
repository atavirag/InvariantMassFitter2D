// runFit.cxx
// Example usage of InvariantMassFitter class

#include "InvariantMassFitter.h"
#include <iostream>
#include <vector>
#include "RooMsgService.h"
#include "Math/MinimizerOptions.h"

std::string osFilename = "~/MyMacros/InvariantMassFitter/data/OS_pairs_mass_1.73_2.05_eta_0.7_pairPt_10_14.root";
std::string lsFilename = "~/MyMacros/InvariantMassFitter/data/LS_pairs_mass_1.73_2.05_eta_0.7_pairPt_10_14.root";

std::vector<double> const massRange = {1.73, 2.05};
std::vector<double> const ptRange = {1.0, 24.0};

std::string const sgnFunction = "gaussian";
std::string const bkgFunction = "exponential";
std::string const reflFunction = "doubleGaussian";

bool useReflections = false;

void runBasicFit() {
    std::cout << "=== Running Basic Fit (no efficiency correction) ===" << std::endl;
    // Disable INFO messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    // Silence Minuit2 warnings
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // Initialize fitter
    InvariantMassFitter fitter(osFilename, lsFilename);
    
    // Configure mass ranges
    fitter.setMassRange(massRange[0], massRange[1], massRange[0], massRange[1]);
    fitter.setPtRange(ptRange[0], ptRange[1]);
    
    // Select fitting functions
    fitter.setFittingFunctions(sgnFunction, bkgFunction, reflFunction);
    
    // Set output file
    fitter.setOutputFile("~/MyMacros/InvariantMassFitter/data/fit_results_eta_0.7_pairPt_10_14.root");
    
    // Load data (will apply efficiency weights automatically)
    fitter.loadData();

    fitter.performPrefit("both", useReflections);  // "OS", "LS", or "both"
    
    // Step 4: Plot prefit results (optional, for validation)
    fitter.plotPrefitResults("~/MyMacros/InvariantMassFitter/plots/prefit_results");
    
    // Build model and fit
    fitter.buildModel(useReflections);
    fitter.setPrefitParameters();
    fitter.fit("both");
    // Print results
    fitter.printFitResults();
    
    // Create plots
    fitter.plotFit("~/MyMacros/InvariantMassFitter/plots/fit_results");
    
    // Save everything to ROOT file
    fitter.saveResultsToFile();
    
    std::cout << "\nBasic fit completed!" << std::endl;
}

void runEfficiencyCorrectedFit() {
    std::cout << "\n=== Running Efficiency-Corrected Fit ===" << std::endl;
    // Disable INFO messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    // Silence Minuit2 warnings
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // Initialize fitter
    InvariantMassFitter fitter(osFilename, lsFilename);
    
    // Configure mass ranges
    fitter.setMassRange(massRange[0], massRange[1], massRange[0], massRange[1]);
    fitter.setPtRange(ptRange[0], ptRange[1]);
    
    // Select fitting functions
    fitter.setFittingFunctions(sgnFunction, bkgFunction, reflFunction);
    
    // Load efficiency map (pt vs y)
    fitter.loadEfficiencyMap("~/MyMacros/InvariantMassFitter/data/Eff_times_Acc_Map_weighted.root", "hEff_new");
    //fitter.loadEfficiencyMap("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/efficiency_map/Eff_times_Acc_Map_weighted.root", "hEfficiencyMap");
    // Set output file
    fitter.setOutputFile("~/MyMacros/InvariantMassFitter/data/fit_results_corrected_eta_0.7_pairPt_0_100.root");
    
    // Load data (will apply efficiency weights automatically)
    fitter.loadData();

    fitter.performPrefit("both", useReflections);  // "OS", "LS", or "both"
    
    // Step 4: Plot prefit results (optional, for validation)
    fitter.plotPrefitResults("~/MyMacros/InvariantMassFitter/plots/prefit_results_corrected");
    
    // Build model and fit
    fitter.buildModel(useReflections);
    fitter.setPrefitParameters();
    fitter.fit("both");
    
    // Print results
    fitter.printFitResults();
    
    // Create plots
    fitter.plotFit("~/MyMacros/InvariantMassFitter/plots/fit_results_corrected", useReflections);
    
    // Save everything to ROOT file
    fitter.saveResultsToFile(useReflections);
    
    std::cout << "\nEfficiency-corrected fit completed!" << std::endl;
}

void runFitWithReflections() {
    std::cout << "\n=== Running Fit with Reflections ===" << std::endl;
    // Disable INFO messages
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    // Silence Minuit2 warnings
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // Initialize fitter
    InvariantMassFitter fitter(osFilename, lsFilename);
    
   // Configure mass ranges
    fitter.setMassRange(massRange[0], massRange[1], massRange[0], massRange[1]);
    fitter.setPtRange(ptRange[0], ptRange[1]);
    
    // Select fitting functions
    fitter.setFittingFunctions(sgnFunction, bkgFunction, reflFunction);
    fitter.setReflOverSgn(0.1);
    
    // Optional: load efficiency map
    fitter.loadEfficiencyMap("~/MyMacros/InvariantMassFitter/data/Eff_times_Acc_Map_weighted.root", "hEfficiencyMap");
    
    fitter.setOutputFile("~/MyMacros/InvariantMassFitter/data/fit_results_refl_eta_0.7_pairPt_0_100.root");

    fitter.loadData();

    fitter.performPrefit("both", true);  // "OS", "LS", or "both"
    
    // Step 4: Plot prefit results (optional, for validation)
    fitter.plotPrefitResults("~/MyMacros/InvariantMassFitter/plots/prefit_results_refl");
    
    // Build model and fit
    fitter.buildModel(true);
    fitter.setPrefitParameters(true);
    fitter.fit("both");
    
    // Output
    fitter.printFitResults();
    fitter.plotFit("~/MyMacros/InvariantMassFitter/plots/fit_results_refl_eta_0.7", true); // true = include reflections
    fitter.saveResultsToFile(true);
    
    std::cout << "\nFit with reflections completed!" << std::endl;
}

int runFitter(int argc, char** argv) {
    
    std::cout << "InvariantMassFitter Example Usage\n" << std::endl;
    
    // Choose which example to run
    if (argc > 1) {
        std::string option = argv[1];
        
        if (option == "basic") {
            runBasicFit();
        } else if (option == "efficiency") {
            runEfficiencyCorrectedFit();
        } else if (option == "reflections") {
            runFitWithReflections();
        } else {
            std::cerr << "Unknown option: " << option << std::endl;
            std::cerr << "Available options: basic, efficiency, reflections" << std::endl;
            return 1;
        }
    } else {
        // Run all examples
        runBasicFit();
        runEfficiencyCorrectedFit();
        //runFitWithReflections();
    }
    
    std::cout << "\n=== All fits completed successfully! ===" << std::endl;
    std::cout << "\nYou can now use Python to create publication-quality plots:" << std::endl;
    std::cout << "  python plot_fit_results.py fit_results.root" << std::endl;
    
    return 0;
}
