// runFit.cxx
// Reads configuration from fitter_config.json and runs the selected fit modes
// Loops over all combinations of mass windows and pair pt bins

#include "InvariantMassFitter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include "RooMsgService.h"
#include "Math/MinimizerOptions.h"

using namespace std;
using json = nlohmann::json;

void silenceLogs() {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
}

// Format a double to match preselect's to_string(...).substr(0,4) convention
// e.g. 1.71 -> "1.71", 2.05 -> "2.05", 2.00 -> "2.00", 0.7 -> "0.70"
// Used only for mass values in filenames
string fmtMass(double v) {
    string s = to_string(v);  // e.g. "1.710000", "2.050000", "2.000000"
    return s.substr(0, 4);    // e.g. "1.71", "2.05", "2.00"
}

// Format pt/eta values (whole numbers stay whole, e.g. 0 -> "0", 14 -> "14")
string fmtPt(double v) {
    ostringstream ss;
    ss << (int)v;
    return ss.str();
}

// Format eta to match preselect's to_string(eta).substr(0,3), e.g. 0.7 -> "0.7"
string fmtEta(double v) {
    return to_string(v).substr(0, 3);
}

// Build the label used in all filenames for a given bin, e.g. "mass_1.71_2.05_eta_0.7_pairPt_0_2"
string binLabel(double mass_min, double mass_max, double eta, double pt_min, double pt_max) {
    return "mass_" + fmtMass(mass_min) + "_" + fmtMass(mass_max)
         + "_eta_"    + fmtEta(eta)
         + "_pairPt_" + fmtPt(pt_min) + "_" + fmtPt(pt_max);
}

// Returns false if the bin could not be processed (missing input, missing efficiency map, ...)
bool runOneBin(const json& cfg, double mass_min, double mass_max,
               double pt_min, double pt_max, const string& mode) {

    double eta        = cfg["eta"];
    string input_dir  = cfg["input_dir"];
    string output_dir = cfg["output_dir"];
    string plots_dir  = cfg["plots_dir"];

    string label = binLabel(mass_min, mass_max, eta, pt_min, pt_max);

    string os_file = input_dir + "/OS_pairs_" + label + ".root";
    string ls_file = input_dir + "/LS_pairs_" + label + ".root";

    InvariantMassFitter fitter(os_file, ls_file);
    if (!fitter.isValid()) {
        cerr << "  >> SKIPPING " << label << " [" << mode << "]: could not open input files/trees" << endl;
        return false;
    }

    fitter.setMassRange(mass_min, mass_max, mass_min, mass_max);
    fitter.setPtRange(cfg["pt_range"][0], cfg["pt_range"][1]);
    if (cfg.contains("min_efficiency")) fitter.setMinEfficiency(cfg["min_efficiency"].get<double>());

    std::string const sgnFunction = cfg["fit_functions"]["signal"];
    std::string const bkgFunction = cfg["fit_functions"]["background"];
    std::string const reflFunction = cfg["fit_functions"]["reflection"];

    fitter.setFittingFunctions(sgnFunction,
                               bkgFunction,
                               reflFunction);

    // "reflections" = efficiency-corrected + reflections; "reflections_raw" = raw (unweighted) + reflections
    bool useRefl = (mode == "reflections" || mode == "reflections_raw");
    string mode_tag;

    if (mode == "basic") {
        mode_tag = "";
    } else if (mode == "efficiency") {
        mode_tag = "corrected";
        if (!fitter.loadEfficiencyMap(cfg["modes"]["efficiency"]["efficiency_map"]["file"],
                                      cfg["modes"]["efficiency"]["efficiency_map"]["hist"])) {
            cerr << "  >> SKIPPING " << label << " [efficiency]: efficiency map not loaded" << endl;
            return false;
        }
    } else if (mode == "reflections") {
        mode_tag = "refl";
        fitter.setReflOverSgn(cfg["modes"]["reflections"]["refl_over_sgn"]);
        if (!fitter.loadEfficiencyMap(cfg["modes"]["reflections"]["efficiency_map"]["file"],
                                      cfg["modes"]["reflections"]["efficiency_map"]["hist"])) {
            cerr << "  >> SKIPPING " << label << " [reflections]: efficiency map not loaded" << endl;
            return false;
        }
    } else if (mode == "efficiency_hybrid") {
        // Raw (unweighted) fit; the 3-sigma yields are corrected with data-driven mean pair weights
        mode_tag = "hybrid";
        if (!fitter.loadEfficiencyMap(cfg["modes"]["efficiency_hybrid"]["efficiency_map"]["file"],
                                      cfg["modes"]["efficiency_hybrid"]["efficiency_map"]["hist"])) {
            cerr << "  >> SKIPPING " << label << " [efficiency_hybrid]: efficiency map not loaded" << endl;
            return false;
        }
        fitter.enableHybrid();
    } else if (mode == "reflections_raw") {
        // No efficiency map: raw (unweighted) fit with reflections
        mode_tag = "refl_raw";
        fitter.setReflOverSgn(cfg["modes"]["reflections_raw"].value("refl_over_sgn", 0.15f));
    }

    // Optionally fix the signal shape (mean, sigma, ...) to the raw (basic) fit of the same bin.
    // The background parameters stay free. Needs the basic mode to have run first.
    if (mode != "basic" && cfg["modes"][mode].value("shape_from_basic", false)) {
        string basic_file = output_dir + "/fit_results_" + bkgFunction + "_" + label + ".root";
        if (!fitter.fixSignalShapeFromFile(basic_file)) {
            cerr << "  >> SKIPPING " << label << " [" << mode << "]: could not read the raw fit shape from "
                 << basic_file << endl;
            return false;
        }
    }

    string suffix = mode_tag.empty() ? "" : "_" + mode_tag;
    string output_file = output_dir + "/fit_results_" + bkgFunction + suffix + "_" + label + ".root";
    string plot_prefit = plots_dir  + "/prefit_results/prefit_results" + suffix + "_" + label;
    string plot_fit    = plots_dir  + "/fit_results/fit_results"    + suffix + "_" + label;

    cout << "  Output: " << output_file << endl;

    fitter.setOutputFile(output_file);
    fitter.loadData();
    fitter.performPrefit("both", useRefl);
    fitter.plotPrefitResults(plot_prefit);
    fitter.buildModel(useRefl);
    fitter.setPrefitParameters(useRefl);
    fitter.fit("both");
    fitter.printFitResults();
    fitter.plotFit(plot_fit, useRefl);
    fitter.saveResultsToFile(useRefl);
    return true;
}

int runFitter(int argc, char** argv) {
    cout << "InvariantMassFitter\n" << endl;
    silenceLogs();

    // Load config
    string config_path = "fitter_config.json";
    if (argc > 1) config_path = argv[1];

    ifstream config_file(config_path);
    if (!config_file.is_open()) {
        cerr << ">> ERROR: Could not open " << config_path << endl;
        return 1;
    }
    json cfg;
    config_file >> cfg;

    vector<double> pt_bins   = cfg["pair_pt_bins"].get<vector<double>>();
    vector<double> mass_mins = cfg["mass_bins"]["min"].get<vector<double>>();
    vector<double> mass_maxs = cfg["mass_bins"]["max"].get<vector<double>>();

    int nFailed = 0;

    // Loop over all combinations
    for (double mass_min : mass_mins) {
        for (double mass_max : mass_maxs) {
            if (mass_min >= mass_max) continue; // meaningless window
            for (size_t i = 0; i + 1 < pt_bins.size(); i++) {
                double pt_min = pt_bins[i];
                double pt_max = pt_bins[i + 1];

                cout << "\n>> mass [" << mass_min << ", " << mass_max << "]"
                     << "  pairPt [" << pt_min << ", " << pt_max << "]" << endl;

                for (const char* mode : {"basic", "efficiency", "efficiency_hybrid", "reflections", "reflections_raw"}) {
                    if (cfg["modes"].contains(mode) && cfg["modes"][mode]["enabled"]) {
                        cout << "   [" << mode << "]" << endl;
                        if (!runOneBin(cfg, mass_min, mass_max, pt_min, pt_max, mode)) nFailed++;
                    }
                }
            }
        }
    }

    if (nFailed == 0) {
        cout << "\n=== All fits completed successfully! ===" << endl;
        return 0;
    }
    cout << "\n=== Finished with " << nFailed << " skipped/failed bin(s) ===" << endl;
    return 2;
}
