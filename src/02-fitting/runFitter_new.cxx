// runFitter.cxx

#include "DatasetCreator.h"
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

string fmtMass(double v) { return to_string(v).substr(0, 4); }
string fmtPt  (double v) { ostringstream ss; ss << (int)v; return ss.str(); }
string fmtEta (double v) { return to_string(v).substr(0, 3); }

string binLabel(double mass_min, double mass_max,
                double eta, double pt_min, double pt_max) {
    return "mass_"    + fmtMass(mass_min) + "_" + fmtMass(mass_max)
         + "_eta_"    + fmtEta(eta)
         + "_pairPt_" + fmtPt(pt_min) + "_" + fmtPt(pt_max);
}

// ── Stage 1 ───────────────────────────────────────────────────────────────────
// Build and save the dataset for one (period, bin, mode) combination.
// Returns the path of the saved file, or "" on failure.
string createDatasetFile(const json& cfg,
                          const json& period,
                          double mass_min, double mass_max,
                          double pt_min,   double pt_max,
                          const string& mode)
{
    double eta         = cfg["eta"];
    string input_dir   = period["input_dir"];   // per-period directory
    string output_dir  = cfg["output_dir"];
    string period_name = period["name"];
    string label       = binLabel(mass_min, mass_max, eta, pt_min, pt_max);
    string bkgFunction = cfg["fit_functions"]["background"];

    // Same filenames as before, different directory
    string os_file = input_dir + "/OS_pairs_" + label + ".root";
    string ls_file = input_dir + "/LS_pairs_" + label + ".root";

    string mode_tag;
    string eff_file, eff_hist;
    bool   useEff = false;

    if (mode == "basic") {
        mode_tag = "";
    } else if (mode == "efficiency") {
        mode_tag = "corrected";
        eff_file = cfg["modes"]["efficiency"]["efficiency_map"]["file"];
        eff_hist = cfg["modes"]["efficiency"]["efficiency_map"]["hist"];
        useEff   = true;
    } else if (mode == "reflections") {
        mode_tag = "refl";
        eff_file = cfg["modes"]["reflections"]["efficiency_map"]["file"];
        eff_hist = cfg["modes"]["reflections"]["efficiency_map"]["hist"];
        useEff   = true;
    }

    string suffix  = mode_tag.empty() ? "" : "_" + mode_tag;
    string ds_file = output_dir + "/datasets_" + bkgFunction
                   + suffix + "_" + period_name + "_" + label + ".root";

    cout << "    [Stage 1] " << period_name << " -> " << ds_file << endl;

    DatasetCreator dc(os_file, ls_file);
    dc.setMassRange(mass_min, mass_max, mass_min, mass_max);
    dc.setPtRange(cfg["pt_range"][0], cfg["pt_range"][1]);
    if (useEff) dc.loadEfficiencyMap(eff_file, eff_hist);
    dc.createDatasets();
    dc.saveDatasets(ds_file);

    return ds_file;
}
// ── Stage 2 ───────────────────────────────────────────────────────────────────
// Load (and merge) all period datasets, then fit.
void fitOneBin(const json& cfg,
               double mass_min, double mass_max,
               double pt_min,   double pt_max,
               const string& mode,
               const vector<string>& ds_files)   // one entry per period
{
    double eta        = cfg["eta"];
    string output_dir = cfg["output_dir"];
    string plots_dir  = cfg["plots_dir"];
    string label      = binLabel(mass_min, mass_max, eta, pt_min, pt_max);

    string bkgFunction = cfg["fit_functions"]["background"];
    string sgnFunction = cfg["fit_functions"]["signal"];
    string reflFunction= cfg["fit_functions"]["reflection"];

    bool   useRefl  = (mode == "reflections");
    string mode_tag;
    if      (mode == "basic")       mode_tag = "";
    else if (mode == "efficiency")  mode_tag = "corrected";
    else if (mode == "reflections") mode_tag = "refl";

    string suffix      = mode_tag.empty() ? "" : "_" + mode_tag;
    string output_file = output_dir + "/fit_results_" + bkgFunction
                       + suffix + "_" + label + ".root";
    string plot_prefit = plots_dir + "/prefit_results/prefit_results"
                       + suffix + "_" + label;
    string plot_fit    = plots_dir + "/fit_results/fit_results"
                       + suffix + "_" + label;

    cout << "    [Stage 2] Fitting " << ds_files.size() << " period(s)"
         << " -> " << output_file << endl;

    InvariantMassFitter fitter;
    fitter.setMassRange(mass_min, mass_max, mass_min, mass_max);
    fitter.setFittingFunctions(sgnFunction, bkgFunction, reflFunction);
    fitter.setOutputFile(output_file);
    if (useRefl)
        fitter.setReflOverSgn(cfg["modes"]["reflections"]["refl_over_sgn"]);

    // Load first period, merge the rest
    fitter.loadDatasets(ds_files[0]);
    for (size_t i = 1; i < ds_files.size(); ++i)
        fitter.mergeDatasets(ds_files[i]);

    fitter.buildModel(useRefl);
    fitter.performPrefit("both", useRefl);
    fitter.plotPrefitResults(plot_prefit);
    fitter.setPrefitParameters(useRefl);
    fitter.fit("both");
    fitter.printFitResults();
    fitter.plotFit(plot_fit, useRefl);
    fitter.saveResultsToFile(useRefl);
}

// ── Main ──────────────────────────────────────────────────────────────────────
int runFitter_new(int argc, char** argv) {
    cout << "InvariantMassFitter\n" << endl;
    silenceLogs();

    string config_path = (argc > 1) ? argv[1] : "fitter_config.json";
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
    auto& periods = cfg["periods"];

    vector<string> enabledModes;
    for (auto const& m : {"basic", "efficiency", "reflections"})
        if (cfg["modes"][m]["enabled"]) enabledModes.push_back(m);

    for (double mass_min : mass_mins) {
        for (double mass_max : mass_maxs) {
            for (size_t i = 0; i + 1 < pt_bins.size(); ++i) {
                double pt_min = pt_bins[i];
                double pt_max = pt_bins[i + 1];

                cout << "\n>> mass [" << mass_min << ", " << mass_max << "]"
                     << "  pairPt [" << pt_min << ", " << pt_max << "]" << endl;

                for (auto const& mode : enabledModes) {
                    cout << "   [" << mode << "]" << endl;

                    // Stage 1: create one dataset file per period
                    vector<string> ds_files;
                    for (auto const& period : periods) {
                        string ds = createDatasetFile(cfg, period,
                                                      mass_min, mass_max,
                                                      pt_min, pt_max, mode);
                        if (ds.empty()) {
                            cerr << "  Dataset creation failed for period ["
                                 << period << "], mode [" << mode << "]. Skipping bin.\n";
                            ds_files.clear();
                            break;
                        }
                        ds_files.push_back(ds);
                    }

                    if (ds_files.empty()) continue;

                    // Stage 2: merge all periods and fit
                    fitOneBin(cfg, mass_min, mass_max,
                              pt_min, pt_max, mode, ds_files);
                }
            }
        }
    }

    cout << "\n=== All fits completed successfully! ===" << endl;
    return 0;
}