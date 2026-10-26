// runSPlot.cxx
// Reads splot_config.json and, for every (mass window, pair-pt bin), runs sPlot on the raw fit
// results written by runFitter. Raw distributions only: no efficiency correction.

#include "SPlotter.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "RooMsgService.h"
#include "TSystem.h"

using namespace std;
using json = nlohmann::json;

namespace {

// Same naming conventions as runFitter.cxx
string fmtMass(double v) { return to_string(v).substr(0, 4); }
string fmtPt(double v)   { ostringstream ss; ss << (int)v; return ss.str(); }
string fmtEta(double v)  { return to_string(v).substr(0, 3); }

string binLabel(double mass_min, double mass_max, double eta, double pt_min, double pt_max) {
    return "mass_" + fmtMass(mass_min) + "_" + fmtMass(mass_max)
         + "_eta_"    + fmtEta(eta)
         + "_pairPt_" + fmtPt(pt_min) + "_" + fmtPt(pt_max);
}

bool runOneBin(const json& cfg, double mass_min, double mass_max, double pt_min, double pt_max) {
    const double eta   = cfg.at("eta").get<double>();
    const string bkg   = cfg.at("background").get<string>();
    const string sfx   = cfg.value("fit_suffix", string(""));   // e.g. "" or "_corrected"
    const string label = binLabel(mass_min, mass_max, eta, pt_min, pt_max);

    const string in_file = cfg.at("fit_results_dir").get<string>()
                         + "/fit_results_" + bkg + sfx + "_" + label + ".root";
    if (gSystem->AccessPathName(in_file.c_str())) {   // returns true if NOT accessible
        cerr << "   Skipped: fit results not found: " << in_file << endl;
        return false;
    }

    const string out_dir = cfg.at("output_dir").get<string>();
    gSystem->mkdir(out_dir.c_str(), true);
    const string out_file = out_dir + "/splot_raw_" + bkg + sfx + "_" + label + ".root";

    SPlotter sp(in_file);
    for (const auto& v : cfg.at("variables"))
        sp.addVariable(v.at("name").get<string>(), v.at("bins").get<int>(),
                       v.at("min").get<double>(),  v.at("max").get<double>());

    const bool okOS = sp.run("os");
    const bool okLS = sp.run("ls");
    if (!okOS && !okLS) {
        cerr << "   sPlot failed for both OS and LS, nothing saved" << endl;
        return false;
    }

    sp.save(out_file);
    cout << "   Output: " << out_file << endl;
    return true;
}

} // namespace

int runSPlot(int argc, char** argv) {
    cout << "SPlotter\n" << endl;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    string config_path = "splot_config.json";
    if (argc > 1) config_path = argv[1];

    ifstream config_file(config_path);
    if (!config_file.is_open()) {
        cerr << ">> ERROR: Could not open " << config_path << endl;
        return 1;
    }

    try {
        json cfg;
        config_file >> cfg;

        vector<double> pt_bins   = cfg.at("pair_pt_bins").get<vector<double>>();
        vector<double> mass_mins = cfg.at("mass_bins").at("min").get<vector<double>>();
        vector<double> mass_maxs = cfg.at("mass_bins").at("max").get<vector<double>>();

        int nOk = 0, nTot = 0;
        for (double mass_min : mass_mins) {
            for (double mass_max : mass_maxs) {
                for (size_t i = 0; i + 1 < pt_bins.size(); i++) {
                    cout << "\n>> mass [" << mass_min << ", " << mass_max << "]"
                         << "  pairPt [" << pt_bins[i] << ", " << pt_bins[i + 1] << "]" << endl;
                    nTot++;
                    if (runOneBin(cfg, mass_min, mass_max, pt_bins[i], pt_bins[i + 1])) nOk++;
                }
            }
        }

        cout << "\n=== sPlot finished: " << nOk << " / " << nTot << " bins ===" << endl;
        return (nOk == nTot) ? 0 : 1;
    } catch (const json::exception& e) {
        cerr << ">> ERROR in " << config_path << ": " << e.what() << endl;
        return 1;
    }
}