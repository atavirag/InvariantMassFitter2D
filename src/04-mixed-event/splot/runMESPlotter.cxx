// Runs the ME sPlot step for all mass windows, using splot_config.json
#include "MESPlotter.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

namespace {
string expandPath(const string& p) {
    TString s(p.c_str());
    gSystem->ExpandPathName(s);
    return string(s.Data());
}
string fmtMass(double v) { return to_string(v).substr(0, 4); }
string fmtEta(double v)  { return to_string(v).substr(0, 3); }

// Must match the label used for the ME fit file names in runFit.cxx
string meLabel(double mass_min, double mass_max, double eta) {
    return "mass_" + fmtMass(mass_min) + "_" + fmtMass(mass_max)
         + "_eta_" + fmtEta(eta);
}

void plotCorrelation(MESPlotter& sp, const string& tag, const string& nx, const string& ny, const string& prefix) {
    const string base = "h2_" + tag + "_" + nx + "_" + ny + "_";
    TH2D* hAll = sp.find2D(base + "All");
    TH2D* hSS  = sp.find2D(base + "SgnSgn");
    TH2D* hN   = sp.find2D(base + "SgnSgn_norm");
    if (!hAll || !hSS) return;

    gStyle->SetOptStat(0);
    TCanvas c("c_me_corr", "ME correlation", 1800, 1000);
    c.Divide(3, 2);
    c.cd(1); hAll->Draw("COLZ");
    c.cd(2); hSS->Draw("COLZ");
    c.cd(3); if (hN) hN->Draw("COLZ");
    c.cd(4);
    TH1D* px = hSS->ProjectionX((base + "px").c_str());
    px->SetTitle(("sWeighted ME (SgnSgn), " + nx + " projection").c_str());
    px->SetLineWidth(2); px->Draw("E");
    c.cd(5);
    TH1D* py = hSS->ProjectionY((base + "py").c_str());
    py->SetTitle(("sWeighted ME (SgnSgn), " + ny + " projection").c_str());
    py->SetLineWidth(2); py->Draw("E");
    c.SaveAs((prefix + "_" + nx + "_" + ny + ".pdf").c_str());
    c.SaveAs((prefix + "_" + nx + "_" + ny + ".png").c_str());
}
}  // namespace

int runMESPlotter(int argc, char** argv) {
    string config_path = "splot_config.json";
    if (argc > 1) config_path = argv[1];
    ifstream config_file(config_path);
    if (!config_file.is_open()) { cerr << ">> ERROR: Could not open " << config_path << endl; return 1; }
    json cfg;
    config_file >> cfg;

    if (!cfg.contains("me")) { cerr << ">> ERROR: no 'me' block in " << config_path << endl; return 1; }
    const json& me = cfg["me"];
    if (!me.value("enabled", true)) return 0;

    const double eta   = cfg["eta"];
    const string bkg   = cfg["background"];
    const string tree  = me.value("tree_name", string("pairs"));
    const string suffix = me.value("fit_suffix", string(""));   // "" or "_corrected" (as in the ME fit file names)
    const string input_dir = expandPath(me["input_dir"].get<string>());
    const string fit_dir   = expandPath(me["fit_results_dir"].get<string>());
    const string out_dir   = expandPath(me["output_dir"].get<string>());
    const string plots_dir = expandPath(me.value("plots_dir", out_dir + "/plots"));
    filesystem::create_directories(out_dir);
    filesystem::create_directories(plots_dir);

    vector<double> mass_mins = cfg["mass_bins"]["min"].get<vector<double>>();
    vector<double> mass_maxs = cfg["mass_bins"]["max"].get<vector<double>>();
    const double normDeta = me.value("norm_deta", 0.2), normDphi = me.value("norm_dphi", 0.4);

    int nFailed = 0;
    for (double mass_min : mass_mins) {
        for (double mass_max : mass_maxs) {
            if (mass_min >= mass_max) continue;
            const string mlabel = meLabel(mass_min, mass_max, eta);
            cout << "\n>> [ME sPlot] " << mlabel << endl;

            MESPlotter sp(input_dir + "/me_tree_" + fmtMass(mass_min) + "_" + fmtMass(mass_max)
                          + "_eta_" + fmtEta(eta) + ".root", tree);

            // Single fit for the whole sample: one "bin" covering all pT
            const string fitFile = fit_dir + "/fit_results_me_" + bkg + suffix + "_" + mlabel + ".root";
            if (!sp.addBin(0., 1e9, fitFile)) {
                cerr << "  >> SKIPPING: could not load " << fitFile << endl;
                nFailed++;
                continue;
            }

            for (const auto& v : cfg["variables"])
                sp.addVariable(v["name"].get<string>(), v["bins"].get<int>(), v["min"].get<double>(), v["max"].get<double>());
            for (const auto& v : cfg["variables2D"])
                sp.addVariable2D(v["x"]["name"].get<string>(), v["x"]["bins"].get<int>(), v["x"]["min"].get<double>(), v["x"]["max"].get<double>(),
                                 v["y"]["name"].get<string>(), v["y"]["bins"].get<int>(), v["y"]["min"].get<double>(), v["y"]["max"].get<double>());

            if (!sp.run("me")) { nFailed++; continue; }

            for (const auto& v : cfg["variables2D"]) {
                const string nx = v["x"]["name"], ny = v["y"]["name"];
                sp.addNormalized2D("h2_me_" + nx + "_" + ny + "_SgnSgn", normDeta, normDphi);
                plotCorrelation(sp, "me", nx, ny, plots_dir + "/me_splot_" + mlabel);
            }
            sp.save(out_dir + "/splot_me_" + mlabel + ".root");
        }
    }
    return nFailed == 0 ? 0 : 2;
}