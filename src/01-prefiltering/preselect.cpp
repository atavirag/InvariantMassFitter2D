#include "TreeReader.h"
#include "TreeSelector.h"
#include <TFile.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

int preselect() {

    // --- Load config ---
    ifstream config_file("preselect_config.json");
    if (!config_file.is_open()) {
        cout << ">> ERROR: Could not open preselect_config.json" << endl;
        return 1;
    }
    json config;
    config_file >> config;

    string input_file = config["input_file"];
    string tree_name  = config["tree_name"];
    string output_dir = config["output_dir"];

    // Fixed parameters
    double pt_min    = config["fixed"]["pt_min"];
    double pt_max    = config["fixed"]["pt_max"];
    double eta_min   = config["fixed"]["eta_min"];
    double eta_max   = config["fixed"]["eta_max"];
    bool   remove_ambiguous = config["fixed"]["remove_ambiguous"];

    vector<double> mass_mins = config["mass_bins"]["min"].get<vector<double>>();
    vector<double> mass_maxs = config["mass_bins"]["max"].get<vector<double>>();


    // --- Open input file ---
    TFile* file = TFile::Open(input_file.c_str(), "read");
    if (!file || file->IsZombie()) {
        cout << ">> ERROR: Could not open file: " << input_file << endl;
        return 1;
    }

    TKey* key = (TKey*)file->GetListOfKeys()->At(0);
    if (!key) {
        cout << ">> ERROR: Key not found" << endl;
        return 1;
    }

    TDirectoryFile* dir = (TDirectoryFile*)file->Get(key->GetName());
    TTree* tree = (TTree*)dir->Get(tree_name.c_str());
    if (!tree) {
        cout << ">> ERROR: Tree '" << tree_name << "' not found in file" << endl;
        return 1;
    }

    TreeReader reader(tree);
    auto rawEntries = reader.getAllEntries();
    reader.countPairEntries();

    // --- Loop over all combinations of mass min and max ---
    vector<double> bins = config["pair_pt_bins"].get<vector<double>>();
    for (double mass_min : mass_mins) {
        for (double mass_max : mass_maxs) {

        cout << "\n>> Mass window: [" << mass_min << ", " << mass_max << "]" << endl;

        // --- Loop over pt bins ---
        for (size_t i = 0; i + 1 < bins.size(); i++) {
            double pair_pt_min = bins[i];
            double pair_pt_max = bins[i + 1];

            // Format labels for output filenames
            auto fmt = [](double v) {
                ostringstream ss;
                ss << (int)v;
                return ss.str();
            };
            string pt_label   = fmt(pair_pt_min) + "_" + fmt(pair_pt_max);
            string mass_label = to_string(mass_min).substr(0, 4) + "_" + to_string(mass_max).substr(0, 4);
            string eta_label  = to_string(eta_max).substr(0, 3);

            string output_os = output_dir + "/OS_pairs_mass_" + mass_label + "_eta_" + eta_label + "_pairPt_" + pt_label + ".root";
            string output_ls = output_dir + "/LS_pairs_mass_" + mass_label + "_eta_" + eta_label + "_pairPt_" + pt_label + ".root";

            cout << "   >> pairPt [" << pair_pt_min << ", " << pair_pt_max << "]" << endl;
            cout << "      OS -> " << output_os << endl;
            cout << "      LS -> " << output_ls << endl;

            TreeSelector selector(rawEntries);
            selector.applyCuts(pt_min, pt_max, eta_min, eta_max, remove_ambiguous);
            selector.applyPairPtCut(pair_pt_min, pair_pt_max);
            selector.applyMassCutsAndDuplicate(mass_min, mass_max);
            selector.saveTreesOSvsLS(output_os.c_str(), output_ls.c_str());
        }  // end pt bin loop
        }  // end mass_max loop
    }  // end mass_min loop

    cout << "\n>> Done! All bins processed." << endl;
    return 0;
}