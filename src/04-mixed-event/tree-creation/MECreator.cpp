#include "METreeReader.h"
#include "MESelector.h"
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

int MECreator() {

    // --- Load config ---
    ifstream config_file("tree_creator_config.json");
    if (!config_file.is_open()) {
        cout << ">> ERROR: Could not open tree_creator_config.json" << endl;
        return 1;
    }
    json config;
    config_file >> config;

    string input_file = config["input_file"];
    string tree_name_cand1  = config["tree_name_cand1"];
    string tree_name_cand2  = config["tree_name_cand2"];
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
    TTree* tree_cand1 = (TTree*)dir->Get(tree_name_cand1.c_str());
    if (!tree_cand1) {
        cout << ">> ERROR: Tree '" << tree_name_cand1 << "' not found in file" << endl;
        return 1;
    }
    TTree* tree_cand2 = (TTree*)dir->Get(tree_name_cand2.c_str());
    if (!tree_cand2) {
        cout << ">> ERROR: Tree '" << tree_name_cand2 << "' not found in file" << endl;
        return 1;
    }

    METreeReader meTreeReader(tree_cand1, tree_cand2);
    auto rawEntries = meTreeReader.getAllEntries();

    // --- Loop over all combinations of mass min and max ---
    for (double mass_min : mass_mins) {
        for (double mass_max : mass_maxs) {

        cout << "\n>> Mass window: [" << mass_min << ", " << mass_max << "]" << endl;

        string mass_label = to_string(mass_min).substr(0, 4) + "_" + to_string(mass_max).substr(0, 4);
        string eta_label  = to_string(eta_max).substr(0, 3);

        string outputname = output_dir + "/me_tree_" + mass_label + "_eta_" + eta_label + ".root";

        cout << "      Saved as -> " << outputname << endl;

        int mix_depth = config["fixed"].value("mix_depth", 10);   // new optional JSON entry

        MESelector selector(rawEntries);
        selector.applyCuts(pt_min, pt_max, eta_min, eta_max, mass_min, mass_max, remove_ambiguous);
        selector.deduplicate();
        selector.buildMixedPairs(mix_depth);
        selector.saveOutputTrees(outputname);

        }  // end mass_max loop
    }  // end mass_min loop

    cout << "\n>> Done! All bins processed." << endl;
    return 0;
}