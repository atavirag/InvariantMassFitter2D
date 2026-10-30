// MESelector.h
#pragma once
#include <vector>
#include <string>
#include "MECandidates.h"
#include "MEPair.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector4D.h"

class MESelector {
public:
    MESelector(const std::tuple<std::vector<MECandidates>, std::vector<MECandidates>>& rawEntries);
    void applyCuts(double ptMin, double ptMax, double etaMin, double etaMax, double massMin, double massMax, bool removeAmbiguous = true);
    void deduplicate();
    void buildMixedPairs(int mixDepth);
    void saveOutputTrees(const std::string& outputname); // stores output as a histo in a root file

private:
    std::tuple<std::vector<MECandidates>, std::vector<MECandidates>> _selectedEntries;
    std::tuple<std::vector<MECandidates>, std::vector<MECandidates>> _rawEntries;
    std::vector<MEPair> _pairs;

};