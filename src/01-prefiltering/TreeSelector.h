// TreeSelector.h
#pragma once
#include <vector>
#include <string>
#include "CandidatePair.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector4D.h"

class TreeSelector {
public:
    TreeSelector(const std::vector<CandidatePair>& rawEntries);
    void applyCuts(double ptMin, double ptMax, double yMin, double yMax, bool removeAmbiguous = true);
    void applyMassCutsAndDuplicate(double massMin, double massMax);
    ROOT::Math::PxPyPzMVector createLorentzVector(double const &phi, double const &y, double const &pt, double const &m);
    void applyPairPtCut(double pairPtMin, double pairPtMax);
    void saveTreesOSvsLS(const std::string& osFilename, const std::string& lsFilename);

private:
    std::vector<CandidatePair> _selectedEntries;
    std::vector<CandidatePair> _rawEntries;
};