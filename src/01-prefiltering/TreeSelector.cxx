// TreeSelector.cpp
#include "TreeSelector.h"
#include "Math/Vector4D.h"
#include <bitset>
#include <iostream>

TreeSelector::TreeSelector(const std::vector<CandidatePair>& rawEntries) : _rawEntries(rawEntries) {}

void TreeSelector::applyCuts(double ptMin, double ptMax, double yMin, double yMax, bool removeAmbiguous) {
    _selectedEntries.clear();
    int counterBeforeAmbiguous = 0, counterAfterAmbiguous = 0;

    for (auto &pair : _rawEntries) {
        if (pair.pt1 < ptMin || pair.pt1 > ptMax) continue;
        if (pair.pt2 < ptMin || pair.pt2 > ptMax) continue;
        if (pair.y1 < yMin || pair.y1 > yMax) continue;
        if (pair.y2 < yMin || pair.y2 > yMax) continue;

        counterBeforeAmbiguous++;
        // Add cuts to avoid ambiguous candidates
        if (removeAmbiguous && (TESTBIT(pair.type1, SelectedD) && TESTBIT(pair.type1, SelectedDbar))) {
            continue;
        }
        if (removeAmbiguous && (TESTBIT(pair.type2, SelectedD) && TESTBIT(pair.type2, SelectedDbar))) {
            continue;
        }
        counterAfterAmbiguous++;

        _selectedEntries.push_back(pair);
    }

    std::cout << "Selected " << _selectedEntries.size() << " entries out of " << _rawEntries.size() << std::endl;
    std::cout << "Counter before ambiguous: " << counterBeforeAmbiguous
              << "; Counter after ambiguous: " << counterAfterAmbiguous << std::endl;
}

void TreeSelector::applyMassCutsAndDuplicate(double massMin, double massMax) {
    // Treat candidates which fall into several pair types as separate events
    std::vector<CandidatePair> expandedEntries;
    int nRejected = 0;

    for (const auto& pair : _selectedEntries) {
        bool added = false;

        // DD
        if (TESTBIT(pair.pairType, DD)) {
            if ((pair.m1 >= massMin && pair.m1 <= massMax) &&
                (pair.m2 >= massMin && pair.m2 <= massMax)) {
                CandidatePair newPair = pair;
                newPair.pairType = 0; 
                expandedEntries.push_back(newPair);
                SETBIT(expandedEntries.back().pairType, DD);
                added = true;
            }
        }

        // DbarDbar
        if (TESTBIT(pair.pairType, DbarDbar)) {
            if ((pair.mbar1 >= massMin && pair.mbar1 <= massMax) &&
                (pair.mbar2 >= massMin && pair.mbar2 <= massMax)) {
                CandidatePair newPair = pair;
                newPair.pairType = 0; 
                expandedEntries.push_back(newPair);
                SETBIT(expandedEntries.back().pairType, DbarDbar);
                added = true;
            }
        }

        // DDbar
        if (TESTBIT(pair.pairType, DDbar)) {
            if ((pair.m1 >= massMin && pair.m1 <= massMax) &&
                (pair.mbar2 >= massMin && pair.mbar2 <= massMax)) {
                CandidatePair newPair = pair;
                newPair.pairType = 0; 
                expandedEntries.push_back(newPair);
                SETBIT(expandedEntries.back().pairType, DDbar);
                added = true;
            }
        }

        // DbarD
        if (TESTBIT(pair.pairType, DbarD)) {
            if ((pair.mbar1 >= massMin && pair.mbar1 <= massMax) &&
                (pair.m2 >= massMin && pair.m2 <= massMax)) {
                CandidatePair newPair = pair;
                newPair.pairType = 0; 
                expandedEntries.push_back(newPair);
                SETBIT(expandedEntries.back().pairType, DbarD);
                added = true;
            }
        }

        if (!added) {
            nRejected++;
        }
    }

    _selectedEntries = std::move(expandedEntries);
    std::cout << "Number of rejected events: " << nRejected 
              << "; Number of accepted events: " << _selectedEntries.size() << std::endl;
}

void TreeSelector::applyPairPtCut(double pairPtMin, double pairPtMax) {
    std::vector<CandidatePair> filtered;
    double pairPt = 0.0;
    
    for (const auto& entry : _selectedEntries) {
        ROOT::Math::PxPyPzMVector lv1, lv2;
        if (TESTBIT(entry.pairType, DD)) {
            lv1 = createLorentzVector(entry.phi1, entry.y1, entry.pt1, entry.m1);
            lv2 = createLorentzVector(entry.phi2, entry.y2, entry.pt2, entry.m2);
        }
        else if (TESTBIT(entry.pairType, DbarDbar)) {
            lv1 = createLorentzVector(entry.phi1, entry.y1, entry.pt1, entry.mbar1);
            lv2 = createLorentzVector(entry.phi2, entry.y2, entry.pt2, entry.mbar2);
        }
        else if (TESTBIT(entry.pairType, DDbar)) {
            lv1 = createLorentzVector(entry.phi1, entry.y1, entry.pt1, entry.m1);
            lv2 = createLorentzVector(entry.phi2, entry.y2, entry.pt2, entry.mbar2);
        }
        else if (TESTBIT(entry.pairType, DbarD)) {
            lv1 = createLorentzVector(entry.phi1, entry.y1, entry.pt1, entry.mbar1);
            lv2 = createLorentzVector(entry.phi2, entry.y2, entry.pt2, entry.m2);
        }

        pairPt = (lv1 + lv2).Pt();
        if (pairPt >= pairPtMin && pairPt <= pairPtMax) {
            filtered.push_back(entry);
        }
    }
    _selectedEntries = std::move(filtered);
}

// Create Lorentz vector to calculate the pt of the pair
ROOT::Math::PxPyPzMVector TreeSelector::createLorentzVector(double const &phi, double const &y, double const &pt, double const &m)
{
  // Calculate momentum components
  double const px = pt * cos(phi);
  double const py = pt * sin(phi);
  double const pz = pt * sinh(y);

  // Create the Lorentz vector
  ROOT::Math::PxPyPzMVector lorentzVector(px, py, pz, m);

  return lorentzVector;
}

void TreeSelector::saveTreesOSvsLS(const std::string& osFilename, const std::string& lsFilename) {
    // Open files
    TFile osFile(osFilename.c_str(), "RECREATE");
    TFile lsFile(lsFilename.c_str(), "RECREATE");

    // Create trees
    osFile.cd();
    TTree osTree("pairs", "OS candidate pairs");

    lsFile.cd();
    TTree lsTree("pairs", "LS candidate pairs");


    // Define a temporary CandidatePair object for branches
    CandidatePair osEntry, lsEntry;
    double pairPt = 0, pairM = 0;

    // Common branches
    auto createBranches = [&](TTree& tree, CandidatePair& entry) {
        tree.Branch("pt1", &entry.pt1);
        tree.Branch("pt2", &entry.pt2);
        tree.Branch("y1", &entry.y1);
        tree.Branch("y2", &entry.y2);
        tree.Branch("phi1", &entry.phi1);
        tree.Branch("phi2", &entry.phi2);
        tree.Branch("m1", &entry.m1);
        tree.Branch("m2", &entry.m2);
        tree.Branch("mbar1", &entry.mbar1);
        tree.Branch("mbar2", &entry.mbar2);
        tree.Branch("type1", &entry.type1);
        tree.Branch("type2", &entry.type2);
        tree.Branch("pairType", &entry.pairType);
        tree.Branch("pairPt", &pairPt); // new branch for the pair transverse momentum
        tree.Branch("pairMass", &pairM); // new branch for the pair invariant mass
    };

    createBranches(osTree, osEntry);
    createBranches(lsTree, lsEntry);

    // Loop over selected entries
    for (const auto& pair : _selectedEntries) {

        // Compute pairPt
        ROOT::Math::PxPyPzMVector lv1, lv2;
        if (TESTBIT(pair.pairType, DD)) {
            lv1 = createLorentzVector(pair.phi1, pair.y1, pair.pt1, pair.m1);
            lv2 = createLorentzVector(pair.phi2, pair.y2, pair.pt2, pair.m2);
        }
        else if (TESTBIT(pair.pairType, DbarDbar)) {
            lv1 = createLorentzVector(pair.phi1, pair.y1, pair.pt1, pair.mbar1);
            lv2 = createLorentzVector(pair.phi2, pair.y2, pair.pt2, pair.mbar2);
        }
        else if (TESTBIT(pair.pairType, DDbar)) {
            lv1 = createLorentzVector(pair.phi1, pair.y1, pair.pt1, pair.m1);
            lv2 = createLorentzVector(pair.phi2, pair.y2, pair.pt2, pair.mbar2);
        }
        else if (TESTBIT(pair.pairType, DbarD)) {
            lv1 = createLorentzVector(pair.phi1, pair.y1, pair.pt1, pair.mbar1);
            lv2 = createLorentzVector(pair.phi2, pair.y2, pair.pt2, pair.m2);
        }

        pairPt = (lv1 + lv2).Pt();
        pairM = (lv1 + lv2).M();

        bool isDD = TESTBIT(pair.pairType, DD);
        bool isDbarDbar = TESTBIT(pair.pairType, DbarDbar);
        bool isDDbar = TESTBIT(pair.pairType, DDbar);
        bool isDbarD = TESTBIT(pair.pairType, DbarD);

        //std::cout << "pairType=" << (int)pair.pairType 
        //  << " DD=" << isDD << " DbarDbar=" << isDbarDbar 
        //  << " DDbar=" << isDDbar << " DbarD=" << isDbarD << std::endl;

        // Determine OS or LS
        if (isDDbar || isDbarD) {
            osEntry = pair;

            //std::cout << "Filling OS: pt1=" << osEntry.pt1 << " pt2=" << osEntry.pt2 
            //            << " pairPt=" << pairPt << " pairType=" << (int)pair.pairType << std::endl;
            osTree.Fill();
        }
        if (isDD || isDbarDbar) {
            lsEntry = pair;
            //std::cout << "Filling LS: phi1=" << lsEntry.phi1 << " pt2=" << lsEntry.pt2 
            //            << " pairPt=" << pairPt << " pairType=" << (int)pair.pairType << std::endl;
            lsTree.Fill();
        }
    }

    std::cout << "OS tree entries: " << osTree.GetEntries() 
          << " LS tree entries: " << lsTree.GetEntries() 
          << " Total selected: " << _selectedEntries.size() << std::endl;

    lsFile.cd();
    lsTree.Write();
    lsFile.Close();

    osFile.cd();
    osTree.Write();
    osFile.Close();


    std::cout << "Saved OS pairs to " << osFilename << ", LS pairs to " << lsFilename << std::endl;
}