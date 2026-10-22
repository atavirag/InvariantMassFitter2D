// TreeReader.cpp
#include "TreeReader.h"
#include <iostream>

TreeReader::TreeReader(TTree* tree) : _tree(tree) {
    if (!_tree) {
        std::cerr << "ERROR: tree not found!" << std::endl;
        return;
    }

    _tree->SetBranchAddress("fPtCand1", &ptCand1);
    _tree->SetBranchAddress("fPtCand2", &ptCand2);
    _tree->SetBranchAddress("fYCand1", &yCand1);
    _tree->SetBranchAddress("fYCand2", &yCand2);
    _tree->SetBranchAddress("fPhiCand1", &phiCand1);
    _tree->SetBranchAddress("fPhiCand2", &phiCand2);
    _tree->SetBranchAddress("fMDCand1", &mDCand1);
    _tree->SetBranchAddress("fMDCand2", &mDCand2);
    _tree->SetBranchAddress("fMDbarCand1", &mDbarCand1);
    _tree->SetBranchAddress("fMDbarCand2", &mDbarCand2);
    _tree->SetBranchAddress("fCandidateType1", &typeCand1);
    _tree->SetBranchAddress("fCandidateType2", &typeCand2);
    _tree->SetBranchAddress("fPairType", &typePair);
}

std::vector<CandidatePair> TreeReader::getAllEntries() {
    std::vector<CandidatePair> entries;

    if (!_tree) return entries;

    Long64_t n = _tree->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        _tree->GetEntry(i);
        CandidatePair pair;
        pair.pt1 = ptCand1;  pair.pt2 = ptCand2;
        pair.y1 = yCand1;    pair.y2 = yCand2;
        pair.phi1 = phiCand1; pair.phi2 = phiCand2;
        pair.m1 = mDCand1;   pair.m2 = mDCand2;
        pair.mbar1 = mDbarCand1;   pair.mbar2 = mDbarCand2;
        pair.type1 = typeCand1; pair.type2 = typeCand2;
        pair.pairType = typePair;
        entries.push_back(pair);
    }

    return entries;
}

void TreeReader::countPairEntries() {
    int nDD = 0, nDbarDbar = 0, nDDbar = 0, nDbarD = 0, nDDbarAll = 0;
    for (Long64_t i = 0; i < _tree->GetEntries(); i++) {
        _tree->GetEntry(i);

        if (TESTBIT(typePair, DD))
            nDD++;
        if (TESTBIT(typePair, DbarDbar))
            nDbarDbar++;
        if (TESTBIT(typePair, DDbar))
            nDDbar++;
        if (TESTBIT(typePair, DbarD))
            nDbarD++;
        if (TESTBIT(typePair, DDbar) | TESTBIT(typePair, DbarD))
            nDDbarAll++;
    }
    cout << " A total of " << _tree->GetEntries() << " pairs analyzed " << endl;
    cout << "   nDD: " << nDD << ", nDbarDbar: " << nDbarDbar << "\n   nDDbar: " << nDDbar 
         << ", nDbarD: " << nDbarD << ", nDDbarAll: " << nDDbarAll << endl;
}
