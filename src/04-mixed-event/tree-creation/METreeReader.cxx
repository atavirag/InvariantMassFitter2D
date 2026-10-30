// METreeReader.cpp
#include "METreeReader.h"
#include <iostream>

METreeReader::METreeReader(TTree* tree1, TTree* tree2) : _tree1(tree1), _tree2(tree2) {
    if (!_tree1 || !_tree2) {
        std::cerr << "ERROR: trees not found!" << std::endl;
        return;
    }

    _tree1->Print();

    _tree1->SetBranchAddress("fPtD", &ptCand1);
    _tree1->SetBranchAddress("fEtaD", &etaCand1);
    _tree1->SetBranchAddress("fPhiD", &phiCand1);
    _tree1->SetBranchAddress("fMD", &mDCand1);
    _tree1->SetBranchAddress("fDType", &typeCand1);
    _tree1->SetBranchAddress("fPoolBinD", &poolBinCand1);
    _tree1->SetBranchAddress("fGIndexColD", &gIndexColCand1);
    _tree1->SetBranchAddress("fTimestampD", &timestampCand1);

    _tree2->Print();

    _tree2->SetBranchAddress("fPtD", &ptCand2);
    _tree2->SetBranchAddress("fEtaD", &etaCand2);
    _tree2->SetBranchAddress("fPhiD", &phiCand2);
    _tree2->SetBranchAddress("fMD", &mDCand2);
    _tree2->SetBranchAddress("fDType", &typeCand2);
    _tree2->SetBranchAddress("fPoolBinD", &poolBinCand2);
    _tree2->SetBranchAddress("fGIndexColD", &gIndexColCand2);
    _tree2->SetBranchAddress("fTimestampD", &timestampCand2);
}

std::tuple<std::vector<MECandidates>, std::vector<MECandidates>> METreeReader::getAllEntries() {
    std::vector<MECandidates> entries1;
    std::vector<MECandidates> entries2;

    if (!_tree1 || !_tree2) return { entries1, entries2 };

    Long64_t n1 = _tree1->GetEntries();
    for (Long64_t i = 0; i < n1; i++) {
        _tree1->GetEntry(i);
        MECandidates cand1;
        cand1.pt = ptCand1;
        cand1.eta = etaCand1;
        cand1.phi = phiCand1;
        cand1.m = mDCand1;
        cand1.type = typeCand1;
        cand1.poolBin = poolBinCand1;
        cand1.gIndexCol = gIndexColCand1;
        cand1.timestamp = timestampCand1;
        entries1.push_back(cand1);
    }

    Long64_t n2 = _tree2->GetEntries();
    for (Long64_t i = 0; i < n2; i++) {
        _tree2->GetEntry(i);
        MECandidates cand2;
        cand2.pt = ptCand2;
        cand2.eta = etaCand2;
        cand2.phi = phiCand2;
        cand2.m = mDCand2;
        cand2.type = typeCand2;
        cand2.poolBin = poolBinCand2;
        cand2.gIndexCol = gIndexColCand2;
        cand2.timestamp = timestampCand2;
        entries2.push_back(cand2);
    }

    return { entries1, entries2 };
}

