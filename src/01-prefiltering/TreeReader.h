// TreeReader.h
#pragma once
#include "TTree.h"
#include <vector>
#include "CandidatePair.h"

class TreeReader {
public:
    TreeReader(TTree* tree);
    std::vector<CandidatePair> getAllEntries();
    void countPairEntries();

private:
    TTree* _tree;

    // Branch variables
    float ptCand1, ptCand2, yCand1, yCand2, phiCand1, phiCand2;
    float mDCand1, mDCand2, mDbarCand1, mDbarCand2;
    UChar_t typeCand1, typeCand2, typePair;
};
