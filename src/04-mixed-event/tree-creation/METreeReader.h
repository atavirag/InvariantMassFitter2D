// TreeReader.h
#pragma once
#include "TTree.h"
#include <tuple>
#include <vector>
#include "MECandidates.h"

class METreeReader {
public:
    METreeReader(TTree* tree1, TTree* tree2);
    std::tuple<std::vector<MECandidates>, std::vector<MECandidates>> getAllEntries(); // returns a tuple of {entriesCand1, entriesCand2}

private:
    TTree* _tree1;
    TTree* _tree2;

    // Branch variables
    float ptCand1, ptCand2, etaCand1, etaCand2, phiCand1, phiCand2;
    float mDCand1, mDCand2;
    UChar_t typeCand1, typeCand2;
    int poolBinCand1, poolBinCand2, gIndexColCand1, gIndexColCand2;
    Long64_t timestampCand1, timestampCand2;
};
