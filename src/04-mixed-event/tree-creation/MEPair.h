#pragma once
#include <RtypesCore.h>

struct MEPair {
    float pt1, eta1, phi1, m1, y1;
    float pt2, eta2, phi2, m2, y2;
    UChar_t type1, type2;
    int poolBin;
    float deltaEta, deltaPhi;   // candidate 1 minus candidate 2, phi in [-pi/2, 3pi/2)
};