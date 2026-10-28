// CandidatePair.h
#pragma once

struct CandidatePair {
    float pt1, pt2;
    float y1, y2;
    float phi1, phi2;
    float eta1, eta2;
    float m1, m2;
    float mbar1, mbar2;
    UChar_t type1, type2, pairType;
    float weight1 = 1.;
    float weight2 = 1.;
};

enum PairTypeSel
{
  DD = 0,   // D0-D0
  DbarDbar, // D0bar-D0bar
  DDbar,
  DbarD
};

enum CandidateType
{
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};