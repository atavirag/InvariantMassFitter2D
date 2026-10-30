// MECandidates.h
#pragma once
#include <RtypesCore.h>
#include <cstdint>

struct MECandidates {
    float pt, eta, phi, m;
    Long64_t timestamp;
    int poolBin;
    int gIndexCol;
    UChar_t type;
    float weight = 1.;
};

enum CandidateType
{
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};