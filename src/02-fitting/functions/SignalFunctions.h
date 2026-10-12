#ifndef SIGNAL_FUNCTIONS_H
#define SIGNAL_FUNCTIONS_H

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <string>

// ---------- Signal models ----------
RooAbsPdf* createGaussianSignal(RooWorkspace* ws, const std::string& name, RooRealVar& mass);
RooAbsPdf* createDoubleGaussianSignal(RooWorkspace* ws, const std::string& name, RooRealVar& mass);
RooAbsPdf* createCrystalBallSignal(RooWorkspace* ws, const std::string& name, RooRealVar& mass);

#endif
