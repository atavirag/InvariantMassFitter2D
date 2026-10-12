#ifndef REFLECTED_FUNCTIONS_H
#define REFLECTED_FUNCTIONS_H

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <string>

// ---------- Reflected models ----------
RooAbsPdf* createGaussianReflected(RooWorkspace* ws, const std::string& name, RooRealVar& mass);
RooAbsPdf* createDoubleGaussianReflected(RooWorkspace* ws, const std::string& name, RooRealVar& mass);

#endif
