#ifndef BACKGROUND_FUNCTIONS_H
#define BACKGROUND_FUNCTIONS_H

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <string>

// ---------- Background models ----------
RooAbsPdf* createChebychevBackground(RooWorkspace* ws, const std::string& name, RooRealVar& mass);
RooAbsPdf* createExponentialBackground(RooWorkspace* ws, const std::string& name, RooRealVar& mass);
RooAbsPdf* createPolynomialBackground(RooWorkspace* ws, const std::string& name, RooRealVar& mass, int order);

#endif