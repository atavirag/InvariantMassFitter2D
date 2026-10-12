#include "BackgroundFunctions.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooPolynomial.h"

using namespace std;

RooAbsPdf* createChebychevBackground(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar c0((name + "_c0").c_str(), "c0", 0, -1, 1);
    RooRealVar c1((name + "_c1").c_str(), "c1", 0, -1, 1);

    RooChebychev bkg((name + "_bkg").c_str(), "Chebychev Background", mass, RooArgList(c0, c1));
    ws->import(bkg);
    return ws->pdf((name + "_bkg").c_str());
}

RooAbsPdf* createExponentialBackground(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar lambda((name + "_lambda").c_str(), "lambda", -1.0, -10.0, 0.0);

    RooExponential bkg((name + "_exp").c_str(), "Exponential Background", mass, lambda);
    ws->import(bkg);
    return ws->pdf((name + "_exp").c_str());
}

RooAbsPdf* createPolynomialBackground(RooWorkspace* ws, const string& name, RooRealVar& mass, int order) {
    RooArgList coeffs;
    for (int i = 0; i < order; ++i) {
        coeffs.addOwned(*(new RooRealVar(Form("%s_c%d", name.c_str(), i), Form("c%d", i), 0, -1, 1)));
    }

    RooPolynomial bkg((name + "_poly").c_str(), "Polynomial Background", mass, coeffs);
    ws->import(bkg);
    return ws->pdf((name + "_poly").c_str());
}

