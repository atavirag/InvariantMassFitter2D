#include "SignalFunctions.h"
#include "BifurcatedCB.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgList.h"   // used by RooAddPdf(..., RooArgList(g1, g2), RooArgList(frac))
#include <string>

using namespace std;

RooAbsPdf* createGaussianSignal(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar mean((name + "_mean").c_str(), "Mean", 1.865, 1.85, 1.88);
    RooRealVar sigma((name + "_sigma").c_str(), "Sigma", 0.01, 0.001, 0.05);

    RooGaussian gauss((name + "_gauss").c_str(), "Signal Gaussian", mass, mean, sigma);
    ws->import(gauss);
    return ws->pdf((name + "_gauss").c_str());
}

RooAbsPdf* createDoubleGaussianSignal(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar mean((name + "_mean").c_str(), "Mean", 1.865, 1.85, 1.88);
    RooRealVar sigma1((name + "_sigma1").c_str(), "Sigma1", 0.010, 0.001, 0.05);
    RooRealVar sigma2((name + "_sigma2").c_str(), "Sigma2", 0.020, 0.001, 0.05);
    RooRealVar frac((name + "_frac").c_str(), "Fraction", 0.5, 0.0, 1.0);

    RooGaussian g1((name + "_g1").c_str(), "Gaussian 1", mass, mean, sigma1);
    RooGaussian g2((name + "_g2").c_str(), "Gaussian 2", mass, mean, sigma2);

    RooAddPdf doubleGauss((name + "_doubleGauss").c_str(), "Double Gaussian",
                          RooArgList(g1, g2), RooArgList(frac));
    ws->import(doubleGauss);
    return ws->pdf((name + "_doubleGauss").c_str());
}

RooAbsPdf* createCrystalBallSignal(RooWorkspace* ws, const std::string& name, RooRealVar& mass)
{
    // Peak position
    RooRealVar* mu = new RooRealVar(
        Form("mu_%s", name.c_str()), "mu",
        mass.getMin() + 0.5 * (mass.getMax() - mass.getMin()),
        mass.getMin(), mass.getMax());

    // Width
    RooRealVar* sigma = new RooRealVar(
        Form("sigma_%s", name.c_str()), "sigma",
        0.01, 0.001, 0.5);

    // Left tail
    RooRealVar* alphaL = new RooRealVar(
        Form("alphaL_%s", name.c_str()), "alphaL", 1.5, 0.1, 10.0);
    RooRealVar* nL = new RooRealVar(
        Form("nL_%s", name.c_str()), "nL", 1.5, 0.1, 50.0);

    // Right tail
    RooRealVar* alphaR = new RooRealVar(
        Form("alphaR_%s", name.c_str()), "alphaR", 1.5, 0.1, 10.0);
    RooRealVar* nR = new RooRealVar(
        Form("nR_%s", name.c_str()), "nR", 1.5, 0.1, 50.0);

    BifurcatedCB* signalPdf = new BifurcatedCB(
        name.c_str(), "Bifurcated Crystal Ball signal",
        mass, *mu, *sigma, *alphaL, *nL, *alphaR, *nR);

    // Register PDF (and its parameters, by dependency) into the workspace
    ws->import(*signalPdf, RooFit::RecycleConflictNodes());

    return ws->pdf(name.c_str());
}