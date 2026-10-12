#include "ReflectedFunctions.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"

using namespace std;

RooAbsPdf* createGaussianReflected(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar mean((name + "_reflmean").c_str(), "Mean", 1.865, 1.7, 2.30);
    RooRealVar sigma((name + "_reflsigma").c_str(), "Sigma", 0.01, 0.001, 0.05);

    RooGaussian gauss((name + "_reflgauss").c_str(), "Reflected Gaussian", mass, mean, sigma);
    ws->import(gauss);
    return ws->pdf((name + "_reflgauss").c_str());
}

RooAbsPdf* createDoubleGaussianReflected(RooWorkspace* ws, const string& name, RooRealVar& mass) {
    RooRealVar mean((name + "_reflmean").c_str(), "Mean", 1.865, 1.7, 2.30);
    RooRealVar sigma1((name + "_reflsigma1").c_str(), "Sigma1", 0.010, 0.001, 0.2);
    RooRealVar sigma2((name + "_reflsigma2").c_str(), "Sigma2", 0.020, 0.001, 0.2);
    RooRealVar frac((name + "_reflfrac").c_str(), "Fraction", 0.5, 0.0, 1.0);

    RooGaussian g1((name + "_reflg1").c_str(), "Reflected Gaussian 1", mass, mean, sigma1);
    RooGaussian g2((name + "_reflg2").c_str(), "Reflected Gaussian 2", mass, mean, sigma2);

    RooAddPdf doubleGauss((name + "_refldoubleGauss").c_str(), "Reflected Double Gaussian",
                          RooArgList(g1, g2), RooArgList(frac));
    ws->import(doubleGauss);
    return ws->pdf((name + "_refldoubleGauss").c_str());
}
