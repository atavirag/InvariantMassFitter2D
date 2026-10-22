// InvariantMassFitter.cxx
#include "InvariantMassFitter.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooCurve.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../01-prefiltering/CandidatePair.h"
#include "functions/BackgroundFunctions.h"
#include "functions/SignalFunctions.h"
#include "functions/ReflectedFunctions.h"

using namespace RooFit;

namespace {

// Fit helper shared by the prefit and the 2D fit
RooFitResult* doFitTo(RooAbsPdf* pdf, RooDataSet& data, bool weighted) {
    if (weighted) {
        return pdf->fitTo(data, RooFit::Save(),
                          RooFit::Extended(true),
                          RooFit::AsymptoticError(true),
                          RooFit::Strategy(2),
                          RooFit::Optimize(2),
                          RooFit::Offset(true),
                          RooFit::Minimizer("Minuit2", "migrad"),
                          RooFit::PrintLevel(1),
                          RooFit::MaxCalls(100000));
    }
    return pdf->fitTo(data, RooFit::Save(),
                      RooFit::Extended(true),
                      RooFit::PrintLevel(1));
}

// Mean and (effective) sigma of the signal pdf with prefix `pre` (e.g. "m1_sig"),
// whatever the signal function is.
bool getSignalMeanSigma(RooWorkspace* ws, const TString& opt, const std::string& pre,
                        double& mean, double& meanErr, double& sigma, double& sigmaErr) {
    if (opt == "gaussian") {
        RooRealVar* m = ws->var((pre + "_mean").c_str());
        RooRealVar* s = ws->var((pre + "_sigma").c_str());
        if (!m || !s) return false;
        mean = m->getVal(); meanErr = m->getError();
        sigma = s->getVal(); sigmaErr = s->getError();
        return true;
    }
    if (opt == "doubleGaussian") {
        RooRealVar* m  = ws->var((pre + "_mean").c_str());
        RooRealVar* s1 = ws->var((pre + "_sigma1").c_str());
        RooRealVar* s2 = ws->var((pre + "_sigma2").c_str());
        RooRealVar* f  = ws->var((pre + "_frac").c_str());
        if (!m || !s1 || !s2 || !f) return false;
        double fr = f->getVal();
        mean = m->getVal(); meanErr = m->getError();
        sigma = std::sqrt(fr * s1->getVal() * s1->getVal() + (1. - fr) * s2->getVal() * s2->getVal());
        sigmaErr = fr * s1->getError() + (1. - fr) * s2->getError(); // approximate
        return true;
    }
    if (opt == "crystalball") {
        RooRealVar* m = ws->var(Form("mu_%s", pre.c_str()));
        RooRealVar* s = ws->var(Form("sigma_%s", pre.c_str()));
        if (!m || !s) return false;
        mean = m->getVal(); meanErr = m->getError();
        sigma = s->getVal(); sigmaErr = s->getError();
        return true;
    }
    return false;
}

// Suffixes appended by the create*Pdf functions to the pdf names
std::string sigPdfSuffix(const TString& opt) {
    if (opt == "gaussian") return "_gauss";
    if (opt == "doubleGaussian") return "_doubleGauss";
    return ""; // crystalball: pdf name == name
}
std::string bkgPdfSuffix(const TString& opt) {
    if (opt == "chebychev") return "_bkg";
    if (opt == "exponential") return "_exp";
    return "_poly";
}

} // namespace

InvariantMassFitter::InvariantMassFitter(const std::string& osFilename, const std::string& lsFilename)
    : _osTree(nullptr), _lsTree(nullptr), _valid(false),
      _m1Min(1.7), _m1Max(2.1), _m2Min(1.7), _m2Max(2.1),
      _ptMin(1.0), _ptMax(24.0), _reflOverSgn(0.2),
      _efficiencyMap(nullptr), _outputRootFile("")
{
    // Workspace first, so it always exists
    _workspace = std::make_unique<RooWorkspace>("w", "Invariant Mass Fit Workspace");

    // Open input files
    _osFile = std::unique_ptr<TFile>(TFile::Open(osFilename.c_str(), "READ"));
    _lsFile = std::unique_ptr<TFile>(TFile::Open(lsFilename.c_str(), "READ"));

    if (!_osFile || _osFile->IsZombie()) {
        std::cerr << "Error: Cannot open OS file " << osFilename << std::endl;
        return;
    }
    if (!_lsFile || _lsFile->IsZombie()) {
        std::cerr << "Error: Cannot open LS file " << lsFilename << std::endl;
        return;
    }

    // Get trees
    _osTree = dynamic_cast<TTree*>(_osFile->Get("pairs"));
    _lsTree = dynamic_cast<TTree*>(_lsFile->Get("pairs"));

    if (!_osTree) std::cerr << "Error: Cannot find tree 'pairs' in OS file" << std::endl;
    if (!_lsTree) std::cerr << "Error: Cannot find tree 'pairs' in LS file" << std::endl;
    if (!_osTree || !_lsTree) return;

    _valid = true;

    std::cout << "InvariantMassFitter initialized" << std::endl;
    std::cout << "  OS entries: " << _osTree->GetEntries() << std::endl;
    std::cout << "  LS entries: " << _lsTree->GetEntries() << std::endl;
}

InvariantMassFitter::~InvariantMassFitter() {
    // Smart pointers will handle cleanup
}

void InvariantMassFitter::setOutputFile(const std::string& filename) {
    _outputRootFile = filename;
    std::cout << "Output ROOT file set to: " << filename << std::endl;
}

bool InvariantMassFitter::loadEfficiencyMap(const std::string& filename, const std::string& histName) {
    std::unique_ptr<TFile> effFile(TFile::Open(filename.c_str(), "READ"));
    if (!effFile || effFile->IsZombie()) {
        std::cerr << "Error: Cannot open efficiency file " << filename << std::endl;
        return false;
    }

    TH2D* effHist = dynamic_cast<TH2D*>(effFile->Get(histName.c_str()));
    if (!effHist) {
        std::cerr << "Error: Cannot find TH2D " << histName << " in efficiency file" << std::endl;
        return false;
    }

    // Clone the histogram so we own it
    _efficiencyMap = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(effHist->Clone("efficiencyMap")));
    _efficiencyMap->SetDirectory(nullptr); // Detach from file

    std::cout << "Efficiency map loaded: " << histName << " from " << filename << std::endl;
    std::cout << "  Map dimensions: " << _efficiencyMap->GetNbinsX() << " x "
              << _efficiencyMap->GetNbinsY() << " bins" << std::endl;
    return true;
}

// Returns 1/efficiency, or -1 if the candidate cannot be corrected.
double InvariantMassFitter::calculateWeights(float const &pt, float const &y)
{
    if (!_efficiencyMap) return 1.0;

    const double epsilon = 1e-4;

    // If the map is in |y| (Y axis starts at 0), use |y|
    double yv = (_efficiencyMap->GetYaxis()->GetXmin() >= 0.) ? std::fabs(y) : y;

    if (pt < _ptMin || pt > _ptMax) return -1.;

    int ptBin = _efficiencyMap->GetXaxis()->FindBin(pt);
    int yBin  = _efficiencyMap->GetYaxis()->FindBin(yv);
    if (ptBin < 1 || ptBin > _efficiencyMap->GetXaxis()->GetNbins() ||
        yBin  < 1 || yBin  > _efficiencyMap->GetYaxis()->GetNbins()) {
        return -1.;
    }

    double efficiency = _efficiencyMap->GetBinContent(ptBin, yBin);
    return (efficiency > epsilon) ? 1.0 / efficiency : -1.;
}

void InvariantMassFitter::setReflOverSgn(float const& relfOverSgn) {
    _reflOverSgn = relfOverSgn;
}

void InvariantMassFitter::setMassRange(double m1Min, double m1Max, double m2Min, double m2Max) {
    _m1Min = m1Min;
    _m1Max = m1Max;
    _m2Min = m2Min;
    _m2Max = m2Max;

    std::cout << "Mass ranges set: m1=[" << m1Min << ", " << m1Max
              << "], m2=[" << m2Min << ", " << m2Max << "]" << std::endl;
}

void InvariantMassFitter::setPtRange(double ptMin, double ptMax) {
    _ptMin = ptMin;
    _ptMax = ptMax;

    std::cout << "pt range set: [" << ptMin << ", " << ptMax << "]" << std::endl;
}

void InvariantMassFitter::setFittingFunctions(TString const &sgnFunc, TString const &bkgFunc, TString const &reflFunc) {
    _sgnFuncOption = sgnFunc;
    _bkgFuncOption = bkgFunc;
    _reflFuncOption = reflFunc;

    std::cout << "Selected signal function: " << sgnFunc << std::endl;
    std::cout << "Selected background function: " << bkgFunc << std::endl;
    std::cout << "Selected reflected function: " << reflFunc << std::endl;
}

void InvariantMassFitter::createVariables() {
    // Create mass variables
    _m1 = std::make_unique<RooRealVar>("m1", "m_{1} [GeV/c^{2}]", _m1Min, _m1Max);
    _m2 = std::make_unique<RooRealVar>("m2", "m_{2} [GeV/c^{2}]", _m2Min, _m2Max);

    _m1->setBins(50);
    _m2->setBins(50);

    // Kinematic variables (ranges wide enough that setVal never clips)
    _y1 = std::make_unique<RooRealVar>("y1", "y_{1}", -2.0, 2.0);
    _y2 = std::make_unique<RooRealVar>("y2", "y_{2}", -2.0, 2.0);
    _pt1 = std::make_unique<RooRealVar>("pt1", "p_{T 1}", 0.0, 1000.0);
    _pt2 = std::make_unique<RooRealVar>("pt2", "p_{T 2}", 0.0, 1000.0);
    _pairPt = std::make_unique<RooRealVar>("pairPt", "p_{T pair}", 0.0, 1000.0);
    _pairMass = std::make_unique<RooRealVar>("pairMass", "m_{pair}", 0.0, 1000.0);
    _phi1 = std::make_unique<RooRealVar>("phi1", "#phi_{1}", -2*TMath::Pi(), 2*TMath::Pi());
    _phi2 = std::make_unique<RooRealVar>("phi2", "#phi_{2}", -2*TMath::Pi(), 2*TMath::Pi());

    _y1->setBins(50);
    _y2->setBins(50);
    _phi1->setBins(50);
    _phi2->setBins(50);
    _pt1->setBins(100);
    _pt2->setBins(100);
    _pairPt->setBins(100);
    _pairMass->setBins(100);

    // Import into workspace
    _workspace->import(*_m1);
    _workspace->import(*_m2);
    _workspace->import(*_y1);
    _workspace->import(*_y2);
    _workspace->import(*_phi1);
    _workspace->import(*_phi2);
    _workspace->import(*_pt1);
    _workspace->import(*_pt2);
    _workspace->import(*_pairPt);
    _workspace->import(*_pairMass);
}

void InvariantMassFitter::loadData() {
    if (!_valid) {
        std::cerr << "Error: loadData() called on an invalid fitter (input files missing)" << std::endl;
        return;
    }
    if (!_m1 || !_m2) {
        createVariables();
    }

    const bool useEfficiency = (_efficiencyMap != nullptr);

    // Builds one dataset (OS or LS) from a tree
    auto buildDataset = [&](TTree* tree, const char* name, const char* title, bool isOS)
        -> std::unique_ptr<RooDataSet> {

        if (!tree || tree->GetEntries() == 0) return nullptr;

        RooArgSet varSet(*_m1, *_m2, *_y1, *_y2, *_phi1, *_phi2, *_pt1, *_pt2, *_pairPt, *_pairMass);
        RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 1e5);
        RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 1e5);

        std::unique_ptr<RooDataSet> data;
        if (useEfficiency) {
            varSet.add(weightCand1);
            varSet.add(weightCand2);
            data = std::make_unique<RooDataSet>(name, title, varSet, RooFit::WeightVar("combinedWeight"));
        } else {
            data = std::make_unique<RooDataSet>(name, title, varSet);
        }

        float m1_val, m2_val, mbar1_val, mbar2_val;
        float pt1, pt2, y1, y2, phi1, phi2;
        double pairPt, pairMass;
        UChar_t pairType;

        tree->SetBranchAddress("m1", &m1_val);
        tree->SetBranchAddress("m2", &m2_val);
        tree->SetBranchAddress("mbar1", &mbar1_val);
        tree->SetBranchAddress("mbar2", &mbar2_val);
        tree->SetBranchAddress("pairType", &pairType);
        tree->SetBranchAddress("y1", &y1);
        tree->SetBranchAddress("y2", &y2);
        tree->SetBranchAddress("phi1", &phi1);
        tree->SetBranchAddress("phi2", &phi2);
        tree->SetBranchAddress("pt1", &pt1);
        tree->SetBranchAddress("pt2", &pt2);
        tree->SetBranchAddress("pairPt", &pairPt);
        tree->SetBranchAddress("pairMass", &pairMass);

        long nSkippedMass = 0, nSkippedEff = 0;

        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);

            double mass1, mass2;

            if (isOS) {
                if (TESTBIT(pairType, DDbar)) {        // D Dbar
                    mass1 = m1_val;
                    mass2 = mbar2_val;
                } else if (TESTBIT(pairType, DbarD)) { // Dbar D
                    mass1 = mbar1_val;
                    mass2 = m2_val;
                } else {
                    continue; // not OS
                }
            } else {
                if (TESTBIT(pairType, DD)) {           // D D
                    mass1 = m1_val;
                    mass2 = m2_val;
                } else if (TESTBIT(pairType, DbarDbar)) { // Dbar Dbar
                    mass1 = mbar1_val;
                    mass2 = mbar2_val;
                } else {
                    continue; // not LS
                }
            }

            // RooRealVar::setVal clips to the range: skip instead of piling up at the edges
            if (mass1 < _m1Min || mass1 > _m1Max || mass2 < _m2Min || mass2 > _m2Max) {
                nSkippedMass++;
                continue;
            }

            double wCand1 = 1., wCand2 = 1.;
            if (useEfficiency) {
                wCand1 = calculateWeights(pt1, y1); // 1/efficiency or -1
                wCand2 = calculateWeights(pt2, y2);
                if (wCand1 < 0. || wCand2 < 0.) { nSkippedEff++; continue; }
            }

            _m1->setVal(mass1);
            _m2->setVal(mass2);
            _y1->setVal(y1);
            _y2->setVal(y2);
            _phi1->setVal(phi1);
            _phi2->setVal(phi2);
            _pt1->setVal(pt1);
            _pt2->setVal(pt2);
            _pairPt->setVal(pairPt);
            _pairMass->setVal(pairMass);

            if (useEfficiency) {
                varSet.setRealValue("weightCand1", wCand1);
                varSet.setRealValue("weightCand2", wCand2);
                data->add(varSet, wCand1 * wCand2);
            } else {
                data->add(varSet);
            }
        }

        tree->ResetBranchAddresses();

        std::cout << "Loaded " << data->numEntries() << " " << (isOS ? "OS" : "LS") << " events";
        if (useEfficiency) {
            std::cout << " (efficiency-corrected, sum of weights = " << data->sumEntries() << ")";
        }
        std::cout << std::endl;
        if (nSkippedMass > 0)
            std::cout << "  Skipped " << nSkippedMass << " candidates outside the mass window" << std::endl;
        if (nSkippedEff > 0)
            std::cout << "  Skipped " << nSkippedEff << " candidates outside the efficiency map / pt range / eff<1e-4" << std::endl;

        return data;
    };

    _osData = buildDataset(_osTree, "osData", "OS Data", true);
    if (_osData) _workspace->import(*_osData);

    _lsData = buildDataset(_lsTree, "lsData", "LS Data", false);
    if (_lsData) _workspace->import(*_lsData);
}

RooAbsPdf* InvariantMassFitter::createSignalPdf(const std::string& name, RooRealVar& mass) {
    if (_sgnFuncOption == "gaussian")
        return createGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "doubleGaussian")
        return createDoubleGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "crystalball")
        return createCrystalBallSignal(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown signal option: " + std::string(_sgnFuncOption.Data()));
}

RooAbsPdf* InvariantMassFitter::createBackgroundPdf(const std::string& name, RooRealVar& mass) {
    if (_bkgFuncOption == "chebychev")
        return createChebychevBackground(_workspace.get(), name, mass);
    else if (_bkgFuncOption == "exponential")
        return createExponentialBackground(_workspace.get(), name, mass);
    else if (_bkgFuncOption == "polynomial")
        return createPolynomialBackground(_workspace.get(), name, mass, 2);
    else if (_bkgFuncOption == "poly3")
        return createPolynomialBackground(_workspace.get(), name, mass, 3);
    else
        throw std::runtime_error("Unknown background option: " + std::string(_bkgFuncOption.Data()));
}

RooAbsPdf* InvariantMassFitter::createReflectedPdf(const std::string& name, RooRealVar& mass) {
    if (_reflFuncOption == "gaussian")
        return createGaussianReflected(_workspace.get(), name, mass);
    else if (_reflFuncOption == "doubleGaussian")
        return createDoubleGaussianReflected(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown reflected option: " + std::string(_reflFuncOption.Data()));
}

void InvariantMassFitter::buildModel(bool const doReflections) {

    if (!_m1 || !_m2) {
        createVariables();
    }

    // 1D PDFs for each mass dimension
    RooAbsPdf* sig1 = createSignalPdf("m1_sig", *_m1);
    RooAbsPdf* bkg1 = createBackgroundPdf("m1_bkg", *_m1);
    RooAbsPdf* sig2 = createSignalPdf("m2_sig", *_m2);
    RooAbsPdf* bkg2 = createBackgroundPdf("m2_bkg", *_m2);

    // Initial values (sum of weights, so that efficiency-corrected data works)
    double nTot = 1000.;
    if (_lsData)      nTot = _lsData->sumEntries();
    else if (_osData) nTot = _osData->sumEntries();
    std::cout << "Initial total yield: " << nTot << std::endl;

    // Four INDEPENDENT yields: no assumption that m1 and m2 are uncorrelated,
    // and no flat direction in the likelihood.
    RooRealVar *nSgnSgn = new RooRealVar("nSgnSgn", "sgn x sgn yield", 0.05 * nTot, 0.0, 10 * nTot);
    RooRealVar *nSgnBkg = new RooRealVar("nSgnBkg", "sgn x bkg yield", 0.15 * nTot, 0.0, 10 * nTot);
    RooRealVar *nBkgSgn = new RooRealVar("nBkgSgn", "bkg x sgn yield", 0.15 * nTot, 0.0, 10 * nTot);
    RooRealVar *nBkgBkg = new RooRealVar("nBkgBkg", "bkg x bkg yield", 0.65 * nTot, 0.0, 10 * nTot);

    // 2D PDFs as products
    RooProdPdf pdfSigSig("pdfSigSig", "Signal-Signal", RooArgList(*sig1, *sig2));
    RooProdPdf pdfSigBkg("pdfSigBkg", "Signal-Background", RooArgList(*sig1, *bkg2));
    RooProdPdf pdfBkgSig("pdfBkgSig", "Background-Signal", RooArgList(*bkg1, *sig2));
    RooProdPdf pdfBkgBkg("pdfBkgBkg", "Background-Background", RooArgList(*bkg1, *bkg2));

    RooAddPdf* totalPdf = nullptr;

    if (doReflections) {
        RooAbsPdf* refl1 = createReflectedPdf("m1_refl", *_m1);
        RooAbsPdf* refl2 = createReflectedPdf("m2_refl", *_m2);

        RooProdPdf* pdfSigRefl  = new RooProdPdf("pdfSigRefl",  "Signal-Reflected",       RooArgList(*sig1, *refl2));
        RooProdPdf* pdfReflSig  = new RooProdPdf("pdfReflSig",  "Reflected-Signal",       RooArgList(*refl1, *sig2));
        RooProdPdf* pdfReflBkg  = new RooProdPdf("pdfReflBkg",  "Reflected-Background",   RooArgList(*refl1, *bkg2));
        RooProdPdf* pdfBkgRefl  = new RooProdPdf("pdfBkgRefl",  "Background-Reflected",   RooArgList(*bkg1, *refl2));
        RooProdPdf* pdfReflRefl = new RooProdPdf("pdfReflRefl", "Reflected-Reflected",    RooArgList(*refl1, *refl2));

        // Reflection yields tied to the signal-containing yields by a fixed ratio r
        RooAbsReal &r = RooConst(_reflOverSgn);
        RooFormulaVar* nSgnRefl  = new RooFormulaVar("nSgnRefl",  "sgn x refl",  "@0*@1",    RooArgList(r, *nSgnSgn));
        RooFormulaVar* nReflSgn  = new RooFormulaVar("nReflSgn",  "refl x sgn",  "@0*@1",    RooArgList(r, *nSgnSgn));
        RooFormulaVar* nBkgRefl  = new RooFormulaVar("nBkgRefl",  "bkg x refl",  "@0*@1",    RooArgList(r, *nBkgSgn));
        RooFormulaVar* nReflBkg  = new RooFormulaVar("nReflBkg",  "refl x bkg",  "@0*@1",    RooArgList(r, *nSgnBkg));
        RooFormulaVar* nReflRefl = new RooFormulaVar("nReflRefl", "refl x refl", "@0*@0*@1", RooArgList(r, *nSgnSgn));

        totalPdf = new RooAddPdf("model", "background + signal + reflection pdf 2D",
                                 RooArgList(pdfSigSig, *pdfSigRefl, *pdfReflSig, pdfSigBkg,
                                            pdfBkgSig, pdfBkgBkg, *pdfBkgRefl, *pdfReflBkg, *pdfReflRefl),
                                 RooArgList(*nSgnSgn, *nSgnRefl, *nReflSgn, *nSgnBkg,
                                            *nBkgSgn, *nBkgBkg, *nBkgRefl, *nReflBkg, *nReflRefl));
    } else {
        totalPdf = new RooAddPdf("model", "background + signal 2D model",
                                 RooArgList(pdfSigSig, pdfSigBkg, pdfBkgSig, pdfBkgBkg),
                                 RooArgList(*nSgnSgn, *nSgnBkg, *nBkgSgn, *nBkgBkg));
    }

    _workspace->import(*totalPdf);

    std::cout << "2D model built" << std::endl;
}

RooFitResult* InvariantMassFitter::performPrefit1D(RooDataSet* data, RooRealVar& mass,
                                                    const std::string& name, bool useReflections) {
    std::cout << "\n=== Performing 1D prefit for " << name << " ===" << std::endl;

    RooAbsPdf* sig = createSignalPdf(name + "_sig", mass);
    RooAbsPdf* bkg = createBackgroundPdf(name + "_bkg", mass);

    // Sum of weights: correct for weighted (efficiency-corrected) data, == numEntries() otherwise
    double nEntries = data->sumEntries();

    RooRealVar* nSig = new RooRealVar((name + "_nSig").c_str(), "Signal yield",
                                      0.5 * nEntries, 0.0, 10 * nEntries);
    RooRealVar* nBkg = new RooRealVar((name + "_nBkg").c_str(), "Background yield",
                                      0.2 * nEntries, 0.0, 10 * nEntries);

    RooAddPdf* model = nullptr;

    if (useReflections) {
        RooAbsPdf* refl = createReflectedPdf(name + "_refl", mass);
        RooRealVar* reflOverSgn = new RooRealVar((name + "_reflOverSgn").c_str(),
                                                 "Reflection/Signal ratio",
                                                 _reflOverSgn, 0.0, 1.0);
        reflOverSgn->setConstant(true); // Keep ratio fixed in prefit

        RooFormulaVar* nRefl = new RooFormulaVar((name + "_nRefl").c_str(),
                                                 "Reflection yield", "@0 * @1",
                                                 RooArgList(*reflOverSgn, *nSig));

        model = new RooAddPdf((name + "_model").c_str(), "Signal + Reflection + Background",
                              RooArgList(*sig, *refl, *bkg),
                              RooArgList(*nSig, *nRefl, *nBkg));
    } else {
        model = new RooAddPdf((name + "_model").c_str(), "Signal + Background",
                              RooArgList(*sig, *bkg),
                              RooArgList(*nSig, *nBkg));
    }

    RooFitResult* result = doFitTo(model, *data, _efficiencyMap != nullptr);
    if (!result) {
        std::cerr << "Error: prefit " << name << " returned no result" << std::endl;
        return nullptr;
    }
    result->Print("v");

    std::cout << "Prefit results for " << name << ":" << std::endl;
    std::cout << "  Signal yield: " << nSig->getVal() << " +- " << nSig->getError() << std::endl;
    std::cout << "  Background yield: " << nBkg->getVal() << " +- " << nBkg->getError() << std::endl;

    _workspace->import(*model);
    _workspace->import(*result, (name + "_prefitResult").c_str());

    return result;
}

void InvariantMassFitter::performPrefit(const std::string& dataType, bool useReflections) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "=== Starting 1D Prefit Procedure ===" << std::endl;
    std::cout << "========================================" << std::endl;

    if (!_osData && !_lsData) {
        std::cerr << "Error: No data loaded. Call loadData() first." << std::endl;
        return;
    }

    if (!_m1 || !_m2) {
        createVariables();
    }

    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n--- OS Data Prefits ---" << std::endl;
        _prefitResult_os_m1.reset(performPrefit1D(_osData.get(), *_m1, "prefit_os_m1", useReflections));
        _prefitResult_os_m2.reset(performPrefit1D(_osData.get(), *_m2, "prefit_os_m2", useReflections));
    }

    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n--- LS Data Prefits ---" << std::endl;
        _prefitResult_ls_m1.reset(performPrefit1D(_lsData.get(), *_m1, "prefit_ls_m1", useReflections));
        _prefitResult_ls_m2.reset(performPrefit1D(_lsData.get(), *_m2, "prefit_ls_m2", useReflections));
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "=== 1D Prefit Procedure Complete ===" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void InvariantMassFitter::extractPrefitParameters(const std::string& prefitName,
                                                   const std::string& targetName) {
    // Transfers all shape parameters of the prefit pdf `prefitName` to the 2D-model pdf `targetName`.
    // The target name is obtained by replacing the prefit prefix, so it works for any
    // signal / background / reflection function (including "mu_<name>" style names).

    std::cout << "  Extracting prefit parameters from " << prefitName << " to " << targetName << std::endl;

    const bool isRefl = (prefitName.find("_refl") != std::string::npos);

    RooArgSet allVars = _workspace->allVars();
    int nParamsSet = 0;
    for (RooAbsArg* arg : allVars) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
        if (!var) continue;
        std::string varName = var->GetName();

        size_t pos = varName.find(prefitName);
        if (pos == std::string::npos) continue;

        // Skip yield / ratio parameters
        if (varName.find("_nBkg") != std::string::npos ||
            varName.find("_nSig") != std::string::npos ||
            varName.find("_nRefl") != std::string::npos ||
            varName.find("_reflOverSgn") != std::string::npos) {
            continue;
        }

        std::string targetParamName = varName;
        targetParamName.replace(pos, prefitName.length(), targetName);

        RooRealVar* targetParam = _workspace->var(targetParamName.c_str());
        if (!targetParam) continue;

        targetParam->setVal(var->getVal());
        if (isRefl) {
            targetParam->setConstant(true);
            std::cout << "    " << targetParamName << " set CONSTANT to " << var->getVal() << std::endl;
        } else {
            std::cout << "    " << targetParamName << " set to " << var->getVal() << std::endl;
        }
        nParamsSet++;
    }

    if (nParamsSet == 0) {
        std::cout << "    Warning: No parameters found to transfer" << std::endl;
    } else {
        std::cout << "    Successfully set " << nParamsSet << " parameter(s)" << std::endl;
    }
}

void InvariantMassFitter::setPrefitParameters(bool const doReflections) {
    std::cout << "\n=== Setting 2D fit initial parameters from prefit ===" << std::endl;

    auto getParam = [this](const std::string& name) -> RooRealVar* {
        return _workspace->var(name.c_str());
    };

    // Shape parameters (taken from the OS prefit)
    if (_prefitResult_os_m1) {
        std::cout << "\n--- Setting m1 shape parameters from OS prefit ---" << std::endl;
        extractPrefitParameters("prefit_os_m1_sig", "m1_sig");
        extractPrefitParameters("prefit_os_m1_bkg", "m1_bkg");
        if (doReflections) extractPrefitParameters("prefit_os_m1_refl", "m1_refl");
    }
    if (_prefitResult_os_m2) {
        std::cout << "\n--- Setting m2 shape parameters from OS prefit ---" << std::endl;
        extractPrefitParameters("prefit_os_m2_sig", "m2_sig");
        extractPrefitParameters("prefit_os_m2_bkg", "m2_bkg");
        if (doReflections) extractPrefitParameters("prefit_os_m2_refl", "m2_refl");
    }

    // Yield start values from the two 1D projections: total N and signal fractions f1, f2
    RooRealVar *s1 = getParam("prefit_os_m1_nSig"), *b1 = getParam("prefit_os_m1_nBkg");
    RooRealVar *s2 = getParam("prefit_os_m2_nSig"), *b2 = getParam("prefit_os_m2_nBkg");
    RooRealVar *nSS = getParam("nSgnSgn"), *nSB = getParam("nSgnBkg");
    RooRealVar *nBS = getParam("nBkgSgn"), *nBB = getParam("nBkgBkg");

    if (s1 && b1 && s2 && b2 && nSS && nSB && nBS && nBB) {
        double t1 = s1->getVal() + b1->getVal();
        double t2 = s2->getVal() + b2->getVal();
        if (t1 > 0. && t2 > 0.) {
            double N  = 0.5 * (t1 + t2);
            double f1 = s1->getVal() / t1;
            double f2 = s2->getVal() / t2;
            nSS->setVal(N * f1 * f2);
            nSB->setVal(N * f1 * (1. - f2));
            nBS->setVal(N * (1. - f1) * f2);
            nBB->setVal(N * (1. - f1) * (1. - f2));
            std::cout << "  Yield starts: SS=" << nSS->getVal() << " SB=" << nSB->getVal()
                      << " BS=" << nBS->getVal() << " BB=" << nBB->getVal() << std::endl;
        }
    }

    std::cout << "\n=== Initial parameters set from prefit ===" << std::endl;
}

void InvariantMassFitter::plotPrefitResults(const std::string& outputPrefix) {
    std::cout << "\n=== Plotting prefit results ===" << std::endl;

    TCanvas* c = new TCanvas("c_prefit", "Prefit Results", 1600, 1200);
    c->Divide(2, 2);

    const std::string sigSuf = sigPdfSuffix(_sgnFuncOption);
    const std::string bkgSuf = bkgPdfSuffix(_bkgFuncOption);

    struct Panel { int pad; bool ok; RooDataSet* data; RooRealVar* var; std::string tag; std::string title; };
    std::vector<Panel> panels = {
        {1, (bool)_prefitResult_os_m1, _osData.get(), _m1.get(), "prefit_os_m1", "OS m_{1} Prefit"},
        {2, (bool)_prefitResult_os_m2, _osData.get(), _m2.get(), "prefit_os_m2", "OS m_{2} Prefit"},
        {3, (bool)_prefitResult_ls_m1, _lsData.get(), _m1.get(), "prefit_ls_m1", "LS m_{1} Prefit"},
        {4, (bool)_prefitResult_ls_m2, _lsData.get(), _m2.get(), "prefit_ls_m2", "LS m_{2} Prefit"},
    };

    for (auto& p : panels) {
        if (!p.ok || !p.data) continue;
        c->cd(p.pad);
        RooPlot* frame = p.var->frame(RooFit::Title(p.title.c_str()));
        p.data->plotOn(frame, RooFit::DrawOption("pz"));

        RooAbsPdf* model = _workspace->pdf((p.tag + "_model").c_str());
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components((p.tag + "_sig" + sigSuf).c_str()),
                          DrawOption("F"), LineWidth(2), LineStyle(1), LineColor(kGreen + 2), FillColor(kGreen - 2));
            model->plotOn(frame, RooFit::Components((p.tag + "_bkg" + bkgSuf).c_str()),
                          RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        }
        frame->Draw();
    }

    std::string filename = outputPrefix + "_prefit.pdf";
    c->SaveAs(filename.c_str());
    c->SaveAs((outputPrefix + "_prefit.png").c_str());

    std::cout << "Prefit plots saved to " << filename << std::endl;

    delete c;
}

void InvariantMassFitter::fit(const std::string& dataType) {
    if (!_workspace->pdf("model")) {
        std::cerr << "Error: Model not built. Call buildModel() first." << std::endl;
        return;
    }

    if (!_osData && !_lsData) {
        loadData();
    }
    if (!_osData && !_lsData) return;

    RooAbsPdf* pdf = _workspace->pdf("model");
    const bool isWeighted = (_efficiencyMap != nullptr);

    // Remember the start values (from the prefit) so that the LS fit does not start from the OS best fit
    {
        RooArgSet* p = pdf->getParameters(_osData ? *_osData : *_lsData);
        _workspace->saveSnapshot("initParams", *p);
        delete p;
    }

    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n=== Fitting OS data ===" << std::endl;
        _workspace->loadSnapshot("initParams");
        RooFitResult* result = doFitTo(pdf, *_osData, isWeighted);
        if (result) {
            result->Print("v");
            _workspace->import(*result, "osFitResult");
            RooArgSet* params = pdf->getParameters(*_osData);
            _workspace->saveSnapshot("osParams", *params);
            delete params;
            delete result;
        }
    }

    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n=== Fitting LS data ===" << std::endl;
        _workspace->loadSnapshot("initParams");
        RooFitResult* result = doFitTo(pdf, *_lsData, isWeighted);
        if (result) {
            result->Print("v");
            _workspace->import(*result, "lsFitResult");
            RooArgSet* params = pdf->getParameters(*_lsData);
            _workspace->saveSnapshot("lsParams", *params);
            delete params;
            delete result;
        }
    }
}

// Sets the 3-sigma window from the CURRENTLY LOADED snapshot, on both the local
// and the workspace m1/m2 variables (the pdfs depend on the latter).
bool InvariantMassFitter::setSigmaWindow() {
    double mean1, mean1Err, sig1, sig1Err, mean2, mean2Err, sig2, sig2Err;
    if (!getSignalMeanSigma(_workspace.get(), _sgnFuncOption, "m1_sig", mean1, mean1Err, sig1, sig1Err) ||
        !getSignalMeanSigma(_workspace.get(), _sgnFuncOption, "m2_sig", mean2, mean2Err, sig2, sig2Err)) {
        std::cerr << "Error: Cannot find signal parameters in workspace" << std::endl;
        return false;
    }
    double lo1 = std::max(mean1 - 3 * sig1, _m1Min), hi1 = std::min(mean1 + 3 * sig1, _m1Max);
    double lo2 = std::max(mean2 - 3 * sig2, _m2Min), hi2 = std::min(mean2 + 3 * sig2, _m2Max);
    for (RooRealVar* v : {_m1.get(), _workspace->var("m1")}) if (v) v->setRange("3sigmaRange", lo1, hi1);
    for (RooRealVar* v : {_m2.get(), _workspace->var("m2")}) if (v) v->setRange("3sigmaRange", lo2, hi2);
    return true;
}

// Loads `snapshot`, sets the 3-sigma window and returns the yield of every 2D component in it.
// Returns an empty vector if the window cannot be set.
auto InvariantMassFitter::computeYields3Sigma(const std::string& snapshot, RooFitResult* res, bool doReflections)
    -> std::vector<YieldEntry> {
    std::vector<YieldEntry> out;
    if (!res) return out;

    _workspace->loadSnapshot(snapshot.c_str());
    if (!setSigmaWindow()) return out;

    std::vector<std::string> comps = {"SgnSgn", "SgnBkg", "BkgSgn", "BkgBkg"};
    if (doReflections) {
        for (const char* s : {"SgnRefl", "ReflSgn", "BkgRefl", "ReflBkg", "ReflRefl"}) comps.push_back(s);
    }

    for (const std::string& comp : comps) {
        // pdfs are named pdfSigSig, pdfSigBkg, ... ("Sig"); yields nSgnSgn, nSgnBkg, ... ("Sgn")
        std::string pdfName = "pdf" + comp;
        std::string nName   = "n"   + comp;
        for (size_t p; (p = pdfName.find("Sgn")) != std::string::npos; ) pdfName.replace(p, 3, "Sig");
        for (size_t p; (p = nName.find("Sig"))   != std::string::npos; ) nName.replace(p, 3, "Sgn");

        YieldEntry e{comp, 0., 0., false};
        RooAbsPdf*  f = _workspace->pdf(pdfName.c_str());
        RooAbsReal* n = dynamic_cast<RooAbsReal*>(_workspace->arg(nName.c_str()));
        if (f && n) {
            RooRealVar* y = getYieldInRange(res, f, n, "3sigmaRange");
            e.value = y->getVal();
            e.error = y->getError();
            e.ok = true;
            delete y;
        }
        out.push_back(e);
    }
    return out;
}

// Text lines with the 3-sigma yields for the plots. Empty if there is no fit result.
std::vector<std::string> InvariantMassFitter::yieldTextLines(const std::string& label,
                                                             const std::string& snapshot,
                                                             const std::string& fitResultName,
                                                             bool doReflections) {
    std::vector<std::string> lines;
    RooFitResult* res = dynamic_cast<RooFitResult*>(_workspace->obj(fitResultName.c_str()));
    if (!res) return lines;

    std::vector<YieldEntry> ys = computeYields3Sigma(snapshot, res, doReflections);
    if (ys.empty()) return lines;

    auto pretty = [](const std::string& s) {
        size_t split = (s.compare(0, 4, "Refl") == 0) ? 4 : 3;     // first token length
        auto tok = [](const std::string& t) { return t == "Sgn" ? std::string("Sgn") : t; };
        return tok(s.substr(0, split)) + " #times " + tok(s.substr(split));
    };

    lines.push_back(label + " yields in 3#sigma range:");
    for (const YieldEntry& e : ys) {
        if (!e.ok) continue;
        lines.push_back(std::string(Form("%s: %.3g #pm %.2g", pretty(e.label).c_str(), e.value, e.error)));
    }
    return lines;
}

void InvariantMassFitter::plotModelOnFrame(RooPlot* frame, RooAbsPdf* pdf, bool const doReflections,
                                           const std::vector<std::string>& textLines, bool showLegend) {
    pdf->plotOn(frame, Name("model_curve"));

    pdf->plotOn(frame, Components("pdfBkgBkg"), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfSigBkg"), LineStyle(3), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfBkgSig"), LineStyle(4), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
    pdf->plotOn(frame, DrawOption("F"), Components("pdfSigSig"), LineWidth(2), LineStyle(1), LineColor(kGreen + 2), FillColor(kGreen - 2), Name("sgnSgnFunc2D_curve"));

    if (doReflections) {
        pdf->plotOn(frame, Components("pdfBkgRefl"), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfSigRefl"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflBkg"), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflSig"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflRefl"), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
    }

    // Legend naming the components, attached to the frame so it is drawn and saved with it
    TLegend* leg = new TLegend(0.55, doReflections ? 0.45 : 0.62, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    auto addEntry = [&](const char* objName, const char* text, const char* opt) {
        if (TObject* o = frame->findObject(objName)) leg->AddEntry(o, text, opt);
    };

    for (const char* dn : {"data_cand1", "data_os_m1", "data_os_m2", "data_ls_m1", "data_ls_m2"}) {
        if (frame->findObject(dn)) { addEntry(dn, "Data", "pe"); break; }
    }
    addEntry("model_curve",          "Total fit",                          "l");
    addEntry("sgnSgnFunc2D_curve",   "Signal #times Signal",               "f");
    addEntry("sgnBkgFunc2D_curve",   "Signal #times Background",           "l");
    addEntry("bkgSgnFunc2D_curve",   "Background #times Signal",           "l");
    addEntry("bkgBkgFunc2D_curve",   "Background #times Background",       "l");
    if (doReflections) {
        addEntry("sgnReflFunc2D_curve",  "Signal #times Reflection",       "l");
        addEntry("reflSgnFunc2D_curve",  "Reflection #times Signal",       "l");
        addEntry("bkgReflFunc2D_curve",  "Background #times Reflection",   "l");
        addEntry("reflBkgFunc2D_curve",  "Reflection #times Background",   "l");
        addEntry("reflReflFunc2D_curve", "Reflection #times Reflection",   "l");
    }
    if (showLegend) frame->addObject(leg);
    else            delete leg;

    // Text box with the 3-sigma yields (top-right; below the legend if it is shown)
    if (!textLines.empty()) {
        const double legBottom = doReflections ? 0.45 : 0.62;
        const double top       = showLegend ? legBottom - 0.01 : 0.89;
        const double height    = 0.042 * textLines.size();
        TPaveText* pave = new TPaveText(0.55, std::max(0.05, top - height), 0.89, top, "NDC");
        pave->SetBorderSize(0);
        pave->SetFillStyle(0);
        pave->SetTextAlign(12);
        pave->SetTextSize(0.037);
        for (const std::string& l : textLines) pave->AddText(l.c_str());
        frame->addObject(pave);
    }

    frame->Draw();
}

double InvariantMassFitter::computeChi2(RooPlot* frame, RooDataSet* data, int nFitParams) {
    // frame->chiSquare() returns chi2/ndf with ndf = nBins - nFitParams.
    // For weighted data RooFit uses sum-of-weights-squared errors on the data points.
    const bool isWeighted = (_efficiencyMap != nullptr);

    double chi2Ndf = frame->chiSquare("model_curve", "data_cand1", nFitParams);

    int nBins = frame->GetXaxis()->GetNbins();
    int ndf   = nBins - nFitParams;
    double chi2 = chi2Ndf * ndf;

    std::cout << "  Chi2 / ndf  = " << chi2 << " / " << ndf
              << " = " << chi2Ndf << std::endl;

    if (isWeighted) {
        std::cout << "  (weighted fit: N_eff = " << data->numEntries()
                  << ", sum of weights = " << data->sumEntries() << ")" << std::endl;
    }

    return chi2Ndf;
}

void InvariantMassFitter::plotFit(const std::string& outputPrefix, bool const doReflections) {
    if (!_workspace->pdf("model")) {
        std::cerr << "Error: No model to plot" << std::endl;
        return;
    }
    if (!_osData && !_lsData) {
        std::cerr << "Error: No data to plot" << std::endl;
        return;
    }

    RooAbsPdf* pdf = _workspace->pdf("model");

    TCanvas* c = new TCanvas("c_fit", "Fit Results", 1200, 800);
    c->Divide(2, 2);

    const double binWidth = 0.005; // Desired bin width
    int nBins = TMath::Nint((_m1Max - _m1Min) / binWidth);
    _m1->setBins(nBins);
    _m2->setBins(nBins);

    // Count only floating (non-constant) parameters for chi2/ndf
    int nFitParams = 0;
    {
        RooArgSet* allParams = pdf->getParameters(_osData ? *_osData : *_lsData);
        for (RooAbsArg* arg : *allParams) {
            RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
            if (var && !var->isConstant()) nFitParams++;
        }
        delete allParams;
    }
    std::cout << "  Floating fit parameters: " << nFitParams << std::endl;

    if (_osData) {
        _workspace->loadSnapshot("osParams");
        const std::vector<std::string> osText = yieldTextLines("OS", "osParams", "osFitResult", doReflections);

        c->cd(1);
        RooPlot* frame1 = _m1->frame(RooFit::Title("OS: m_{1} projection"));
        _osData->plotOn(frame1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame1, pdf, doReflections, {}, true);          // legend only
        std::cout << "OS m1: "; computeChi2(frame1, _osData.get(), nFitParams);

        c->cd(2);
        RooPlot* frame2 = _m2->frame(RooFit::Title("OS: m_{2} projection"));
        _osData->plotOn(frame2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame2, pdf, doReflections, osText, false);     // numbers only
        std::cout << "OS m2: "; computeChi2(frame2, _osData.get(), nFitParams);
    }

    if (_lsData) {
        _workspace->loadSnapshot("lsParams");
        const std::vector<std::string> lsText = yieldTextLines("LS", "lsParams", "lsFitResult", doReflections);

        c->cd(3);
        RooPlot* frame3 = _m1->frame(RooFit::Title("LS: m_{1} projection"));
        _lsData->plotOn(frame3, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame3, pdf, doReflections, {}, true);          // legend only
        std::cout << "LS m1: "; computeChi2(frame3, _lsData.get(), nFitParams);

        c->cd(4);
        RooPlot* frame4 = _m2->frame(RooFit::Title("LS: m_{2} projection"));
        _lsData->plotOn(frame4, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame4, pdf, doReflections, lsText, false);     // numbers only
        std::cout << "LS m2: "; computeChi2(frame4, _lsData.get(), nFitParams);
    }

    c->SaveAs((outputPrefix + ".pdf").c_str());
    c->SaveAs((outputPrefix + ".png").c_str());

    std::cout << "Plots saved to " << outputPrefix << ".pdf/png" << std::endl;

    delete c;
}

void InvariantMassFitter::saveResultsToFile(bool const doReflections) {
    if (_outputRootFile.empty()) {
        std::cerr << "Error: Output file not set. Call setOutputFile() first." << std::endl;
        return;
    }

    TFile* outFile = TFile::Open(_outputRootFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << _outputRootFile << std::endl;
        return;
    }

    std::cout << "Saving results to " << _outputRootFile << std::endl;

    // Save the workspace (contains all PDFs, parameters and snapshots)
    _workspace->Write();

    RooAbsPdf* pdf = _workspace->pdf("model");

    // ---- 1D projections (data histogram + RooPlot frame with the fit curves) ----
    auto writeProjections = [&](RooDataSet* data, const std::string& tag, const std::string& snapshot,
                                const std::string& label) {
        if (!data || !pdf) return;
        _workspace->loadSnapshot(snapshot.c_str());
        const std::vector<std::string> textLines =
            yieldTextLines(label, snapshot, tag + "FitResult", doReflections); // also loads the snapshot

        struct Proj { RooRealVar* var; std::string idx; };
        for (const Proj& pr : {Proj{_m1.get(), "m1"}, Proj{_m2.get(), "m2"}}) {
            RooPlot* frame = pr.var->frame(RooFit::Title((label + ": " + pr.idx + " projection").c_str()));
            data->plotOn(frame, RooFit::Name(("data_" + tag + "_" + pr.idx).c_str()), RooFit::DrawOption("pz"));
            // m1 projection: legend only; m2 projection: 3-sigma yields only
            if (pr.idx == "m1") plotModelOnFrame(frame, pdf, doReflections, {}, true);
            else                plotModelOnFrame(frame, pdf, doReflections, textLines, false);

            TH1* h = data->createHistogram(("h_data_" + tag + "_" + pr.idx).c_str(), *pr.var);
            h->SetTitle((label + " Data " + pr.idx + " projection;" + pr.idx + " [GeV/c^{2}];Entries").c_str());
            h->Write();

            frame->Write(("frame_" + tag + "_" + pr.idx).c_str());
        }
    };

    writeProjections(_osData.get(), "os", "osParams", "OS");
    writeProjections(_lsData.get(), "ls", "lsParams", "LS");

    // ---- 2D data histograms ----
    if (_osData) {
        TH2* h2_os = _osData->createHistogram(*_m1, *_m2);
        h2_os->SetName("h2_data_os");
        h2_os->SetTitle("OS Data 2D;m_{1} [GeV/c^{2}];m_{2} [GeV/c^{2}]");
        h2_os->Write();
    }
    if (_lsData) {
        TH2* h2_ls = _lsData->createHistogram(*_m1, *_m2);
        h2_ls->SetName("h2_data_ls");
        h2_ls->SetTitle("LS Data 2D;m_{1} [GeV/c^{2}];m_{2} [GeV/c^{2}]");
        h2_ls->Write();
    }

    // ---- Fit results ----
    RooFitResult* osRes = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
    RooFitResult* lsRes = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));
    if (osRes) osRes->Write("osFitResult");
    if (lsRes) lsRes->Write("lsFitResult");

    // ---- Yields in the 3-sigma window (computeYields3Sigma is shared with plotFit) ----
    auto writeYields = [&](RooFitResult* res, const std::string& snapshot, const std::string& tag,
                           const std::string& label) {
        if (!res) {
            std::cerr << "[" << label << " yields] skipped: no fit result in workspace" << std::endl;
            return;
        }
        std::vector<YieldEntry> ys = computeYields3Sigma(snapshot, res, doReflections);
        if (ys.empty()) {
            std::cerr << "[" << label << " yields] skipped: 3-sigma window could not be set" << std::endl;
            return;
        }

        std::cout << "\n=== " << label << " Yields in 3-sigma range ===" << std::endl;

        TH1D* h = new TH1D(("hYields_" + tag).c_str(), (label + " Yields in 3#sigma range").c_str(),
                           ys.size(), 0, ys.size());
        for (size_t i = 0; i < ys.size(); i++) {
            h->GetXaxis()->SetBinLabel(i + 1, ys[i].label.c_str());
            if (!ys[i].ok) {
                std::cerr << "Warning: component " << ys[i].label << " not found in workspace" << std::endl;
                continue;
            }
            std::cout << " " << ys[i].label << " yield: " << ys[i].value << " +- " << ys[i].error << std::endl;
            h->SetBinContent(i + 1, ys[i].value);
            h->SetBinError(i + 1, ys[i].error);
        }
        h->Write();
    };

    writeYields(osRes, "osParams", "os", "OS");
    writeYields(lsRes, "lsParams", "ls", "LS");

    std::cout << "\n ---------------------------------------------- " << std::endl;

    // ---- Fit parameters (mean, sigma, chi2) for OS and LS as single-bin TH1D ----
    auto makeFitParamsHisto = [&](const std::string& name, RooDataSet* data,
                                  const std::string& snapshotName) -> bool {

        _workspace->loadSnapshot(snapshotName.c_str());

        double mean1, mean1Err, sig1, sig1Err, mean2, mean2Err, sig2, sig2Err;
        if (!getSignalMeanSigma(_workspace.get(), _sgnFuncOption, "m1_sig", mean1, mean1Err, sig1, sig1Err) ||
            !getSignalMeanSigma(_workspace.get(), _sgnFuncOption, "m2_sig", mean2, mean2Err, sig2, sig2Err)) {
            std::cerr << "Warning: signal parameters not found for " << name << std::endl;
            return false;
        }

        int nBinsFrame = TMath::Nint((_m1Max - _m1Min) / 0.005);
        _m1->setBins(nBinsFrame);
        _m2->setBins(nBinsFrame);

        RooAbsPdf* pdfModel = _workspace->pdf("model");

        int nFloating = 0;
        {
            RooArgSet* allParams = pdfModel->getParameters(*data);
            for (RooAbsArg* arg : *allParams) {
                RooRealVar* v = dynamic_cast<RooRealVar*>(arg);
                if (v && !v->isConstant()) nFloating++;
            }
            delete allParams;
        }

        RooPlot* f1 = _m1->frame();
        data->plotOn(f1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f1, RooFit::Name("model_curve"));
        double chi2_m1 = f1->chiSquare("model_curve", "data_cand1", nFloating);

        RooPlot* f2 = _m2->frame();
        data->plotOn(f2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f2, RooFit::Name("model_curve"));
        double chi2_m2 = f2->chiSquare("model_curve", "data_cand1", nFloating);

        auto makeHisto = [](const std::string& hname, const std::string& title,
                            double val, double err) -> TH1D* {
            TH1D* h = new TH1D(hname.c_str(), title.c_str(), 1, 0, 1);
            h->SetBinContent(1, val);
            h->SetBinError(1, err);
            return h;
        };

        std::vector<TH1D*> histos = {
            makeHisto(name + "_m1Mean",     "m1 mean [GeV/c^{2}]",   mean1, mean1Err),
            makeHisto(name + "_m1Sigma",    "m1 sigma [GeV/c^{2}]",  sig1,  sig1Err),
            makeHisto(name + "_m2Mean",     "m2 mean [GeV/c^{2}]",   mean2, mean2Err),
            makeHisto(name + "_m2Sigma",    "m2 sigma [GeV/c^{2}]",  sig2,  sig2Err),
            makeHisto(name + "_chi2ndf_m1", "chi2/ndf m1 projection", chi2_m1, 0.0),
            makeHisto(name + "_chi2ndf_m2", "chi2/ndf m2 projection", chi2_m2, 0.0),
        };
        for (TH1D* h : histos) { h->Write(); delete h; }

        delete f1;
        delete f2;
        return true;
    };

    if (_osData && osRes) makeFitParamsHisto("hOS", _osData.get(), "osParams");
    if (_lsData && lsRes) makeFitParamsHisto("hLS", _lsData.get(), "lsParams");

    if (_efficiencyMap) {
        _efficiencyMap->Write("efficiencyMap");
    }

    outFile->Close();
    delete outFile;

    std::cout << "Results saved successfully to " << _outputRootFile << std::endl;
    std::cout << "File contains:" << std::endl;
    std::cout << "  - Workspace 'w' with all PDFs, parameters and snapshots" << std::endl;
    std::cout << "  - 1D data histograms: h_data_os_m1, h_data_os_m2, h_data_ls_m1, h_data_ls_m2" << std::endl;
    std::cout << "  - 2D data histograms: h2_data_os, h2_data_ls" << std::endl;
    std::cout << "  - RooPlot frames (data + fit curves): frame_os_m1, frame_os_m2, frame_ls_m1, frame_ls_m2" << std::endl;
    std::cout << "  - Fit results: osFitResult, lsFitResult" << std::endl;
    std::cout << "  - Yield histograms (3-sigma range): hYields_os, hYields_ls" << std::endl;
}

void InvariantMassFitter::saveWorkspace(const std::string& filename) {
    _workspace->writeToFile(filename.c_str());
    std::cout << "Workspace saved to " << filename << std::endl;
}

void InvariantMassFitter::printFitResults() {
    std::cout << "\n=== Fit Results Summary ===" << std::endl;

    RooFitResult* osRes = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
    RooFitResult* lsRes = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));

    if (osRes) {
        std::cout << "\n--- OS Fit ---" << std::endl;
        osRes->Print("v");
    }
    if (lsRes) {
        std::cout << "\n--- LS Fit ---" << std::endl;
        lsRes->Print("v");
    }

    auto printYields = [&](const std::string& snapshot, const char* label) {
        _workspace->loadSnapshot(snapshot.c_str());
        RooAbsReal* ss = dynamic_cast<RooAbsReal*>(_workspace->arg("nSgnSgn"));
        RooAbsReal* sb = dynamic_cast<RooAbsReal*>(_workspace->arg("nSgnBkg"));
        RooAbsReal* bs = dynamic_cast<RooAbsReal*>(_workspace->arg("nBkgSgn"));
        RooAbsReal* bb = dynamic_cast<RooAbsReal*>(_workspace->arg("nBkgBkg"));
        if (!ss || !sb || !bs || !bb) return;

        auto err = [](RooAbsReal* v) {
            RooRealVar* rv = dynamic_cast<RooRealVar*>(v);
            return rv ? rv->getError() : 0.;
        };

        std::cout << "\n=== " << label << " Yield Summary (full range) ===" << std::endl;
        std::cout << "Signal-Signal:         " << ss->getVal() << " +- " << err(ss) << std::endl;
        std::cout << "Signal-Background:     " << sb->getVal() << " +- " << err(sb) << std::endl;
        std::cout << "Background-Signal:     " << bs->getVal() << " +- " << err(bs) << std::endl;
        std::cout << "Background-Background: " << bb->getVal() << " +- " << err(bb) << std::endl;
        if (_efficiencyMap) std::cout << "(Yields are efficiency-corrected)" << std::endl;
    };

    if (osRes) printYields("osParams", "OS");
    if (lsRes) printYields("lsParams", "LS");
}

RooRealVar *InvariantMassFitter::getYieldInRange(RooFitResult *fitResult, RooAbsPdf* function,
                                                 RooAbsReal* nCands, TString range)
{
    // Integrate over the WORKSPACE variables, which carry the named range
    RooArgSet obs(*_workspace->var("m1"), *_workspace->var("m2"));

    RooAbsReal *intRange = function->createIntegral(obs, Range(range));
    RooAbsReal *intTotal = function->createIntegral(obs);

    double frac = intRange->getVal() / intTotal->getVal();
    double yieldRange = nCands->getVal() * frac;

    // Error propagation on the yield only (shape-parameter uncertainty on `frac` is neglected)
    double nCandsErr = nCands->getPropagatedError(*fitResult);
    double errorRange = nCandsErr * frac;

    delete intRange;
    delete intTotal;

    RooRealVar *yield = new RooRealVar("yield", "Yield in 3-sigma range", yieldRange);
    yield->setError(errorRange);
    return yield;
}
