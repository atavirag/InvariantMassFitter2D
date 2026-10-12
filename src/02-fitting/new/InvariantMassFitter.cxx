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
#include "RooPlot.h"
#include "RooArgSet.h"
#include "TMath.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>

#include "functions/BackgroundFunctions.h"
#include "functions/SignalFunctions.h"
#include "functions/ReflectedFunctions.h"

using namespace RooFit;

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

InvariantMassFitter::InvariantMassFitter()
    : _m1Min(1.7), _m1Max(2.1), _m2Min(1.7), _m2Max(2.1),
      _reflOverSgn(0.2), _datasetsAreWeighted(false),
      _outputRootFile("")
{
    _workspace = std::make_unique<RooWorkspace>("w", "Invariant Mass Fit Workspace");
    std::cout << "[InvariantMassFitter] Initialised (no data yet)" << std::endl;
}

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::setOutputFile(const std::string& filename) {
    _outputRootFile = filename;
    std::cout << "Output ROOT file set to: " << filename << std::endl;
}

void InvariantMassFitter::setReflOverSgn(float const reflOverSgn) {
    _reflOverSgn = reflOverSgn;
}

void InvariantMassFitter::setMassRange(double m1Min, double m1Max,
                                        double m2Min, double m2Max) {
    _m1Min = m1Min; _m1Max = m1Max;
    _m2Min = m2Min; _m2Max = m2Max;
    std::cout << "Mass ranges set: m1=[" << m1Min << ", " << m1Max
              << "], m2=[" << m2Min << ", " << m2Max << "]" << std::endl;
}

void InvariantMassFitter::setFittingFunctions(TString const& sgnFunc,
                                               TString const& bkgFunc,
                                               TString const& reflFunc) {
    _sgnFuncOption  = sgnFunc;
    _bkgFuncOption  = bkgFunc;
    _reflFuncOption = reflFunc;
    std::cout << "Selected signal function:     " << sgnFunc  << std::endl;
    std::cout << "Selected background function: " << bkgFunc  << std::endl;
    std::cout << "Selected reflected function:  " << reflFunc << std::endl;
}

// ─────────────────────────────────────────────────────────────────────────────
// Dataset loading
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::readDatasetsFromFile(TFile* f, bool append)
{
    if (!f || f->IsZombie()) return;
    if (!_m1) createVariables();

    auto loadOne = [&](const char* key, std::unique_ptr<RooDataSet>& target) {
        RooDataSet* ds = dynamic_cast<RooDataSet*>(f->Get(key));
        if (!ds) {
            std::cerr << "[InvariantMassFitter] Key '" << key
                      << "' not found in file\n";
            return;
        }
        if (ds->isWeighted()) _datasetsAreWeighted = true;

        if (!append || !target) {
            // Deep-copy so the object survives after f is closed
            target.reset(static_cast<RooDataSet*>(
                ds->reduce(RooArgSet(*ds->get()))));
            target->SetName(key);
        } else {
            target->append(*ds);
        }
        std::cout << "[InvariantMassFitter] "
                  << (append ? "Merged" : "Loaded")
                  << " '" << key << "': "
                  << target->numEntries() << " entries";
        if (target->isWeighted())
            std::cout << "  sumW=" << target->sumEntries();
        std::cout << std::endl;
    };

    loadOne("osData", _osData);
    loadOne("lsData", _lsData);
}

void InvariantMassFitter::loadDatasets(const std::string& filename)
{
    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[InvariantMassFitter] Cannot open: " << filename << std::endl;
        return;
    }
    _osData.reset();
    _lsData.reset();
    readDatasetsFromFile(f, false);
    f->Close(); delete f;

    if (_osData) _workspace->import(*_osData);
    if (_lsData) _workspace->import(*_lsData);
}

void InvariantMassFitter::mergeDatasets(const std::string& filename)
{
    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[InvariantMassFitter] Cannot open: " << filename << std::endl;
        return;
    }
    readDatasetsFromFile(f, true);
    f->Close(); delete f;

    // Re-import merged datasets, overwriting the previous copies in workspace
    if (_osData) _workspace->import(*_osData, RooFit::RecycleConflictNodes());
    if (_lsData) _workspace->import(*_lsData, RooFit::RecycleConflictNodes());
}

// ─────────────────────────────────────────────────────────────────────────────
// Variables
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::createVariables() {
    _m1 = std::make_unique<RooRealVar>("m1", "m_{1} [GeV/c^{2}]", _m1Min, _m1Max);
    _m2 = std::make_unique<RooRealVar>("m2", "m_{2} [GeV/c^{2}]", _m2Min, _m2Max);
    _m1->setBins(50);
    _m2->setBins(50);

    _y1   = std::make_unique<RooRealVar>("y1",   "y_{1}",    -2.0,  2.0);
    _y2   = std::make_unique<RooRealVar>("y2",   "y_{2}",    -2.0,  2.0);
    _pt1  = std::make_unique<RooRealVar>("pt1",  "p_{T 1}",   0.0, 50.0);
    _pt2  = std::make_unique<RooRealVar>("pt2",  "p_{T 2}",   0.0, 50.0);
    _pairPt   = std::make_unique<RooRealVar>("pairPt",   "p_{T pair}", 0.0, 100.0);
    _pairMass = std::make_unique<RooRealVar>("pairMass", "m_{pair}",   0.0, 100.0);
    _phi1 = std::make_unique<RooRealVar>("phi1", "#phi_{1}",
                                          -2*TMath::Pi(), 2*TMath::Pi());
    _phi2 = std::make_unique<RooRealVar>("phi2", "#phi_{2}",
                                          -2*TMath::Pi(), 2*TMath::Pi());

    _y1->setBins(50);   _y2->setBins(50);
    _phi1->setBins(50); _phi2->setBins(50);
    _pt1->setBins(100); _pt2->setBins(100);
    _pairPt->setBins(100); _pairMass->setBins(100);

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

// ─────────────────────────────────────────────────────────────────────────────
// PDF factories
// ─────────────────────────────────────────────────────────────────────────────

RooAbsPdf* InvariantMassFitter::createSignalPdf(const std::string& name,
                                                  RooRealVar& mass) {
    if (_sgnFuncOption == "gaussian")
        return createGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "doubleGaussian")
        return createDoubleGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "crystalball")
        return createCrystalBallSignal(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown signal option: " + _sgnFuncOption);
}

RooAbsPdf* InvariantMassFitter::createBackgroundPdf(const std::string& name,
                                                      RooRealVar& mass) {
    if (_bkgFuncOption == "chebychev")
        return createChebychevBackground(_workspace.get(), name, mass);
    else if (_bkgFuncOption == "exponential")
        return createExponentialBackground(_workspace.get(), name, mass);
    else if (_bkgFuncOption == "polynomial")
        return createPolynomialBackground(_workspace.get(), name, mass, 2);
    else if (_bkgFuncOption == "poly3")
        return createPolynomialBackground(_workspace.get(), name, mass, 3);
    else
        throw std::runtime_error("Unknown background option: " + _bkgFuncOption);
}

RooAbsPdf* InvariantMassFitter::createReflectedPdf(const std::string& name,
                                                     RooRealVar& mass) {
    if (_reflFuncOption == "gaussian")
        return createGaussianReflected(_workspace.get(), name, mass);
    else if (_reflFuncOption == "doubleGaussian")
        return createDoubleGaussianReflected(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown reflected option: " + _reflFuncOption);
}

// ─────────────────────────────────────────────────────────────────────────────
// Model
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::buildModel(bool const doReflections) {
    if (!_m1 || !_m2) createVariables();

    RooAbsPdf* sig1 = createSignalPdf("m1_sig", *_m1);
    RooAbsPdf* bkg1 = createBackgroundPdf("m1_bkg", *_m1);
    RooAbsPdf* sig2 = createSignalPdf("m2_sig", *_m2);
    RooAbsPdf* bkg2 = createBackgroundPdf("m2_bkg", *_m2);

    int n_events = _lsData->numEntries();
    double initBkg = 0.7 * std::sqrt(n_events);
    double initSgn = 0.2 * std::sqrt(n_events);
    std::cout << "Initial parameters, bkg: " << initBkg
              << " sgn: " << initSgn << std::endl;

    RooRealVar* nBkg1 = new RooRealVar("nBkg1", "background yield of cand 1",
                                        initBkg, 0.0, 1e6);
    RooRealVar* nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1",
                                        initSgn, 0.0, 1e5);
    RooFormulaVar* nRefl1 = new RooFormulaVar("nRefl1",
        "reflected signal yield of cand1", "@0 * @1",
        RooArgList(_reflOverSgn, *nSgn1));

    RooRealVar* nBkg2 = new RooRealVar("nBkg2", "background yield of cand2",
                                        initBkg, 0.0, 1e6);
    RooRealVar* nSgn2 = new RooRealVar("nSgn2", "signal yield of cand2",
                                        initSgn, 0.0, 1e5);
    RooFormulaVar* nRefl2 = new RooFormulaVar("nRefl2",
        "reflected signal yield of cand 2", "@0 * @1",
        RooArgList(_reflOverSgn, *nSgn2));

    RooFormulaVar nBkgBkg("nBkgBkg", "fraction of bkg x bkg",
                           "@0 * @1", RooArgList(*nBkg1, *nBkg2));
    RooFormulaVar nBkgSgn("nBkgSgn", "fraction of bkg x sgn",
                           "@0 * @1", RooArgList(*nBkg1, *nSgn2));
    RooFormulaVar nSgnBkg("nSgnBkg", "fraction of sgn x bkg",
                           "@0 * @1", RooArgList(*nSgn1, *nBkg2));
    RooFormulaVar nBkgRefl("nBkgRefl", "fraction of bkg x refl",
                            "@0 * @1", RooArgList(*nBkg1, *nRefl2));
    RooFormulaVar nReflBkg("nReflBkg", "fraction of refl x bkg",
                            "@0 * @1", RooArgList(*nRefl1, *nBkg2));
    RooFormulaVar nSgnSgn("nSgnSgn", "fraction of sgn x sgn",
                           "@0 * @1", RooArgList(*nSgn1, *nSgn2));
    RooFormulaVar nReflRefl("nReflRefl", "fraction of refl x refl",
                             "@0 * @1", RooArgList(*nRefl1, *nRefl2));
    RooFormulaVar nSgnRefl("nSgnRefl", "fraction of sgn x refl",
                            "@0 * @1", RooArgList(*nSgn1, *nRefl2));
    RooFormulaVar nReflSgn("nReflSgn", "fraction of refl x sgn",
                            "@0 * @1", RooArgList(*nRefl1, *nSgn2));

    RooProdPdf pdfSigSig("pdfSigSig", "Signal-Signal",     RooArgList(*sig1, *sig2));
    RooProdPdf pdfSigBkg("pdfSigBkg", "Signal-Background", RooArgList(*sig1, *bkg2));
    RooProdPdf pdfBkgSig("pdfBkgSig", "Background-Signal", RooArgList(*bkg1, *sig2));
    RooProdPdf pdfBkgBkg("pdfBkgBkg", "Background-Background", RooArgList(*bkg1, *bkg2));

    RooAbsPdf* refl1 = nullptr;
    RooAbsPdf* refl2 = nullptr;
    RooProdPdf* pdfSigRefl  = nullptr;
    RooProdPdf* pdfReflSig  = nullptr;
    RooProdPdf* pdfReflBkg  = nullptr;
    RooProdPdf* pdfBkgRefl  = nullptr;
    RooProdPdf* pdfReflRefl = nullptr;

    if (doReflections) {
        refl1 = createReflectedPdf("m1_refl", *_m1);
        refl2 = createReflectedPdf("m2_refl", *_m2);
        pdfSigRefl  = new RooProdPdf("pdfSigRefl",  "Signal-Reflected",
                                      RooArgList(*sig1,  *refl2));
        pdfReflSig  = new RooProdPdf("pdfReflSig",  "Reflected-Signal",
                                      RooArgList(*refl1, *sig2));
        pdfReflBkg  = new RooProdPdf("pdfReflBkg",  "Reflected-Background",
                                      RooArgList(*refl1, *bkg2));
        pdfBkgRefl  = new RooProdPdf("pdfBkgRefl",  "Background-Reflected",
                                      RooArgList(*bkg1,  *refl2));
        pdfReflRefl = new RooProdPdf("pdfReflRefl", "Reflected-Reflected",
                                      RooArgList(*refl1, *refl2));
    }

    RooAddPdf* totalPdf = nullptr;
    if (doReflections) {
        totalPdf = new RooAddPdf("model",
            "background + signal + reflection pdf 2D",
            RooArgList(pdfSigSig, *pdfSigRefl, *pdfReflSig, pdfSigBkg,
                       pdfBkgSig, pdfBkgBkg, *pdfBkgRefl, *pdfReflBkg, *pdfReflRefl),
            RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg,
                       nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));
    } else {
        totalPdf = new RooAddPdf("model", "background + signal 2D model",
            RooArgList(pdfSigSig, pdfSigBkg, pdfBkgSig, pdfBkgBkg),
            RooArgList(nSgnSgn, nSgnBkg, nBkgSgn, nBkgBkg));
    }

    _workspace->import(*totalPdf);
    std::cout << "2D model built" << std::endl;
}

// ─────────────────────────────────────────────────────────────────────────────
// Prefit
// ─────────────────────────────────────────────────────────────────────────────

RooFitResult* InvariantMassFitter::performPrefit1D(RooDataSet* data,
                                                    RooRealVar& mass,
                                                    const std::string& name,
                                                    bool useReflections) {
    std::cout << "\n=== Performing 1D prefit for " << name << " ===" << std::endl;

    RooAbsPdf* sig = createSignalPdf(name + "_sig", mass);
    RooAbsPdf* bkg = createBackgroundPdf(name + "_bkg", mass);
    double nEntries = data->numEntries();

    RooRealVar* nSig = new RooRealVar((name + "_nSig").c_str(),
                                      "Signal yield", 0.1*nEntries, 0.0, 1e10);
    RooRealVar* nBkg = new RooRealVar((name + "_nBkg").c_str(),
                                      "Background yield", 0.5*nEntries, 0.0, 1e10);

    RooAddPdf* model = nullptr;

    if (useReflections) {
        RooAbsPdf* refl = createReflectedPdf(name + "_refl", mass);
        RooRealVar* reflOverSgn = new RooRealVar((name + "_reflOverSgn").c_str(),
                                                  "Reflection/Signal ratio",
                                                  _reflOverSgn, 0.0, 1.0);
        reflOverSgn->setConstant(true);
        RooFormulaVar* nRefl = new RooFormulaVar((name + "_nRefl").c_str(),
                                                  "Reflection yield", "@0 * @1",
                                                  RooArgList(*reflOverSgn, *nSig));
        model = new RooAddPdf((name + "_model").c_str(),
                              "Signal + Reflection + Background",
                              RooArgList(*sig, *refl, *bkg),
                              RooArgList(*nSig, *nRefl, *nBkg));
    } else {
        model = new RooAddPdf((name + "_model").c_str(),
                              "Signal + Background",
                              RooArgList(*sig, *bkg),
                              RooArgList(*nSig, *nBkg));
    }

    RooFitResult* result = nullptr;
    if (_datasetsAreWeighted) {
        result = model->fitTo(*data, RooFit::Save(),
                              Extended(true),
                              AsymptoticError(true),
                              Strategy(2),
                              Optimize(2),
                              Offset(true),
                              Minimizer("Minuit2", "migrad"),
                              PrintLevel(1),
                              RooFit::MaxCalls(100000));
    } else {
        result = model->fitTo(*data, RooFit::Save(),
                              RooFit::Extended(true));
    }
    result->Print("v");

    std::cout << "Prefit results for " << name << ":\n";
    std::cout << "  Signal yield:     " << nSig->getVal()
              << " +/- " << nSig->getError() << "\n";
    std::cout << "  Background yield: " << nBkg->getVal()
              << " +/- " << nBkg->getError() << "\n";

    _workspace->import(*model);
    _workspace->import(*result, (name + "_prefitResult").c_str());

    return result;
}

void InvariantMassFitter::performPrefit(const std::string& dataType,
                                         bool useReflections) {
    std::cout << "\n========================================\n"
              << "=== Starting 1D Prefit Procedure ===\n"
              << "========================================\n";

    if (!_osData && !_lsData) {
        std::cerr << "Error: No data loaded. Call loadDatasets() first.\n";
        return;
    }
    if (!_m1 || !_m2) createVariables();

    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n--- OS Data Prefits ---\n";
        _prefitResult_os_m1.reset(
            performPrefit1D(_osData.get(), *_m1, "prefit_os_m1", useReflections));
        _prefitResult_os_m2.reset(
            performPrefit1D(_osData.get(), *_m2, "prefit_os_m2", useReflections));
    }

    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n--- LS Data Prefits ---\n";
        _prefitResult_ls_m1.reset(
            performPrefit1D(_lsData.get(), *_m1, "prefit_ls_m1", useReflections));
        _prefitResult_ls_m2.reset(
            performPrefit1D(_lsData.get(), *_m2, "prefit_ls_m2", useReflections));
    }

    std::cout << "\n========================================\n"
              << "=== 1D Prefit Procedure Complete ===\n"
              << "========================================\n";
}

void InvariantMassFitter::extractPrefitParameters(const std::string& prefitName,
                                                   const std::string& targetName) {
    std::cout << "  Extracting prefit parameters from "
              << prefitName << " to " << targetName << std::endl;

    RooArgSet allVars = _workspace->allVars();
    int nParamsSet = 0;

    for (RooAbsArg* arg : allVars) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
        if (!var) continue;

        std::string varName = var->GetName();
        if (varName.find(prefitName) == std::string::npos) continue;

        size_t pos = varName.find(prefitName);
        std::string suffix = varName.substr(pos + prefitName.length());

        // Skip yield parameters
        if (suffix.empty()
            || varName.find("_nBkg")       != std::string::npos
            || varName.find("_nSig")       != std::string::npos
            || varName.find("_nRefl")      != std::string::npos
            || varName.find("_reflOverSgn")!= std::string::npos)
            continue;

        std::string targetParamName = targetName + suffix;
        std::cout << targetParamName << std::endl;

        RooRealVar* targetParam = _workspace->var(targetParamName.c_str());
        if (targetParam) {
            targetParam->setVal(var->getVal());
            if (varName.find("_refl") != std::string::npos) {
                targetParam->setConstant(true);
                std::cout << "    " << targetParamName
                          << " set CONSTANT (_refl detected)\n";
            } else {
                std::cout << "    " << targetParamName
                          << " set to: " << var->getVal() << "\n";
            }
            ++nParamsSet;
        }
    }

    if (nParamsSet == 0)
        std::cout << "    Warning: No background parameters found to transfer\n";
    else
        std::cout << "    Successfully set " << nParamsSet
                  << " background parameter(s)\n";
}

void InvariantMassFitter::setPrefitParameters(bool const doReflections) {
    std::cout << "\n=== Setting 2D fit initial parameters from prefit ===\n";

    auto getParam = [this](const std::string& name) -> RooRealVar* {
        return _workspace->var(name.c_str());
    };

    if (_prefitResult_os_m1) {
        std::cout << "\n--- Setting m1 parameters from OS prefit ---\n";
        extractPrefitParameters("prefit_os_m1_sig", "m1_sig");
        extractPrefitParameters("prefit_os_m1_bkg", "m1_bkg");
        if (doReflections)
            extractPrefitParameters("prefit_os_m1_refl", "m1_refl");

        RooRealVar* prefit_nSig = getParam("prefit_os_m1_nSig");
        RooRealVar* prefit_nBkg = getParam("prefit_os_m1_nBkg");
        if (prefit_nSig && prefit_nBkg) {
            RooRealVar* nSgn1 = getParam("nSgn1");
            RooRealVar* nBkg1 = getParam("nBkg1");
            if (nSgn1) {
                double v = std::sqrt(prefit_nSig->getVal());
                nSgn1->setVal(v);
                std::cout << "  nSgn1 set to: " << v << "\n";
            }
            if (nBkg1) {
                double v = std::sqrt(prefit_nBkg->getVal());
                nBkg1->setVal(v);
                std::cout << "  nBkg1 set to: " << v << "\n";
            }
        }
    }

    if (_prefitResult_os_m2) {
        std::cout << "\n--- Setting m2 parameters from OS prefit ---\n";
        extractPrefitParameters("prefit_os_m2_sig", "m2_sig");
        extractPrefitParameters("prefit_os_m2_bkg", "m2_bkg");
        if (doReflections)
            extractPrefitParameters("prefit_os_m2_refl", "m2_refl");

        RooRealVar* prefit_nSig = getParam("prefit_os_m2_nSig");
        RooRealVar* prefit_nBkg = getParam("prefit_os_m2_nBkg");
        if (prefit_nSig && prefit_nBkg) {
            RooRealVar* nSgn2 = getParam("nSgn2");
            RooRealVar* nBkg2 = getParam("nBkg2");
            if (nSgn2) {
                double v = std::sqrt(prefit_nSig->getVal());
                nSgn2->setVal(v);
                std::cout << "  nSgn2 set to: " << v << "\n";
            }
            if (nBkg2) {
                double v = std::sqrt(prefit_nBkg->getVal());
                nBkg2->setVal(v);
                std::cout << "  nBkg2 set to: " << v << "\n";
            }
        }
    }

    std::cout << "\n=== Initial parameters set from prefit ===\n";
}

void InvariantMassFitter::plotPrefitResults(const std::string& outputPrefix) {
    std::cout << "\n=== Plotting prefit results ===\n";

    TCanvas* c = new TCanvas("c_prefit", "Prefit Results", 1600, 1200);
    c->Divide(2, 2);

    auto plotOne = [&](int pad, RooDataSet* data, RooRealVar& mass,
                       const std::string& modelName,
                       const std::string& sigComp,
                       const std::string& bkgComp,
                       const std::string& title) {
        if (!data) return;
        c->cd(pad);
        RooPlot* frame = mass.frame(RooFit::Title(title.c_str()));
        data->plotOn(frame, RooFit::DrawOption("pz"));
        RooAbsPdf* model = _workspace->pdf(modelName.c_str());
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components(sigComp.c_str()),
                          DrawOption("F"), LineWidth(2), LineStyle(1),
                          LineColor(kGreen+2), FillColor(kGreen-2));
            model->plotOn(frame, RooFit::Components(bkgComp.c_str()),
                          RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        }
        frame->Draw();
    };

    if (_prefitResult_os_m1)
        plotOne(1, _osData.get(), *_m1, "prefit_os_m1_model",
                "prefit_os_m1_sig_gauss", "prefit_os_m1_bkg_exp",
                "OS m_{1} Prefit");
    if (_prefitResult_os_m2)
        plotOne(2, _osData.get(), *_m2, "prefit_os_m2_model",
                "prefit_os_m2_sig_gauss", "prefit_os_m2_bkg_exp",
                "OS m_{2} Prefit");
    if (_prefitResult_ls_m1)
        plotOne(3, _lsData.get(), *_m1, "prefit_ls_m1_model",
                "prefit_ls_m1_sig_gauss", "prefit_ls_m1_bkg_exp",
                "LS m_{1} Prefit");
    if (_prefitResult_ls_m2)
        plotOne(4, _lsData.get(), *_m2, "prefit_ls_m2_model",
                "prefit_ls_m2_sig_gauss", "prefit_ls_m2_bkg_exp",
                "LS m_{2} Prefit");

    c->SaveAs((outputPrefix + "_prefit.pdf").c_str());
    c->SaveAs((outputPrefix + "_prefit.png").c_str());
    std::cout << "Prefit plots saved to " << outputPrefix << "_prefit.pdf/png\n";
    delete c;
}

// ─────────────────────────────────────────────────────────────────────────────
// Main fit
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::fit(const std::string& dataType) {
    if (!_workspace->pdf("model")) {
        std::cerr << "Error: Model not built. Call buildModel() first.\n";
        return;
    }
    if (!_osData && !_lsData) {
        std::cerr << "Error: No data loaded. Call loadDatasets() first.\n";
        return;
    }

    RooAbsPdf* pdf = _workspace->pdf("model");

    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n=== Fitting OS data ===\n";
        RooFitResult* result = nullptr;
        if (_datasetsAreWeighted) {
            result = pdf->fitTo(*_osData, RooFit::Save(),
                                Extended(true),
                                AsymptoticError(true),
                                Strategy(2),
                                Optimize(2),
                                Offset(true),
                                Minimizer("Minuit2", "migrad"),
                                PrintLevel(1),
                                RooFit::MaxCalls(100000));
        } else {
            result = pdf->fitTo(*_osData, RooFit::Save(),
                                RooFit::Extended(true),
                                RooFit::PrintLevel(1));
        }
        result->Print("v");
        _workspace->import(*result, "osFitResult");
        RooArgSet* params = pdf->getParameters(*_osData);
        _workspace->saveSnapshot("osParams", *params);
        delete params;
        delete result;
    }

    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n=== Fitting LS data ===\n";
        RooFitResult* result = nullptr;
        if (_datasetsAreWeighted) {
            result = pdf->fitTo(*_lsData, RooFit::Save(),
                                Extended(true),
                                AsymptoticError(true),
                                Strategy(2),
                                Optimize(2),
                                Offset(true),
                                Minimizer("Minuit2", "migrad"),
                                PrintLevel(1),
                                RooFit::MaxCalls(100000));
        } else {
            result = pdf->fitTo(*_lsData, RooFit::Save(),
                                RooFit::Extended(true),
                                RooFit::PrintLevel(1));
        }
        result->Print("v");
        _workspace->import(*result, "lsFitResult");
        RooArgSet* params = pdf->getParameters(*_lsData);
        _workspace->saveSnapshot("lsParams", *params);
        delete params;
        delete result;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Plotting helpers
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::plotModelOnFrame(RooPlot* frame, RooAbsPdf* pdf,
                                            bool const doReflections) {
    pdf->plotOn(frame, Name("model_curve"));
    pdf->plotOn(frame, Components("pdfBkgBkg"),
                LineWidth(2), LineStyle(kDashed), LineColor(kRed),
                Name("bkgBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfSigBkg"),
                LineStyle(3), LineWidth(2), LineColor(kOrange+6),
                Name("sgnBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfBkgSig"),
                LineStyle(4), LineWidth(2), LineColor(kOrange),
                Name("bkgSgnFunc2D_curve"));
    pdf->plotOn(frame, DrawOption("F"), Components("pdfSigSig"),
                LineWidth(2), LineStyle(1), LineColor(kGreen+2),
                FillColor(kGreen-2), Name("sgnSgnFunc2D_curve"));

    if (doReflections) {
        pdf->plotOn(frame, Components("pdfBkgRefl"),
                    LineWidth(2), LineStyle(2), LineColor(kViolet),
                    Name("bkgReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfSigRefl"),
                    LineWidth(2), LineStyle(2), LineColor(kMagenta-7),
                    Name("sgnReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflBkg"),
                    LineWidth(2), LineStyle(2), LineColor(kCyan-3),
                    Name("reflBkgFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflSig"),
                    LineWidth(2), LineStyle(2), LineColor(kMagenta-5),
                    Name("reflSgnFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflRefl"),
                    LineWidth(2), LineStyle(2), LineColor(kAzure+3),
                    Name("reflReflFunc2D_curve"));
    }
    frame->Draw();
}

double InvariantMassFitter::computeChi2(RooPlot* frame, RooDataSet* data,
                                         int nFitParams) {
    double chi2Ndf = frame->chiSquare("model_curve", "data_cand1", nFitParams);
    int nBins = frame->GetXaxis()->GetNbins();
    int ndf   = nBins - nFitParams;
    double chi2 = chi2Ndf * ndf;

    std::cout << "  Chi2 / ndf = " << chi2 << " / " << ndf
              << " = " << chi2Ndf << std::endl;

    if (_datasetsAreWeighted) {
        std::cout << "  (weighted fit: N_eff = " << data->numEntries()
                  << ", sumW = " << data->sumEntries() << ")\n";
    }
    return chi2Ndf;
}

void InvariantMassFitter::plotFit(const std::string& outputPrefix,
                                   bool const doReflections) {
    if (!_workspace->pdf("model")) {
        std::cerr << "Error: No model to plot\n";
        return;
    }

    RooAbsPdf* pdf = _workspace->pdf("model");
    TCanvas* c = new TCanvas("c_fit", "Fit Results", 1600, 800);
    c->Divide(2, 2);

    double binWidth = 0.005;
    int nBins = (int)((_m1Max - _m1Min) / binWidth);
    _m1->setBins(nBins);
    _m2->setBins(nBins);

    // Count floating parameters
    RooArgSet* allParams = pdf->getParameters(*_osData);
    int nFitParams = 0;
    for (RooAbsArg* arg : *allParams) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
        if (var && !var->isConstant()) ++nFitParams;
    }
    delete allParams;
    std::cout << "  Floating fit parameters: " << nFitParams << std::endl;

    if (_osData) {
        _workspace->loadSnapshot("osParams");
        c->cd(1);
        RooPlot* f1 = _m1->frame(RooFit::Title("OS: m_{1} projection"));
        _osData->plotOn(f1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(f1, pdf, doReflections);
        std::cout << "OS m1: "; computeChi2(f1, _osData.get(), nFitParams);

        c->cd(2);
        RooPlot* f2 = _m2->frame(RooFit::Title("OS: m_{2} projection"));
        _osData->plotOn(f2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(f2, pdf, doReflections);
        std::cout << "OS m2: "; computeChi2(f2, _osData.get(), nFitParams);
    }

    if (_lsData) {
        _workspace->loadSnapshot("lsParams");
        c->cd(3);
        RooPlot* f3 = _m1->frame(RooFit::Title("LS: m_{1} projection"));
        _lsData->plotOn(f3, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(f3, pdf, doReflections);
        std::cout << "LS m1: "; computeChi2(f3, _lsData.get(), nFitParams);

        c->cd(4);
        RooPlot* f4 = _m2->frame(RooFit::Title("LS: m_{2} projection"));
        _lsData->plotOn(f4, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(f4, pdf, doReflections);
        std::cout << "LS m2: "; computeChi2(f4, _lsData.get(), nFitParams);
    }

    c->SaveAs((outputPrefix + ".pdf").c_str());
    c->SaveAs((outputPrefix + ".png").c_str());
    std::cout << "Plots saved to " << outputPrefix << ".pdf/png\n";
    delete c;
}

// ─────────────────────────────────────────────────────────────────────────────
// Save results
// ─────────────────────────────────────────────────────────────────────────────

void InvariantMassFitter::saveResultsToFile(bool const doReflections) {
    if (_outputRootFile.empty()) {
        std::cerr << "Error: Output file not set. Call setOutputFile() first.\n";
        return;
    }

    TFile* outFile = TFile::Open(_outputRootFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << _outputRootFile << "\n";
        return;
    }
    std::cout << "Saving results to " << _outputRootFile << "\n";

    _workspace->Write();

    RooAbsPdf* pdf = _workspace->pdf("model");

    if (_osData && pdf) {
        _workspace->loadSnapshot("osParams");

        RooPlot* frame_os_m1 = _m1->frame(RooFit::Title("OS: m_{1} projection"));
        _osData->plotOn(frame_os_m1, RooFit::Name("data_os_m1"),
                        RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_os_m1, pdf, doReflections);
        TH1* h1 = _osData->createHistogram("h_data_os_m1", *_m1);
        h1->SetTitle("OS Data m_{1};m_{1} [GeV/c^{2}];Entries");
        h1->Write(); frame_os_m1->Write("frame_os_m1");

        RooPlot* frame_os_m2 = _m2->frame(RooFit::Title("OS: m_{2} projection"));
        _osData->plotOn(frame_os_m2, RooFit::Name("data_os_m2"),
                        RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_os_m2, pdf, doReflections);
        TH1* h2 = _osData->createHistogram("h_data_os_m2", *_m2);
        h2->SetTitle("OS Data m_{2};m_{2} [GeV/c^{2}];Entries");
        h2->Write(); frame_os_m2->Write("frame_os_m2");
    }

    if (_lsData && pdf) {
        _workspace->loadSnapshot("lsParams");

        RooPlot* frame_ls_m1 = _m1->frame(RooFit::Title("LS: m_{1} projection"));
        _lsData->plotOn(frame_ls_m1, RooFit::Name("data_ls_m1"),
                        RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_ls_m1, pdf, doReflections);
        TH1* h3 = _lsData->createHistogram("h_data_ls_m1", *_m1);
        h3->SetTitle("LS Data m_{1};m_{1} [GeV/c^{2}];Entries");
        h3->Write(); frame_ls_m1->Write("frame_ls_m1");

        RooPlot* frame_ls_m2 = _m2->frame(RooFit::Title("LS: m_{2} projection"));
        _lsData->plotOn(frame_ls_m2, RooFit::Name("data_ls_m2"),
                        RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_ls_m2, pdf, doReflections);
        TH1* h4 = _lsData->createHistogram("h_data_ls_m2", *_m2);
        h4->SetTitle("LS Data m_{2};m_{2} [GeV/c^{2}];Entries");
        h4->Write(); frame_ls_m2->Write("frame_ls_m2");
    }

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

    RooFitResult* osRes = nullptr;
    RooFitResult* lsRes = nullptr;
    if (_workspace->obj("osFitResult")) {
        osRes = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
        if (osRes) osRes->Write("osFitResult");
    }
    if (_workspace->obj("lsFitResult")) {
        lsRes = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));
        if (lsRes) lsRes->Write("lsFitResult");
    }

    // Signal parameters for 3-sigma range
    RooRealVar* m1_mean  = _workspace->var("m1_sig_mean");
    RooRealVar* m1_sigma = _workspace->var("m1_sig_sigma");
    RooRealVar* m2_mean  = _workspace->var("m2_sig_mean");
    RooRealVar* m2_sigma = _workspace->var("m2_sig_sigma");

    if (!m1_mean || !m1_sigma || !m2_mean || !m2_sigma) {
        std::cerr << "Error: Cannot find signal parameters in workspace\n";
        outFile->Close(); delete outFile; return;
    }

    _m1->setRange("3sigmaRange",
                  m1_mean->getVal() - 3*m1_sigma->getVal(),
                  m1_mean->getVal() + 3*m1_sigma->getVal());
    _m2->setRange("3sigmaRange",
                  m2_mean->getVal() - 3*m2_sigma->getVal(),
                  m2_mean->getVal() + 3*m2_sigma->getVal());

    RooProdPdf*    pdfSigSig  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfSigSig"));
    RooProdPdf*    pdfSigBkg  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfSigBkg"));
    RooProdPdf*    pdfBkgSig  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfBkgSig"));
    RooProdPdf*    pdfBkgBkg  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfBkgBkg"));
    RooFormulaVar* nSgnSgn    = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnSgn"));
    RooFormulaVar* nSgnBkg    = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnBkg"));
    RooFormulaVar* nBkgSgn    = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgSgn"));
    RooFormulaVar* nBkgBkg_fv = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgBkg"));

    RooProdPdf*    pdfSigRefl  = nullptr;
    RooProdPdf*    pdfReflSig  = nullptr;
    RooProdPdf*    pdfBkgRefl  = nullptr;
    RooProdPdf*    pdfReflBkg  = nullptr;
    RooProdPdf*    pdfReflRefl = nullptr;
    RooFormulaVar* nSgnRefl    = nullptr;
    RooFormulaVar* nReflSgn    = nullptr;
    RooFormulaVar* nBkgRefl    = nullptr;
    RooFormulaVar* nReflBkg    = nullptr;
    RooFormulaVar* nReflRefl   = nullptr;

    if (doReflections) {
        pdfSigRefl  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfSigRefl"));
        pdfReflSig  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfReflSig"));
        pdfBkgRefl  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfBkgRefl"));
        pdfReflBkg  = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfReflBkg"));
        pdfReflRefl = dynamic_cast<RooProdPdf*>  (_workspace->pdf("pdfReflRefl"));
        nSgnRefl    = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnRefl"));
        nReflSgn    = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflSgn"));
        nBkgRefl    = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgRefl"));
        nReflBkg    = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflBkg"));
        nReflRefl   = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflRefl"));
    }

    // Helper: compute and print yields, return them
    auto calcYields = [&](RooFitResult* res, const std::string& label)
        -> std::array<RooRealVar*, 9>
    {
        std::array<RooRealVar*, 9> y{};
        if (!res || !pdfSigSig || !nSgnSgn) return y;

        std::cout << "\n=== " << label << " Yields in 3-sigma range ===\n";
        y[0] = getYieldInRange(res, pdfSigSig, nSgnSgn,    "3sigmaRange");
        y[1] = getYieldInRange(res, pdfSigBkg, nSgnBkg,    "3sigmaRange");
        y[2] = getYieldInRange(res, pdfBkgSig, nBkgSgn,    "3sigmaRange");
        y[3] = getYieldInRange(res, pdfBkgBkg, nBkgBkg_fv, "3sigmaRange");
        std::cout << " SgnSgn: " << y[0]->getVal() << " +- " << y[0]->getError() << "\n";
        std::cout << " SgnBkg: " << y[1]->getVal() << " +- " << y[1]->getError() << "\n";
        std::cout << " BkgSgn: " << y[2]->getVal() << " +- " << y[2]->getError() << "\n";
        std::cout << " BkgBkg: " << y[3]->getVal() << " +- " << y[3]->getError() << "\n";

        if (doReflections && pdfSigRefl) {
            y[4] = getYieldInRange(res, pdfSigRefl,  nSgnRefl,  "3sigmaRange");
            y[5] = getYieldInRange(res, pdfReflSig,  nReflSgn,  "3sigmaRange");
            y[6] = getYieldInRange(res, pdfBkgRefl,  nBkgRefl,  "3sigmaRange");
            y[7] = getYieldInRange(res, pdfReflBkg,  nReflBkg,  "3sigmaRange");
            y[8] = getYieldInRange(res, pdfReflRefl, nReflRefl, "3sigmaRange");
            std::cout << " SgnRefl:  " << y[4]->getVal() << " +- " << y[4]->getError() << "\n";
            std::cout << " ReflSgn:  " << y[5]->getVal() << " +- " << y[5]->getError() << "\n";
            std::cout << " BkgRefl:  " << y[6]->getVal() << " +- " << y[6]->getError() << "\n";
            std::cout << " ReflBkg:  " << y[7]->getVal() << " +- " << y[7]->getError() << "\n";
            std::cout << " ReflRefl: " << y[8]->getVal() << " +- " << y[8]->getError() << "\n";
        }
        return y;
    };

    if (osRes) _workspace->loadSnapshot("osParams");
    auto yOS = calcYields(osRes, "OS");

    if (lsRes) _workspace->loadSnapshot("lsParams");
    auto yLS = calcYields(lsRes, "LS");

    // --- Fit-parameter histograms ---
    auto makeFitParamsHisto = [&](const std::string& name,
                                  RooDataSet* data,
                                  const std::string& snapshotName) {
        _workspace->loadSnapshot(snapshotName.c_str());
        RooRealVar* m1m  = _workspace->var("m1_sig_mean");
        RooRealVar* m1s  = _workspace->var("m1_sig_sigma");
        RooRealVar* m2m  = _workspace->var("m2_sig_mean");
        RooRealVar* m2s  = _workspace->var("m2_sig_sigma");
        if (!m1m || !m1s || !m2m || !m2s) {
            std::cerr << "Warning: signal parameters not found for " << name << "\n";
            return;
        }

        int nBinsF = (int)((_m1Max - _m1Min) / 0.005);
        _m1->setBins(nBinsF); _m2->setBins(nBinsF);
        RooAbsPdf* pdfModel = _workspace->pdf("model");

        // Count floating params
        RooArgSet* ap = pdfModel->getParameters(*data);
        int nFl = 0;
        for (RooAbsArg* a : *ap) {
            RooRealVar* v = dynamic_cast<RooRealVar*>(a);
            if (v && !v->isConstant()) ++nFl;
        }
        delete ap;

        RooPlot* f1 = _m1->frame();
        data->plotOn(f1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f1, RooFit::Name("model_curve"));
        double chi2_m1 = f1->chiSquare("model_curve", "data_cand1", nFl);

        RooPlot* f2 = _m2->frame();
        data->plotOn(f2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f2, RooFit::Name("model_curve"));
        double chi2_m2 = f2->chiSquare("model_curve", "data_cand1", nFl);

        auto makeH = [](const std::string& hname, const std::string& title,
                        double val, double err) -> TH1D* {
            TH1D* h = new TH1D(hname.c_str(), title.c_str(), 1, 0, 1);
            h->SetBinContent(1, val); h->SetBinError(1, err);
            return h;
        };
        std::vector<TH1D*> hs = {
            makeH(name+"_m1Mean",    "m1 mean [GeV/c^{2}]",   m1m->getVal(), m1m->getError()),
            makeH(name+"_m1Sigma",   "m1 sigma [GeV/c^{2}]",  m1s->getVal(), m1s->getError()),
            makeH(name+"_m2Mean",    "m2 mean [GeV/c^{2}]",   m2m->getVal(), m2m->getError()),
            makeH(name+"_m2Sigma",   "m2 sigma [GeV/c^{2}]",  m2s->getVal(), m2s->getError()),
            makeH(name+"_chi2ndf_m1","chi2/ndf m1 projection", chi2_m1,       0.0),
            makeH(name+"_chi2ndf_m2","chi2/ndf m2 projection", chi2_m2,       0.0),
        };
        for (TH1D* h : hs) { h->Write(); delete h; }
        delete f1; delete f2;
    };

    if (_osData && osRes) makeFitParamsHisto("hOS", _osData.get(), "osParams");
    if (_lsData && lsRes) makeFitParamsHisto("hLS", _lsData.get(), "lsParams");

    // --- Yield histograms ---
    auto writeYieldHisto = [&](const std::string& hname, const std::string& title,
                                std::array<RooRealVar*, 9>& y) {
        if (!y[0]) return;
        int nBinsY = doReflections ? 9 : 4;
        TH1D* h = new TH1D(hname.c_str(), title.c_str(), nBinsY, 0, nBinsY);
        const char* labels[9] = {"SgnSgn","SgnBkg","BkgSgn","BkgBkg",
                                  "SgnRefl","ReflSgn","BkgRefl","ReflBkg","ReflRefl"};
        for (int i = 0; i < nBinsY; ++i) {
            h->GetXaxis()->SetBinLabel(i+1, labels[i]);
            if (y[i]) {
                h->SetBinContent(i+1, y[i]->getVal());
                h->SetBinError  (i+1, y[i]->getError());
            }
        }
        h->Write();
        delete h;
    };

    writeYieldHisto("hYields_os", "OS Yields in 3#sigma range", yOS);
    writeYieldHisto("hYields_ls", "LS Yields in 3#sigma range", yLS);

    // Cleanup yield vars
    for (auto& arr : {yOS, yLS})
        for (RooRealVar* v : arr) delete v;

    outFile->Close(); delete outFile;
    std::cout << "Results saved successfully to " << _outputRootFile << "\n";
}

void InvariantMassFitter::saveWorkspace(const std::string& filename) {
    _workspace->writeToFile(filename.c_str());
    std::cout << "Workspace saved to " << filename << "\n";
}

void InvariantMassFitter::printFitResults() {
    std::cout << "\n=== Fit Results Summary ===\n";

    if (_workspace->obj("osFitResult")) {
        std::cout << "\n--- OS Fit ---\n";
        auto* r = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
        if (r) r->Print("v");
    }
    if (_workspace->obj("lsFitResult")) {
        std::cout << "\n--- LS Fit ---\n";
        auto* r = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));
        if (r) r->Print("v");
    }

    RooRealVar* nSigSig = _workspace->var("nSigSig");
    RooRealVar* nSigBkg = _workspace->var("nSigBkg");
    RooRealVar* nBkgSig = _workspace->var("nBkgSig");
    RooRealVar* nBkgBkg = _workspace->var("nBkgBkg");
    if (nSigSig && nSigBkg && nBkgSig && nBkgBkg) {
        std::cout << "\n=== Yield Summary ===\n";
        std::cout << "Signal-Signal:         " << nSigSig->getVal()
                  << " +/- " << nSigSig->getError() << "\n";
        std::cout << "Signal-Background:     " << nSigBkg->getVal()
                  << " +/- " << nSigBkg->getError() << "\n";
        std::cout << "Background-Signal:     " << nBkgSig->getVal()
                  << " +/- " << nBkgSig->getError() << "\n";
        std::cout << "Background-Background: " << nBkgBkg->getVal()
                  << " +/- " << nBkgBkg->getError() << "\n";
        if (_datasetsAreWeighted)
            std::cout << "\n(Yields are efficiency-corrected)\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Yield integration helper
// ─────────────────────────────────────────────────────────────────────────────

RooRealVar* InvariantMassFitter::getYieldInRange(RooFitResult* fitResult,
                                                  RooProdPdf* function,
                                                  RooFormulaVar* nCands,
                                                  TString range)
{
    RooAbsReal* intRange = function->createIntegral(RooArgSet(*_m1, *_m2),
                                                     Range(range));
    RooAbsReal* intTotal = function->createIntegral(RooArgSet(*_m1, *_m2));

    double yieldRange = nCands->getVal() * intRange->getVal() / intTotal->getVal();
    double nCandsErr  = nCands->getPropagatedError(*fitResult);
    double errorRange = nCandsErr * intRange->getVal() / intTotal->getVal();

    RooRealVar* yield = new RooRealVar("yield", "Yield in 3-sigma range", yieldRange);
    yield->setError(errorRange);

    delete intRange;
    delete intTotal;
    return yield;
}