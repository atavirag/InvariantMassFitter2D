// InvariantMassFitter.cxx
//TODO:
// - Add prefit stage
// - Add more fitting functions

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

#include "../01-prefiltering/CandidatePair.h"
#include "functions/BackgroundFunctions.h"
#include "functions/SignalFunctions.h"
#include "functions/ReflectedFunctions.h"

using namespace RooFit;

InvariantMassFitter::InvariantMassFitter(const std::string& osFilename, const std::string& lsFilename)
    : _osTree(nullptr), _lsTree(nullptr),
      _m1Min(1.7), _m1Max(2.1), _m2Min(1.7), _m2Max(2.1),
      _ptMin(1.0), _ptMax(24.0), _reflOverSgn(0.2),
      _efficiencyMap(nullptr), _outputRootFile("")
{
    // Open input files
    _osFile = std::make_unique<TFile>(osFilename.c_str(), "READ");
    _lsFile = std::make_unique<TFile>(lsFilename.c_str(), "READ");
    
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
    
    if (!_osTree) {
        std::cerr << "Error: Cannot find tree 'pairs' in OS file" << std::endl;
    }
    if (!_lsTree) {
        std::cerr << "Error: Cannot find tree 'pairs' in LS file" << std::endl;
    }
    
    // Create workspace
    _workspace = std::make_unique<RooWorkspace>("w", "Invariant Mass Fit Workspace");
    
    std::cout << "InvariantMassFitter initialized" << std::endl;
    if (_osTree) std::cout << "  OS entries: " << _osTree->GetEntries() << std::endl;
    if (_lsTree) std::cout << "  LS entries: " << _lsTree->GetEntries() << std::endl;
}

InvariantMassFitter::~InvariantMassFitter() {
    // Smart pointers will handle cleanup
}

void InvariantMassFitter::setOutputFile(const std::string& filename) {
    _outputRootFile = filename;
    std::cout << "Output ROOT file set to: " << filename << std::endl;
}

void InvariantMassFitter::loadEfficiencyMap(const std::string& filename, const std::string& histName) {
    TFile* effFile = TFile::Open(filename.c_str(), "READ");
    if (!effFile || effFile->IsZombie()) {
        std::cerr << "Error: Cannot open efficiency file " << filename << std::endl;
        return;
    }
    
    TH2D* effHist = dynamic_cast<TH2D*>(effFile->Get(histName.c_str()));
    if (!effHist) {
        std::cerr << "Error: Cannot find histogram " << histName << " in efficiency file" << std::endl;
        effFile->Close();
        return;
    }
    
    // Clone the histogram so we own it
    _efficiencyMap = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(effHist->Clone("efficiencyMap")));
    _efficiencyMap->SetDirectory(nullptr); // Detach from file
    
    effFile->Close();
    delete effFile;
    
    std::cout << "Efficiency map loaded: " << histName << " from " << filename << std::endl;
    std::cout << "  Map dimensions: " << _efficiencyMap->GetNbinsX() << " x " 
              << _efficiencyMap->GetNbinsY() << " bins" << std::endl;
}

double InvariantMassFitter::calculateWeights(float const &pt, float const &y)
{
    if (!_efficiencyMap) return 1.0;
    
    const float epsilon = 1e-4;
    
    // Check if we're within the histogram bounds
    if (pt < _ptMin || 
        pt > _ptMax ||
        y < _efficiencyMap->GetYaxis()->GetXmin() || 
        y > _efficiencyMap->GetYaxis()->GetXmax()) {
        return 1.0; // Outside histogram range
    }

    // Directly find the bin for this (pt, y) point
    int bin = _efficiencyMap->FindBin(pt, y);

    // Extra safety: ensure bin is not overflow/underflow
    int ptBin = _efficiencyMap->GetXaxis()->FindBin(pt);
    int yBin = _efficiencyMap->GetYaxis()->FindBin(y);
    if (ptBin < 1 || ptBin > _efficiencyMap->GetXaxis()->GetNbins() ||
        yBin < 1 || yBin > _efficiencyMap->GetYaxis()->GetNbins()) {
        return 1.0;
    }
    
    double efficiency = _efficiencyMap->GetBinContent(bin);
    
    if (efficiency > epsilon) {
        return 1.0 / efficiency;
    }
    
    return 1.0; // Return 1 for invalid efficiency
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
    
    // Create kinematic variables (adjust ranges as needed)
    _y1 = std::make_unique<RooRealVar>("y1", "y_{1}", -2.0, 2.0);
    _y2 = std::make_unique<RooRealVar>("y2", "y_{2}", -2.0, 2.0);
    _pt1 = std::make_unique<RooRealVar>("pt1", "p_{T 1}", 0.0, 50.0);
    _pt2 = std::make_unique<RooRealVar>("pt2", "p_{T 2}", 0.0, 50.0);
    _pairPt = std::make_unique<RooRealVar>("pairPt", "p_{T pair}", 0.0, 100.0);
    _pairMass = std::make_unique<RooRealVar>("pairMass", "m_{pair}", 0.0, 100.0);
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
    if (!_m1 || !_m2) {
        createVariables();
    }
    
    bool useEfficiency = (_efficiencyMap != nullptr);
    
    // Create RooDataSet from OS tree
    if (_osTree && _osTree->GetEntries() > 0) {
        std::cout << "loading data part 1" << std::endl;
        RooArgSet varSet(*_m1, *_m2, *_y1, *_y2, *_phi1, *_phi2, *_pt1, *_pt2, *_pairPt, *_pairMass);
        RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 1000.);
        RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 1000.);
        
        // Create weighted or unweighted dataset
        if (useEfficiency) {
            
            varSet.add(weightCand1);
            varSet.add(weightCand2);
            RooFormulaVar combinedWeight("combinedWeight", "combined weight", "weightCand1 * weightCand2", RooArgList(weightCand1, weightCand2));
            _osData = std::make_unique<RooDataSet>("osData", "OS Data", varSet, RooFit::WeightVar(combinedWeight.GetName()));
        } else {
            _osData = std::make_unique<RooDataSet>("osData", "OS Data", varSet);
        }
        
        float m1_val, m2_val, mbar1_val, mbar2_val;
        float pt1, pt2, y1, y2, phi1, phi2;
        double pairPt, pairMass;
        UChar_t pairType;

        double totalWeight = 0.; // Variable to store the total weight for normalization
        
        _osTree->SetBranchAddress("m1", &m1_val);
        _osTree->SetBranchAddress("m2", &m2_val);
        _osTree->SetBranchAddress("mbar1", &mbar1_val);
        _osTree->SetBranchAddress("mbar2", &mbar2_val);
        _osTree->SetBranchAddress("pairType", &pairType);
        _osTree->SetBranchAddress("y1", &y1);
        _osTree->SetBranchAddress("y2", &y2);
        _osTree->SetBranchAddress("phi1", &phi1);
        _osTree->SetBranchAddress("phi2", &phi2);
        _osTree->SetBranchAddress("pt1", &pt1);
        _osTree->SetBranchAddress("pt2", &pt2);
        _osTree->SetBranchAddress("pairPt", &pairPt);
        _osTree->SetBranchAddress("pairMass", &pairMass);
        
        for (Long64_t i = 0; i < _osTree->GetEntries(); i++) {
            _osTree->GetEntry(i);
            
            double mass1, mass2;
            
            // DDbar (bit 2) or DbarD (bit 3)
            if (TESTBIT(pairType, DDbar)) { // DDbar
                mass1 = m1_val;
                mass2 = mbar2_val;
            } else if (TESTBIT(pairType, DbarD)) { // DbarD
                mass1 = mbar1_val;
                mass2 = m2_val;
            } else {
                continue; // Skip if not OS
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
                double wCand1  = calculateWeights(pt1, y1); // Weight = 1/efficiency
                varSet.setRealValue("weightCand1", wCand1);

                double wCand2  = calculateWeights(pt2, y2); // Weight = 1/efficiency
                varSet.setRealValue("weightCand2", wCand2);

                double combinedWeight = wCand1 * wCand2;

                _osData->add(varSet, combinedWeight);
            } else {
                _osData->add(varSet);
            }
        }
        
        std::cout << "Loaded " << _osData->numEntries() << " OS events";
        if (useEfficiency) {
            std::cout << " (efficiency-corrected, sum of weights = " 
                      << _osData->sumEntries() << ")";
        }
        std::cout << std::endl;
        _workspace->import(*_osData);
    }

    // Create RooDataSet from LS tree
    if (_lsTree && _lsTree->GetEntries() > 0) {
        RooArgSet varSet(*_m1, *_m2, *_y1, *_y2, *_phi1, *_phi2, *_pt1, *_pt2, *_pairPt, *_pairMass);
        RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 1000.);
        RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 1000.);
        
        if (useEfficiency) {

            varSet.add(weightCand1);
            varSet.add(weightCand2);
            RooFormulaVar combinedWeight("combinedWeight", "combined weight", "weightCand1 * weightCand2", RooArgList(weightCand1, weightCand2));
            _lsData = std::make_unique<RooDataSet>("lsData", "LS Data", varSet, RooFit::WeightVar(combinedWeight.GetName()));
        } else {
            _lsData = std::make_unique<RooDataSet>("lsData", "LS Data", varSet);
        }
        
        float m1_val, m2_val, mbar1_val, mbar2_val;
        float pt1, pt2, y1, y2, phi1, phi2;
        double pairPt, pairMass;
        UChar_t pairType;
        
        _lsTree->SetBranchAddress("m1", &m1_val);
        _lsTree->SetBranchAddress("m2", &m2_val);
        _lsTree->SetBranchAddress("mbar1", &mbar1_val);
        _lsTree->SetBranchAddress("mbar2", &mbar2_val);
        _lsTree->SetBranchAddress("pairType", &pairType);
        _lsTree->SetBranchAddress("y1", &y1);
        _lsTree->SetBranchAddress("y2", &y2);
        _lsTree->SetBranchAddress("phi1", &phi1);
        _lsTree->SetBranchAddress("phi2", &phi2);
        _lsTree->SetBranchAddress("pt1", &pt1);
        _lsTree->SetBranchAddress("pt2", &pt2);
        _lsTree->SetBranchAddress("pairPt", &pairPt);
        _lsTree->SetBranchAddress("pairMass", &pairMass);
        
        for (Long64_t i = 0; i < _lsTree->GetEntries(); i++) {
            _lsTree->GetEntry(i);
            
            double mass1, mass2;

            // DD (bit 0) or DbarDbar (bit 1)
            if (TESTBIT(pairType, DD)) { // DD
                mass1 = m1_val;
                mass2 = m2_val;
            } else if (TESTBIT(pairType, DbarDbar)) { // DbarDbar
                mass1 = mbar1_val;
                mass2 = mbar2_val;
            } else {
                continue; // Skip if not LS
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
                double wCand1  = calculateWeights(pt1, y1); // Weight = 1/efficiency
                varSet.setRealValue("weightCand1", wCand1);

                double wCand2  = calculateWeights(pt2, y2); // Weight = 1/efficiency
                varSet.setRealValue("weightCand2", wCand2);

                double combinedWeight = wCand1 * wCand2;

                _lsData->add(varSet, combinedWeight);
            } else {
                _lsData->add(varSet);
            }
        }
        
        std::cout << "Loaded " << _lsData->numEntries() << " LS events";
        if (useEfficiency) {
            std::cout << " (efficiency-corrected, sum of weights = " 
                      << _lsData->sumEntries() << ")";
        }
        std::cout << std::endl;
        _workspace->import(*_lsData);
    }
}

RooAbsPdf* InvariantMassFitter::createSignalPdf(const std::string& name, RooRealVar& mass) {
    if (_sgnFuncOption == "gaussian")
        return createGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "doubleGaussian")
        return createDoubleGaussianSignal(_workspace.get(), name, mass);
    else if (_sgnFuncOption == "crystalball")
        return createCrystalBallSignal(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown signal option: " + _sgnFuncOption);
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
        throw std::runtime_error("Unknown background option: " + _bkgFuncOption);
}

RooAbsPdf* InvariantMassFitter::createReflectedPdf(const std::string& name, RooRealVar& mass) {
    if (_reflFuncOption == "gaussian")
        return createGaussianReflected(_workspace.get(), name, mass);
    else if (_reflFuncOption == "doubleGaussian")
        return createDoubleGaussianReflected(_workspace.get(), name, mass);
    else
        throw std::runtime_error("Unknown reflected option: " + _reflFuncOption);
}

void InvariantMassFitter::buildModel(bool const doReflections) {

    if (!_m1 || !_m2) {
        createVariables();
    }
    
    // Create 1D PDFs for each mass dimension
    // m1 direction
    RooAbsPdf* sig1 = createSignalPdf("m1_sig", *_m1);
    RooAbsPdf* bkg1 = createBackgroundPdf("m1_bkg", *_m1);
    
    // m2 direction
    RooAbsPdf* sig2 = createSignalPdf("m2_sig", *_m2);
    RooAbsPdf* bkg2 = createBackgroundPdf("m2_bkg", *_m2);

    // Initial values
    int n_events = _lsData->numEntries();
    double initBkg = 0.7 * std::sqrt(n_events);
    double initSgn = 0.2 * std::sqrt(n_events);
    cout << "Initial parameters, bkg: " << initBkg << " sgn: " << initSgn << endl;

    // RAW PARAMETERS
    RooRealVar *nBkg1 = new RooRealVar("nBkg1", "background yield of cand 1", initBkg, 0.0, 1e6);
    RooRealVar *nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1", initSgn, 0.0, 1e5);
    RooFormulaVar *nRefl1 = new RooFormulaVar("nRefl1", "reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgn1));

    RooRealVar *nBkg2 = new RooRealVar("nBkg2", "background yield of cand2", initBkg, 0.0, 1e6);
    RooRealVar *nSgn2 = new RooRealVar("nSgn2", "signal yield of cand2", initSgn, 0.0, 1e5);
    RooFormulaVar *nRefl2 = new RooFormulaVar("nRefl2", "reflected signal yield of cand 2", "@0 * @1", RooArgList(_reflOverSgn, *nSgn2));

    RooFormulaVar nBkgBkg("nBkgBkg", "fraction of bkg x bkg", "@0 * @1", RooArgList(*nBkg1, *nBkg2));
    RooFormulaVar nBkgSgn("nBkgSgn", "fraction of bkg x sgn", "@0 * @1", RooArgList(*nBkg1, *nSgn2));
    RooFormulaVar nSgnBkg("nSgnBkg", "fraction of sgn x bkg", "@0 * @1", RooArgList(*nSgn1, *nBkg2));
    RooFormulaVar nBkgRefl("nBkgRefl", "fraction of bkg x refl", "@0 * @1", RooArgList(*nBkg1, *nRefl2));
    RooFormulaVar nReflBkg("nReflBkg", "fraction of refl x bkg", "@0 * @1", RooArgList(*nRefl1, *nBkg2));
    RooFormulaVar nSgnSgn("nSgnSgn", "fraction of sgn x sgn", "@0 * @1", RooArgList(*nSgn1, *nSgn2));
    RooFormulaVar nReflRefl("nReflRefl", "fraction of refl x refl", "@0 * @1", RooArgList(*nRefl1, *nRefl2));
    RooFormulaVar nSgnRefl("nSgnRefl", "fraction of sgn x refl", "@0 * @1", RooArgList(*nSgn1, *nRefl2));
    RooFormulaVar nReflSgn("nReflSgn", "fraction of refl x sgn", "@0 * @1", RooArgList(*nRefl1, *nSgn2));

    // Create 2D PDFs by taking products
    RooProdPdf pdfSigSig("pdfSigSig", "Signal-Signal", RooArgList(*sig1, *sig2));
    RooProdPdf pdfSigBkg("pdfSigBkg", "Signal-Background", RooArgList(*sig1, *bkg2));
    RooProdPdf pdfBkgSig("pdfBkgSig", "Background-Signal", RooArgList(*bkg1, *sig2));
    RooProdPdf pdfBkgBkg("pdfBkgBkg", "Background-Background", RooArgList(*bkg1, *bkg2));

    RooAbsPdf* refl1 = nullptr;
    RooAbsPdf* refl2 = nullptr;

    RooProdPdf* pdfSigRefl = nullptr;
    RooProdPdf* pdfReflSig = nullptr;
    RooProdPdf* pdfReflBkg = nullptr;
    RooProdPdf* pdfBkgRefl = nullptr;
    RooProdPdf* pdfReflRefl = nullptr;

    if (doReflections) {
        refl1 = createReflectedPdf("m1_refl", *_m1);
        refl2 = createReflectedPdf("m2_refl", *_m2);

        pdfSigRefl = new RooProdPdf("pdfSigRefl", "Signal-Reflected", RooArgList(*sig1, *refl2));
        pdfReflSig = new RooProdPdf("pdfReflSig", "Reflected-Signal", RooArgList(*refl1, *sig2));
        pdfReflBkg = new RooProdPdf("pdfReflBkg", "Reflected-Background", RooArgList(*refl1, *bkg2));
        pdfBkgRefl = new RooProdPdf("pdfBkgRefl", "Background-Reflected", RooArgList(*bkg1, *refl2));
        pdfReflRefl = new RooProdPdf("pdfReflRefl", "Reflected-Reflected", RooArgList(*refl1, *refl2));
    }

    // Create total PDF as sum of components
    RooAddPdf* totalPdf = nullptr;
    if (doReflections) {
        totalPdf = new RooAddPdf("model", "background + signal + reflection pdf 2D",
                                RooArgList(pdfSigSig, *pdfSigRefl, *pdfReflSig, pdfSigBkg,
                                           pdfBkgSig, pdfBkgBkg, *pdfBkgRefl, *pdfReflBkg, *pdfReflRefl),
                                RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));

    } else {
        totalPdf = new RooAddPdf("model", "background + signal 2D model",
                       RooArgList(pdfSigSig, pdfSigBkg, pdfBkgSig, pdfBkgBkg),
                       RooArgList(nSgnSgn, nSgnBkg, nBkgSgn, nBkgBkg));
    }

    
    _workspace->import(*totalPdf);
    
    std::cout << "2D model built" << std::endl;
}

RooFitResult* InvariantMassFitter::performPrefit1D(RooDataSet* data, RooRealVar& mass, 
                                                    const std::string& name, bool useReflections) {
    std::cout << "\n=== Performing 1D prefit for " << name << " ===" << std::endl;
    
    // Create 1D PDFs for this mass variable
    RooAbsPdf* sig = createSignalPdf(name + "_sig", mass);
    RooAbsPdf* bkg = createBackgroundPdf(name + "_bkg", mass);

    double nEntries = data->numEntries();
    
    // Create yield variables
    RooRealVar* nSig = new RooRealVar((name + "_nSig").c_str(), 
                                      "Signal yield", 
                                      0.5*nEntries, 0.0, 10*nEntries);
    RooRealVar* nBkg = new RooRealVar((name + "_nBkg").c_str(), 
                                      "Background yield", 
                                      0.2*nEntries, 0.0, 10*nEntries);
    
    RooAddPdf* model = nullptr;
    
    if (useReflections) {
        // Add reflection component
        RooAbsPdf* refl = createReflectedPdf(name + "_refl", mass);
        RooRealVar* reflOverSgn = new RooRealVar((name + "_reflOverSgn").c_str(),
                                                 "Reflection/Signal ratio",
                                                 _reflOverSgn, 0.0, 1.0);
        reflOverSgn->setConstant(true); // Keep ratio fixed in prefit
        
        RooFormulaVar* nRefl = new RooFormulaVar((name + "_nRefl").c_str(),
                                                 "Reflection yield",
                                                 "@0 * @1",
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
    
    // Perform the fit
    RooFitResult* result = nullptr;

    bool isWeighted = (_efficiencyMap != nullptr);
    
    if (isWeighted) {
        result = model->fitTo(*data, RooFit::Save(),
                      Extended(true),
                      AsymptoticError(true),  // Changed from SumW2Error
                      Strategy(2),
                      Optimize(2),
                      Offset(true),
                      Minimizer("Minuit2", "migrad"),
                      PrintLevel(1),
                      RooFit::MaxCalls(100000)
                    );
    } else {
        result = model->fitTo(*data, RooFit::Save(),
                            RooFit::Extended(true));
    }
    result->Print("v");
    
    // Print results
    std::cout << "Prefit results for " << name << ":" << std::endl;
    std::cout << "  Signal yield: " << nSig->getVal() << " ± " << nSig->getError() << std::endl;
    std::cout << "  Background yield: " << nBkg->getVal() << " ± " << nBkg->getError() << std::endl;
    
    // Store the model and result in workspace for later use
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
    
    // Perform prefits for OS data
    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n--- OS Data Prefits ---" << std::endl;
        
        // Prefit m1 projection
        _prefitResult_os_m1 = std::unique_ptr<RooFitResult>(
            performPrefit1D(_osData.get(), *_m1, "prefit_os_m1", useReflections)
        );
        
        // Prefit m2 projection
        _prefitResult_os_m2 = std::unique_ptr<RooFitResult>(
            performPrefit1D(_osData.get(), *_m2, "prefit_os_m2", useReflections)
        );
    }
    
    // Perform prefits for LS data
    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n--- LS Data Prefits ---" << std::endl;
        
        // Prefit m1 projection
        _prefitResult_ls_m1 = std::unique_ptr<RooFitResult>(
            performPrefit1D(_lsData.get(), *_m1, "prefit_ls_m1", useReflections)
        );
        
        // Prefit m2 projection
        _prefitResult_ls_m2 = std::unique_ptr<RooFitResult>(
            performPrefit1D(_lsData.get(), *_m2, "prefit_ls_m2", useReflections)
        );
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "=== 1D Prefit Procedure Complete ===" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void InvariantMassFitter::extractPrefitParameters(const std::string& prefitName, 
                                                   const std::string& targetName) {
    // This method automatically transfers all background parameters from prefit to 2D model
    // Works for any background function (exponential, polynomial, chebychev, etc.)
    
    std::cout << "  Extracting prefit parameters from " << prefitName << " to " << targetName << std::endl;
    
    // Get all parameters from the workspace
    RooArgSet allVars = _workspace->allVars();
    int nParamsSet = 0;
    for (RooAbsArg* arg : allVars) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
        if (!var) continue;
        std::string varName = var->GetName();
        
        // Check if this variable belongs to the prefit background
        if (varName.find(prefitName) != std::string::npos) {
            // Extract the parameter suffix (e.g., "_lambda", "_c0", "_c1", etc.)
            size_t pos = varName.find(prefitName);
            std::string suffix = varName.substr(pos + prefitName.length());
            
            // Skip the yield parameter
            if (suffix.empty() || varName.find("_nBkg") != std::string::npos || 
                varName.find("_nSig") != std::string::npos ||
                varName.find("_nRefl") != std::string::npos ||
                varName.find("_reflOverSgn") != std::string::npos) {
                continue;
            }
            
            // Construct the target parameter name
            std::string targetParamName = targetName + suffix;
            std::cout << targetParamName << std::endl;
            
            // Get the target parameter
            RooRealVar* targetParam = _workspace->var(targetParamName.c_str());
            
            if (targetParam) {
                // If the parameter name contains "_refl" → set constant
                if (varName.find("_refl") != std::string::npos) {
                    targetParam->setVal(var->getVal());
                    targetParam->setConstant(true);
                    std::cout << "    " << targetParamName 
                                << " set CONSTANT (_refl detected)" << std::endl;
                } 
                // Otherwise, set the value normally
                else {
                    targetParam->setVal(var->getVal());
                    std::cout << "    " << targetParamName 
                            << " set to: " << var->getVal() << std::endl;
                }
                nParamsSet++;
            }
        }
    }
    
    if (nParamsSet == 0) {
        std::cout << "    Warning: No background parameters found to transfer" << std::endl;
    } else {
        std::cout << "    Successfully set " << nParamsSet << " background parameter(s)" << std::endl;
    }
}


void InvariantMassFitter::setPrefitParameters(bool const doReflections) {
    std::cout << "\n=== Setting 2D fit initial parameters from prefit ===" << std::endl;
    
    // Helper lambda to get parameter value from workspace
    auto getParam = [this](const std::string& name) -> RooRealVar* {
        return _workspace->var(name.c_str());
    };
    
    // Set parameters for m1 based on OS prefit (you can choose OS or LS or average)
    if (_prefitResult_os_m1) {
        std::cout << "\n--- Setting m1 parameters from OS prefit ---" << std::endl;
        
        // Set all signal parameters automatically
        extractPrefitParameters("prefit_os_m1_sig", "m1_sig");
        
        // Set all background parameters automatically
        extractPrefitParameters("prefit_os_m1_bkg", "m1_bkg");
        
        // Set all reflection parameters if they exist
        if (doReflections) {
            extractPrefitParameters("prefit_os_m1_refl", "m1_refl");
        }
        
        // Set yield parameters
        RooRealVar* prefit_nSig = getParam("prefit_os_m1_nSig");
        RooRealVar* prefit_nBkg = getParam("prefit_os_m1_nBkg");
        
        if (prefit_nSig && prefit_nBkg) {
            RooRealVar* nSgn1 = getParam("nSgn1");
            RooRealVar* nBkg1 = getParam("nBkg1");
            
            if (nSgn1) {
                // Use sqrt of prefit yield as starting point for marginal yield
                double init_nSgn1 = std::sqrt(prefit_nSig->getVal());
                nSgn1->setVal(init_nSgn1);
                std::cout << "  nSgn1 set to: " << init_nSgn1 << std::endl;
            }
            
            if (nBkg1) {
                double init_nBkg1 = std::sqrt(prefit_nBkg->getVal());
                nBkg1->setVal(init_nBkg1);
                std::cout << "  nBkg1 set to: " << init_nBkg1 << std::endl;
            }
        }
    }
    
    // Set parameters for m2 based on OS prefit
    if (_prefitResult_os_m2) {
        std::cout << "\n--- Setting m2 parameters from OS prefit ---" << std::endl;
        
        // Set all signal parameters automatically
        extractPrefitParameters("prefit_os_m2_sig", "m2_sig");
        
        // Set all background parameters automatically
        extractPrefitParameters("prefit_os_m2_bkg", "m2_bkg");
        
        // Set all reflection parameters if they exist
        if (doReflections) {
            extractPrefitParameters("prefit_os_m2_refl", "m2_refl");
        }
        
        // Set yield parameters
        RooRealVar* prefit_nSig = getParam("prefit_os_m2_nSig");
        RooRealVar* prefit_nBkg = getParam("prefit_os_m2_nBkg");
        
        if (prefit_nSig && prefit_nBkg) {
            RooRealVar* nSgn2 = getParam("nSgn2");
            RooRealVar* nBkg2 = getParam("nBkg2");
            
            if (nSgn2) {
                double init_nSgn2 = std::sqrt(prefit_nSig->getVal());
                nSgn2->setVal(init_nSgn2);
                std::cout << "  nSgn2 set to: " << init_nSgn2 << std::endl;
            }
            
            if (nBkg2) {
                double init_nBkg2 = std::sqrt(prefit_nBkg->getVal());
                nBkg2->setVal(init_nBkg2);
                std::cout << "  nBkg2 set to: " << init_nBkg2 << std::endl;
            }
        }
    }
    
    std::cout << "\n=== Initial parameters set from prefit ===" << std::endl;
}

// Optional: Add a method to plot the prefit results
void InvariantMassFitter::plotPrefitResults(const std::string& outputPrefix) {
    std::cout << "\n=== Plotting prefit results ===" << std::endl;
    
    TCanvas* c = new TCanvas("c_prefit", "Prefit Results", 1600, 1200);
    c->Divide(2, 2);
    
    // Plot OS m1 prefit
    if (_prefitResult_os_m1 && _osData) {
        c->cd(1);
        RooPlot* frame = _m1->frame(RooFit::Title("OS m_{1} Prefit"));
        _osData->plotOn(frame, RooFit::DrawOption("pz"));
        
        RooAbsPdf* model = _workspace->pdf("prefit_os_m1_model");
        
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components("prefit_os_m1_sig_gauss"), 
                         DrawOption("F"), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2));
            model->plotOn(frame, RooFit::Components("prefit_os_m1_bkg_exp"), 
                         RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        }
        frame->Draw();
    }

    // Plot OS m2 prefit
    if (_prefitResult_os_m2 && _osData) {
        c->cd(2);
        RooPlot* frame = _m2->frame(RooFit::Title("OS m_{2} Prefit"));
        _osData->plotOn(frame, RooFit::DrawOption("pz"));
        
        RooAbsPdf* model = _workspace->pdf("prefit_os_m2_model");
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components("prefit_os_m2_sig_gauss"), 
                         DrawOption("F"), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2));
            model->plotOn(frame, RooFit::Components("prefit_os_m2_bkg_exp"), 
                         RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        }
        frame->Draw();
    }
    
    // Plot LS m1 prefit
    if (_prefitResult_ls_m1 && _lsData) {
        c->cd(3);
        RooPlot* frame = _m1->frame(RooFit::Title("LS m_{1} Prefit"));
        _lsData->plotOn(frame, RooFit::DrawOption("pz"));
        
        RooAbsPdf* model = _workspace->pdf("prefit_ls_m1_model");
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components("prefit_ls_m1_sig_gauss"), 
                         DrawOption("F"), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2));
            model->plotOn(frame, RooFit::Components("prefit_ls_m1_bkg_exp"), 
                         RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        }
        frame->Draw();
    }
    
    // Plot LS m2 prefit
    if (_prefitResult_ls_m2 && _lsData) {
        c->cd(4);
        RooPlot* frame = _m2->frame(RooFit::Title("LS m_{2} Prefit"));
        _lsData->plotOn(frame, RooFit::DrawOption("pz"));
        
        RooAbsPdf* model = _workspace->pdf("prefit_ls_m2_model");
        if (model) {
            model->plotOn(frame, RooFit::Name("prefit_model"));
            model->plotOn(frame, RooFit::Components("prefit_ls_m2_sig_gauss"), 
                         DrawOption("F"), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2));
            model->plotOn(frame, RooFit::Components("prefit_ls_m2_bkg_exp"), 
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
    
    RooAbsPdf* pdf = _workspace->pdf("model");
    
    // Determine if we're using weighted data
    bool isWeighted = (_efficiencyMap != nullptr);
    
    // Fit to OS data
    if ((dataType == "OS" || dataType == "both") && _osData) {
        std::cout << "\n=== Fitting OS data ===" << std::endl;
        RooFitResult* result = nullptr;
        if (isWeighted) {
            result = pdf->fitTo(*_osData, RooFit::Save(), 
                                Extended(true),
                                AsymptoticError(true),  // Changed from SumW2Error
                                Strategy(2),
                                Optimize(2),
                                Offset(true),
                                Minimizer("Minuit2", "migrad"),
                                PrintLevel(1),
                                RooFit::MaxCalls(100000)
                            );
        } else {
            result = pdf->fitTo(*_osData, RooFit::Save(), 
                              RooFit::Extended(true),
                              RooFit::PrintLevel(1));
        }
        result->Print("v");
        _workspace->import(*result, "osFitResult");

        // Save parameter snapshot
        RooArgSet* params = pdf->getParameters(*_osData);
        _workspace->saveSnapshot("osParams", *params); 

        delete result;
    }
    
    // Fit to LS data
    if ((dataType == "LS" || dataType == "both") && _lsData) {
        std::cout << "\n=== Fitting LS data ===" << std::endl;
        RooFitResult* result = nullptr;
        if (isWeighted) {
            result = pdf->fitTo(*_lsData, RooFit::Save(),
                              Extended(true),
                                AsymptoticError(true),  // Changed from SumW2Error
                                Strategy(2),
                                Optimize(2),
                                Offset(true),
                                Minimizer("Minuit2", "migrad"),
                                PrintLevel(1),
                                RooFit::MaxCalls(100000)
                            );
        } else {
            result = pdf->fitTo(*_lsData, RooFit::Save(),
                              RooFit::Extended(true),
                              RooFit::PrintLevel(1));
        }
        result->Print("v");
        _workspace->import(*result, "lsFitResult");

        // Save parameter snapshot
        RooArgSet* params = pdf->getParameters(*_lsData);
        _workspace->saveSnapshot("lsParams", *params); 

        delete result;
    }
}
void InvariantMassFitter::plotModelOnFrame(RooPlot* frame, RooAbsPdf* pdf, bool const doReflections) {
    pdf->plotOn(frame, Name("model_curve"));

    // Plot the components of the model and keep references to the returned RooCurve objects
    pdf->plotOn(frame, Components("pdfBkgBkg"), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfSigBkg"), LineStyle(3), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
    pdf->plotOn(frame, Components("pdfBkgSig"), LineStyle(4), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
    pdf->plotOn(frame, DrawOption("F"), Components("pdfSigSig"), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2), Name("sgnSgnFunc2D_curve"));

    if (doReflections) {
        pdf->plotOn(frame, Components("pdfBkgRefl"), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfSigRefl"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflBkg"), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflSig"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
        pdf->plotOn(frame, Components("pdfReflRefl"), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
    }

    frame->Draw();
}

double InvariantMassFitter::computeChi2(RooPlot* frame, RooDataSet* data, int nFitParams) {
    // frame->chiSquare() computes chi2/ndf where:
    //   - numerator:   sum over bins of (data - model)^2 / error^2
    //   - denominator: nBins - nFitParams
    //
    // For WEIGHTED data, RooFit uses sumW2 errors on each bin (i.e. sqrt(sum of w^2)),
    // which is the correct Poisson-like error on the weighted histogram.
    // We pass the number of unweighted entries as N_eff so the ndf is meaningful.
    //
    // For UNWEIGHTED data, this reduces to the standard chi2/ndf.

    bool isWeighted = (_efficiencyMap != nullptr);

    double chi2Ndf = frame->chiSquare("model_curve", "data_cand1", nFitParams);

    int nBins = frame->GetXaxis()->GetNbins();
    int ndf   = nBins - nFitParams;
    double chi2 = chi2Ndf * ndf;

    std::cout << "  Chi2 / ndf  = " << chi2 << " / " << ndf
              << " = " << chi2Ndf << std::endl;

    if (isWeighted) {
        // For weighted fits, also report the effective number of entries
        // (unweighted count), which sets the statistical power of the test
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
    
    RooAbsPdf* pdf = _workspace->pdf("model");
    
    // Create canvas
    TCanvas* c = new TCanvas("c_fit", "Fit Results", 1600, 800);
    c->Divide(2, 2);

    double binWidth = 0.005; // Desired bin width
    int nBins = (int)((_m1Max - _m1Min) / binWidth);

    // Set the calculated number of bins
    _m1->setBins(nBins);
    _m2->setBins(nBins);
    
    // Count only floating (non-constant) parameters for chi2/ndf
    RooArgSet* allParams = pdf->getParameters(*_osData);
    RooArgSet floatingParams;

    for (RooAbsArg* arg : *allParams) {
        RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
        if (var && !var->isConstant()) floatingParams.add(*var);
    }

    delete allParams;
    int nFitParams = floatingParams.getSize();
    std::cout << "  Floating fit parameters: " << nFitParams << std::endl;

    // Plot OS data
    if (_osData) {
        _workspace->loadSnapshot("osParams");

        // m1 projection
        c->cd(1);
        RooPlot* frame1 = _m1->frame(RooFit::Title("OS: m_{1} projection"));
        _osData->plotOn(frame1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame1, pdf, doReflections);
        std::cout << "OS m1: "; computeChi2(frame1, _osData.get(), nFitParams);

        // m2 projection
        c->cd(2);
        RooPlot* frame2 = _m2->frame(RooFit::Title("OS: m_{2} projection"));
        _osData->plotOn(frame2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame2, pdf, doReflections);
        std::cout << "OS m2: "; computeChi2(frame2, _osData.get(), nFitParams);
    }
    
    // Plot LS data
    if (_lsData) {
        _workspace->loadSnapshot("lsParams");

        // m1 projection
        c->cd(3);
        RooPlot* frame3 = _m1->frame(RooFit::Title("LS: m_{1} projection"));
        _lsData->plotOn(frame3, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame3, pdf, doReflections);
        std::cout << "LS m1: "; computeChi2(frame3, _lsData.get(), nFitParams);
        
        // m2 projection
        c->cd(4);
        RooPlot* frame4 = _m2->frame(RooFit::Title("LS: m_{2} projection"));
        _lsData->plotOn(frame4, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame4, pdf, doReflections);
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
    
    // Save the workspace (contains all PDFs and parameters)
    _workspace->Write();
    
    // Create and save 1D projection histograms for easy plotting in Python
    RooAbsPdf* pdf = _workspace->pdf("model");
    
    if (_osData && pdf) {
        // Load OS parameters
        _workspace->loadSnapshot("osParams");
        
        // OS m1 projection
        RooPlot* frame_os_m1 = _m1->frame(RooFit::Title("OS: m_{1} projection"));
        _osData->plotOn(frame_os_m1, RooFit::Name("data_os_m1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_os_m1, pdf, doReflections);

        // Convert RooPlot to TH1 - this is just the data histogram
        TH1* h_data_os_m1 = _osData->createHistogram("h_data_os_m1", *_m1);
        h_data_os_m1->SetTitle("OS Data m_{1} projection;m_{1} [GeV/c^{2}];Entries");
        h_data_os_m1->Write();
        
        // Save the RooPlot frame (contains data + fit curves)
        frame_os_m1->Write("frame_os_m1");
        
        // OS m2 projection
        RooPlot* frame_os_m2 = _m2->frame(RooFit::Title("OS: m_{2} projection"));
        _osData->plotOn(frame_os_m2, RooFit::Name("data_os_m2"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_os_m2, pdf, doReflections);
        
        TH1* h_data_os_m2 = _osData->createHistogram("h_data_os_m2", *_m2);
        h_data_os_m2->SetTitle("OS Data m_{2} projection;m_{2} [GeV/c^{2}];Entries");
        h_data_os_m2->Write();
        
        // Save the RooPlot frame
        frame_os_m2->Write("frame_os_m2");
    }
    
    if (_lsData && pdf) {
        // Load LS parameters
        _workspace->loadSnapshot("lsParams");
        
        // LS m1 projection
        RooPlot* frame_ls_m1 = _m1->frame(RooFit::Title("LS: m_{1} projection"));
        _lsData->plotOn(frame_ls_m1, RooFit::Name("data_ls_m1"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_ls_m1, pdf, doReflections);
        
        TH1* h_data_ls_m1 = _lsData->createHistogram("h_data_ls_m1", *_m1);
        h_data_ls_m1->SetTitle("LS Data m_{1} projection;m_{1} [GeV/c^{2}];Entries");
        h_data_ls_m1->Write();
        
        // Save the RooPlot frame
        frame_ls_m1->Write("frame_ls_m1");
        
        // LS m2 projection
        RooPlot* frame_ls_m2 = _m2->frame(RooFit::Title("LS: m_{2} projection"));
        _lsData->plotOn(frame_ls_m2, RooFit::Name("data_ls_m2"), RooFit::DrawOption("pz"));
        plotModelOnFrame(frame_ls_m2, pdf, doReflections);
        
        TH1* h_data_ls_m2 = _lsData->createHistogram("h_data_ls_m2", *_m2);
        h_data_ls_m2->SetTitle("LS Data m_{2} projection;m_{2} [GeV/c^{2}];Entries");
        h_data_ls_m2->Write();
        
        // Save the RooPlot frame
        frame_ls_m2->Write("frame_ls_m2");
    }
    
    // Create and save 2D histograms
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
    
    // Save fit results
    RooFitResult* osRes = nullptr;
    RooFitResult* lsRes = nullptr;
    
    if (_workspace->obj("osFitResult")) {
        osRes = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
        if (osRes) {
            osRes->Write("osFitResult");
        }
    }
    
    if (_workspace->obj("lsFitResult")) {
        lsRes = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));
        if (lsRes) {
            lsRes->Write("lsFitResult");
        }
    }
    
    // Get signal PDF parameters for 3-sigma range
    RooRealVar* m1_mean = _workspace->var("m1_sig_mean");
    RooRealVar* m1_sigma = _workspace->var("m1_sig_sigma");
    RooRealVar* m2_mean = _workspace->var("m2_sig_mean");
    RooRealVar* m2_sigma = _workspace->var("m2_sig_sigma");
    
    if (!m1_mean || !m1_sigma || !m2_mean || !m2_sigma) {
        std::cerr << "Error: Cannot find signal parameters in workspace" << std::endl;
        outFile->Close();
        delete outFile;
        return;
    }
    
    // Set the 3-sigma range for m1 and m2
    _m1->setRange("3sigmaRange", m1_mean->getVal() - 3 * m1_sigma->getVal(),
                                 m1_mean->getVal() + 3 * m1_sigma->getVal());
    _m2->setRange("3sigmaRange", m2_mean->getVal() - 3 * m2_sigma->getVal(),
                                 m2_mean->getVal() + 3 * m2_sigma->getVal());
    
    // Get the 2D PDF components
    RooProdPdf* pdfSigSig = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfSigSig"));
    RooProdPdf* pdfSigBkg = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfSigBkg"));
    RooProdPdf* pdfBkgSig = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfBkgSig"));
    RooProdPdf* pdfBkgBkg = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfBkgBkg"));
    
    // Get yield variables
    RooFormulaVar* nSgnSgn = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnSgn"));
    RooFormulaVar* nSgnBkg = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnBkg"));
    RooFormulaVar* nBkgSgn = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgSgn"));
    RooFormulaVar* nBkgBkg = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgBkg"));
    
    // Initialize reflection variables
    RooProdPdf* pdfSigRefl = nullptr;
    RooProdPdf* pdfReflSig = nullptr;
    RooProdPdf* pdfBkgRefl = nullptr;
    RooProdPdf* pdfReflBkg = nullptr;
    RooProdPdf* pdfReflRefl = nullptr;
    RooFormulaVar* nSgnRefl = nullptr;
    RooFormulaVar* nReflSgn = nullptr;
    RooFormulaVar* nBkgRefl = nullptr;
    RooFormulaVar* nReflBkg = nullptr;
    RooFormulaVar* nReflRefl = nullptr;
    
    if (doReflections) {
        pdfSigRefl = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfSigRefl"));
        pdfReflSig = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfReflSig"));
        pdfBkgRefl = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfBkgRefl"));
        pdfReflBkg = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfReflBkg"));
        pdfReflRefl = dynamic_cast<RooProdPdf*>(_workspace->pdf("pdfReflRefl"));
        
        nSgnRefl = dynamic_cast<RooFormulaVar*>(_workspace->function("nSgnRefl"));
        nReflSgn = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflSgn"));
        nBkgRefl = dynamic_cast<RooFormulaVar*>(_workspace->function("nBkgRefl"));
        nReflBkg = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflBkg"));
        nReflRefl = dynamic_cast<RooFormulaVar*>(_workspace->function("nReflRefl"));
    }
    
    // Calculate OS yields if fit result exists
    RooRealVar* yieldSgnSgn_os = nullptr;
    RooRealVar* yieldSgnBkg_os = nullptr;
    RooRealVar* yieldBkgSgn_os = nullptr;
    RooRealVar* yieldBkgBkg_os = nullptr;
    RooRealVar* yieldSgnRefl_os = nullptr;
    RooRealVar* yieldReflSgn_os = nullptr;
    RooRealVar* yieldBkgRefl_os = nullptr;
    RooRealVar* yieldReflBkg_os = nullptr;
    RooRealVar* yieldReflRefl_os = nullptr;
    
    if (osRes && pdfSigSig && nSgnSgn) {
        std::cout << "\n=== OS Yields in 3-sigma range ===" << std::endl;
        _workspace->loadSnapshot("osParams");
        
        yieldSgnSgn_os = getYieldInRange(osRes, pdfSigSig, nSgnSgn, "3sigmaRange");
        yieldSgnBkg_os = getYieldInRange(osRes, pdfSigBkg, nSgnBkg, "3sigmaRange");
        yieldBkgSgn_os = getYieldInRange(osRes, pdfBkgSig, nBkgSgn, "3sigmaRange");
        yieldBkgBkg_os = getYieldInRange(osRes, pdfBkgBkg, nBkgBkg, "3sigmaRange");
        
        std::cout << " SgnSgn yield: " << yieldSgnSgn_os->getVal() << " +- " << yieldSgnSgn_os->getError() << std::endl;
        std::cout << " SgnBkg yield: " << yieldSgnBkg_os->getVal() << " +- " << yieldSgnBkg_os->getError() << std::endl;
        std::cout << " BkgSgn yield: " << yieldBkgSgn_os->getVal() << " +- " << yieldBkgSgn_os->getError() << std::endl;
        std::cout << " BkgBkg yield: " << yieldBkgBkg_os->getVal() << " +- " << yieldBkgBkg_os->getError() << std::endl;
        
        if (doReflections && pdfSigRefl) {
            yieldSgnRefl_os = getYieldInRange(osRes, pdfSigRefl, nSgnRefl, "3sigmaRange");
            yieldReflSgn_os = getYieldInRange(osRes, pdfReflSig, nReflSgn, "3sigmaRange");
            yieldBkgRefl_os = getYieldInRange(osRes, pdfBkgRefl, nBkgRefl, "3sigmaRange");
            yieldReflBkg_os = getYieldInRange(osRes, pdfReflBkg, nReflBkg, "3sigmaRange");
            yieldReflRefl_os = getYieldInRange(osRes, pdfReflRefl, nReflRefl, "3sigmaRange");
            
            std::cout << " SgnRefl yield: " << yieldSgnRefl_os->getVal() << " +- " << yieldSgnRefl_os->getError() << std::endl;
            std::cout << " ReflSgn yield: " << yieldReflSgn_os->getVal() << " +- " << yieldReflSgn_os->getError() << std::endl;
            std::cout << " BkgRefl yield: " << yieldBkgRefl_os->getVal() << " +- " << yieldBkgRefl_os->getError() << std::endl;
            std::cout << " ReflBkg yield: " << yieldReflBkg_os->getVal() << " +- " << yieldReflBkg_os->getError() << std::endl;
            std::cout << " ReflRefl yield: " << yieldReflRefl_os->getVal() << " +- " << yieldReflRefl_os->getError() << std::endl;
        }
    }
    
    // Calculate LS yields if fit result exists
    RooRealVar* yieldSgnSgn_ls = nullptr;
    RooRealVar* yieldSgnBkg_ls = nullptr;
    RooRealVar* yieldBkgSgn_ls = nullptr;
    RooRealVar* yieldBkgBkg_ls = nullptr;
    RooRealVar* yieldSgnRefl_ls = nullptr;
    RooRealVar* yieldReflSgn_ls = nullptr;
    RooRealVar* yieldBkgRefl_ls = nullptr;
    RooRealVar* yieldReflBkg_ls = nullptr;
    RooRealVar* yieldReflRefl_ls = nullptr;
    
    if (lsRes && pdfSigSig && nSgnSgn) {
        std::cout << "\n=== LS Yields in 3-sigma range ===" << std::endl;
        _workspace->loadSnapshot("lsParams");
        
        yieldSgnSgn_ls = getYieldInRange(lsRes, pdfSigSig, nSgnSgn, "3sigmaRange");
        yieldSgnBkg_ls = getYieldInRange(lsRes, pdfSigBkg, nSgnBkg, "3sigmaRange");
        yieldBkgSgn_ls = getYieldInRange(lsRes, pdfBkgSig, nBkgSgn, "3sigmaRange");
        yieldBkgBkg_ls = getYieldInRange(lsRes, pdfBkgBkg, nBkgBkg, "3sigmaRange");
        
        std::cout << " SgnSgn yield: " << yieldSgnSgn_ls->getVal() << " +- " << yieldSgnSgn_ls->getError() << std::endl;
        std::cout << " SgnBkg yield: " << yieldSgnBkg_ls->getVal() << " +- " << yieldSgnBkg_ls->getError() << std::endl;
        std::cout << " BkgSgn yield: " << yieldBkgSgn_ls->getVal() << " +- " << yieldBkgSgn_ls->getError() << std::endl;
        std::cout << " BkgBkg yield: " << yieldBkgBkg_ls->getVal() << " +- " << yieldBkgBkg_ls->getError() << std::endl;
        
        if (doReflections && pdfSigRefl) {
            yieldSgnRefl_ls = getYieldInRange(lsRes, pdfSigRefl, nSgnRefl, "3sigmaRange");
            yieldReflSgn_ls = getYieldInRange(lsRes, pdfReflSig, nReflSgn, "3sigmaRange");
            yieldBkgRefl_ls = getYieldInRange(lsRes, pdfBkgRefl, nBkgRefl, "3sigmaRange");
            yieldReflBkg_ls = getYieldInRange(lsRes, pdfReflBkg, nReflBkg, "3sigmaRange");
            yieldReflRefl_ls = getYieldInRange(lsRes, pdfReflRefl, nReflRefl, "3sigmaRange");
            
            std::cout << " SgnRefl yield: " << yieldSgnRefl_ls->getVal() << " +- " << yieldSgnRefl_ls->getError() << std::endl;
            std::cout << " ReflSgn yield: " << yieldReflSgn_ls->getVal() << " +- " << yieldReflSgn_ls->getError() << std::endl;
            std::cout << " BkgRefl yield: " << yieldBkgRefl_ls->getVal() << " +- " << yieldBkgRefl_ls->getError() << std::endl;
            std::cout << " ReflBkg yield: " << yieldReflBkg_ls->getVal() << " +- " << yieldReflBkg_ls->getError() << std::endl;
            std::cout << " ReflRefl yield: " << yieldReflRefl_ls->getVal() << " +- " << yieldReflRefl_ls->getError() << std::endl;
        }
    }
    
    std::cout << "\n ---------------------------------------------- " << std::endl;
    std::cout << " ---------------------------------------------- \n" << std::endl;

    // --- Save fit parameters (mean, sigma, chi2) for OS and LS ---
    // Stored as labeled TH1D so they can be accessed in Python like:
    //   data[pt_bin]["OS"]["hFitParams"]["values"][i]
    // Bins: m1_mean, m1_mean_err, m1_sigma, m1_sigma_err,
    //       m2_mean, m2_mean_err, m2_sigma, m2_sigma_err,
    //       chi2_m1, chi2_m2

    auto makeFitParamsHisto = [&](const std::string& name,
                                  RooDataSet* data,
                                  const std::string& snapshotName) -> bool {

        _workspace->loadSnapshot(snapshotName.c_str());

        RooRealVar* m1_mean_snap  = _workspace->var("m1_sig_mean");
        RooRealVar* m1_sigma_snap = _workspace->var("m1_sig_sigma");
        RooRealVar* m2_mean_snap  = _workspace->var("m2_sig_mean");
        RooRealVar* m2_sigma_snap = _workspace->var("m2_sig_sigma");

        if (!m1_mean_snap || !m1_sigma_snap || !m2_mean_snap || !m2_sigma_snap) {
            std::cerr << "Warning: signal parameters not found for " << name << std::endl;
            return false;
        }

        // Compute chi2 for m1 and m2 projections using the saved frames
        // Re-build temporary frames to get chi2
        int nBinsFrame = (int)((_m1Max - _m1Min) / 0.005);
        _m1->setBins(nBinsFrame);
        _m2->setBins(nBinsFrame);

        RooAbsPdf* pdfModel = _workspace->pdf("model");

        RooArgSet* allParams = pdfModel->getParameters(*data);
        int nFloating = 0;

        for (RooAbsArg* arg : *allParams) {
            RooRealVar* v = dynamic_cast<RooRealVar*>(arg);
            if (v && !v->isConstant()) nFloating++;
        }
        delete allParams;

        RooPlot* f1 = _m1->frame();
        data->plotOn(f1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f1, RooFit::Name("model_curve"));
        double chi2_m1 = f1->chiSquare("model_curve", "data_cand1", nFloating);

        RooPlot* f2 = _m2->frame();
        data->plotOn(f2, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
        pdfModel->plotOn(f2, RooFit::Name("model_curve"));
        double chi2_m2 = f2->chiSquare("model_curve", "data_cand1", nFloating);

        // Helper to make a single-bin TH1D with value and error
        auto makeHisto = [](const std::string& hname, const std::string& title,
                            double val, double err) -> TH1D* {
            TH1D* h = new TH1D(hname.c_str(), title.c_str(), 1, 0, 1);
            h->SetBinContent(1, val);
            h->SetBinError(1, err);
            return h;
        };

        // Save individual histograms
        std::vector<TH1D*> histos = {
            makeHisto(name + "_m1Mean",    "m1 mean [GeV/c^{2}]",    m1_mean_snap->getVal(),  m1_mean_snap->getError()),
            makeHisto(name + "_m1Sigma",   "m1 sigma [GeV/c^{2}]",   m1_sigma_snap->getVal(), m1_sigma_snap->getError()),
            makeHisto(name + "_m2Mean",    "m2 mean [GeV/c^{2}]",    m2_mean_snap->getVal(),  m2_mean_snap->getError()),
            makeHisto(name + "_m2Sigma",   "m2 sigma [GeV/c^{2}]",   m2_sigma_snap->getVal(), m2_sigma_snap->getError()),
            makeHisto(name + "_chi2ndf_m1","chi2/ndf m1 projection",  chi2_m1,                 0.0),
            makeHisto(name + "_chi2ndf_m2","chi2/ndf m2 projection",  chi2_m2,                 0.0),
        };
        for (TH1D* h : histos) { h->Write(); delete h; }

        delete f1;
        delete f2;
        return true;
    };

    if (_osData && osRes) makeFitParamsHisto("hOS", _osData.get(), "osParams");
    if (_lsData && lsRes) makeFitParamsHisto("hLS", _lsData.get(), "lsParams");

    // Save OS yields in histogram
    if (yieldSgnSgn_os) {
        int nBins = doReflections ? 9 : 4;
        TH1D* hYields_os = new TH1D("hYields_os", "OS Yields in 3#sigma range", nBins, 0, nBins);
        hYields_os->GetXaxis()->SetBinLabel(1, "SgnSgn");
        hYields_os->GetXaxis()->SetBinLabel(2, "SgnBkg");
        hYields_os->GetXaxis()->SetBinLabel(3, "BkgSgn");
        hYields_os->GetXaxis()->SetBinLabel(4, "BkgBkg");
        
        hYields_os->SetBinContent(1, yieldSgnSgn_os->getVal());
        hYields_os->SetBinError(1, yieldSgnSgn_os->getError());
        hYields_os->SetBinContent(2, yieldSgnBkg_os->getVal());
        hYields_os->SetBinError(2, yieldSgnBkg_os->getError());
        hYields_os->SetBinContent(3, yieldBkgSgn_os->getVal());
        hYields_os->SetBinError(3, yieldBkgSgn_os->getError());
        hYields_os->SetBinContent(4, yieldBkgBkg_os->getVal());
        hYields_os->SetBinError(4, yieldBkgBkg_os->getError());
        
        if (doReflections && yieldSgnRefl_os) {
            hYields_os->GetXaxis()->SetBinLabel(5, "SgnRefl");
            hYields_os->GetXaxis()->SetBinLabel(6, "ReflSgn");
            hYields_os->GetXaxis()->SetBinLabel(7, "BkgRefl");
            hYields_os->GetXaxis()->SetBinLabel(8, "ReflBkg");
            hYields_os->GetXaxis()->SetBinLabel(9, "ReflRefl");
            
            hYields_os->SetBinContent(5, yieldSgnRefl_os->getVal());
            hYields_os->SetBinError(5, yieldSgnRefl_os->getError());
            hYields_os->SetBinContent(6, yieldReflSgn_os->getVal());
            hYields_os->SetBinError(6, yieldReflSgn_os->getError());
            hYields_os->SetBinContent(7, yieldBkgRefl_os->getVal());
            hYields_os->SetBinError(7, yieldBkgRefl_os->getError());
            hYields_os->SetBinContent(8, yieldReflBkg_os->getVal());
            hYields_os->SetBinError(8, yieldReflBkg_os->getError());
            hYields_os->SetBinContent(9, yieldReflRefl_os->getVal());
            hYields_os->SetBinError(9, yieldReflRefl_os->getError());
        }
        
        hYields_os->Write();
    }
    
    // Save LS yields in histogram
    if (yieldSgnSgn_ls) {
        int nBins = doReflections ? 9 : 4;
        TH1D* hYields_ls = new TH1D("hYields_ls", "LS Yields in 3#sigma range", nBins, 0, nBins);
        hYields_ls->GetXaxis()->SetBinLabel(1, "SgnSgn");
        hYields_ls->GetXaxis()->SetBinLabel(2, "SgnBkg");
        hYields_ls->GetXaxis()->SetBinLabel(3, "BkgSgn");
        hYields_ls->GetXaxis()->SetBinLabel(4, "BkgBkg");
        
        hYields_ls->SetBinContent(1, yieldSgnSgn_ls->getVal());
        hYields_ls->SetBinError(1, yieldSgnSgn_ls->getError());
        hYields_ls->SetBinContent(2, yieldSgnBkg_ls->getVal());
        hYields_ls->SetBinError(2, yieldSgnBkg_ls->getError());
        hYields_ls->SetBinContent(3, yieldBkgSgn_ls->getVal());
        hYields_ls->SetBinError(3, yieldBkgSgn_ls->getError());
        hYields_ls->SetBinContent(4, yieldBkgBkg_ls->getVal());
        hYields_ls->SetBinError(4, yieldBkgBkg_ls->getError());
        
        if (doReflections && yieldSgnRefl_ls) {
            hYields_ls->GetXaxis()->SetBinLabel(5, "SgnRefl");
            hYields_ls->GetXaxis()->SetBinLabel(6, "ReflSgn");
            hYields_ls->GetXaxis()->SetBinLabel(7, "BkgRefl");
            hYields_ls->GetXaxis()->SetBinLabel(8, "ReflBkg");
            hYields_ls->GetXaxis()->SetBinLabel(9, "ReflRefl");
            
            hYields_ls->SetBinContent(5, yieldSgnRefl_ls->getVal());
            hYields_ls->SetBinError(5, yieldSgnRefl_ls->getError());
            hYields_ls->SetBinContent(6, yieldReflSgn_ls->getVal());
            hYields_ls->SetBinError(6, yieldReflSgn_ls->getError());
            hYields_ls->SetBinContent(7, yieldBkgRefl_ls->getVal());
            hYields_ls->SetBinError(7, yieldBkgRefl_ls->getError());
            hYields_ls->SetBinContent(8, yieldReflBkg_ls->getVal());
            hYields_ls->SetBinError(8, yieldReflBkg_ls->getError());
            hYields_ls->SetBinContent(9, yieldReflRefl_ls->getVal());
            hYields_ls->SetBinError(9, yieldReflRefl_ls->getError());
        }
        
        hYields_ls->Write();
    }
    
    // Save efficiency map if it exists
    if (_efficiencyMap) {
        _efficiencyMap->Write("efficiencyMap");
    }
    
    // Clean up dynamically allocated RooRealVar objects
    delete yieldSgnSgn_os;
    delete yieldSgnBkg_os;
    delete yieldBkgSgn_os;
    delete yieldBkgBkg_os;
    delete yieldSgnRefl_os;
    delete yieldReflSgn_os;
    delete yieldBkgRefl_os;
    delete yieldReflBkg_os;
    delete yieldReflRefl_os;
    
    delete yieldSgnSgn_ls;
    delete yieldSgnBkg_ls;
    delete yieldBkgSgn_ls;
    delete yieldBkgBkg_ls;
    delete yieldSgnRefl_ls;
    delete yieldReflSgn_ls;
    delete yieldBkgRefl_ls;
    delete yieldReflBkg_ls;
    delete yieldReflRefl_ls;
    
    outFile->Close();
    delete outFile;
    
    std::cout << "Results saved successfully to " << _outputRootFile << std::endl;
    std::cout << "File contains:" << std::endl;
    std::cout << "  - Workspace 'w' with all PDFs and parameters" << std::endl;
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
    
    if (_workspace->obj("osFitResult")) {
        std::cout << "\n--- OS Fit ---" << std::endl;
        RooFitResult* osRes = dynamic_cast<RooFitResult*>(_workspace->obj("osFitResult"));
        if (osRes) osRes->Print("v");
    }
    
    if (_workspace->obj("lsFitResult")) {
        std::cout << "\n--- LS Fit ---" << std::endl;
        RooFitResult* lsRes = dynamic_cast<RooFitResult*>(_workspace->obj("lsFitResult"));
        if (lsRes) lsRes->Print("v");
    }
    
    // Print yield parameters
    RooRealVar* nSigSig = _workspace->var("nSigSig");
    RooRealVar* nSigBkg = _workspace->var("nSigBkg");
    RooRealVar* nBkgSig = _workspace->var("nBkgSig");
    RooRealVar* nBkgBkg = _workspace->var("nBkgBkg");
    
    if (nSigSig && nSigBkg && nBkgSig && nBkgBkg) {
        std::cout << "\n=== Yield Summary ===" << std::endl;
        std::cout << "Signal-Signal:      " << nSigSig->getVal() << " ± " << nSigSig->getError() << std::endl;
        std::cout << "Signal-Background:  " << nSigBkg->getVal() << " ± " << nSigBkg->getError() << std::endl;
        std::cout << "Background-Signal:  " << nBkgSig->getVal() << " ± " << nBkgSig->getError() << std::endl;
        std::cout << "Background-Background: " << nBkgBkg->getVal() << " ± " << nBkgBkg->getError() << std::endl;
        
        if (_efficiencyMap) {
            std::cout << "\n(Yields are efficiency-corrected)" << std::endl;
        }
    }
}

RooRealVar *InvariantMassFitter::getYieldInRange(RooFitResult *fitResult, RooProdPdf* function, RooFormulaVar* nCands, TString range)
{
  // Calculate the integrals over the specified range and the full range
  RooAbsReal *intRange = function->createIntegral(RooArgSet(*_m1, *_m2), Range(range));
  RooAbsReal *intTotal = function->createIntegral(RooArgSet(*_m1, *_m2));

  double yieldRange = nCands->getVal() * intRange->getVal() / intTotal->getVal();

  // Propagate the uncertainties
  double nCandsErr = nCands->getPropagatedError(*fitResult);

  // Calculate the uncertainty using error propagation
  double errorRange = nCandsErr * intRange->getVal() / intTotal->getVal();

  // Initialize the RooRealVar to hold the yield and error
  RooRealVar *yield = new RooRealVar("yield", "Yield in 3-sigma range", yieldRange);
  yield->setError(errorRange);
  return yield;
}
