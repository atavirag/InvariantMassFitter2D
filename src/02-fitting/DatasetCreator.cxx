// DatasetCreator.cxx

#include "DatasetCreator.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"
#include "RooFit.h"
#include <iostream>
#include <functional>

#include "../01-prefiltering/CandidatePair.h"  // TESTBIT, pairType enum

using namespace RooFit;

DatasetCreator::DatasetCreator(const std::string& osFilename,
                               const std::string& lsFilename)
{
    _osFile = std::make_unique<TFile>(osFilename.c_str(), "READ");
    _lsFile = std::make_unique<TFile>(lsFilename.c_str(), "READ");

    if (!_osFile || _osFile->IsZombie())
        std::cerr << "[DatasetCreator] Cannot open OS file: " << osFilename << std::endl;
    else
        _osTree = dynamic_cast<TTree*>(_osFile->Get("pairs"));

    if (!_lsFile || _lsFile->IsZombie())
        std::cerr << "[DatasetCreator] Cannot open LS file: " << lsFilename << std::endl;
    else
        _lsTree = dynamic_cast<TTree*>(_lsFile->Get("pairs"));

    if (!_osTree) std::cerr << "[DatasetCreator] No 'pairs' tree in OS file\n";
    if (!_lsTree) std::cerr << "[DatasetCreator] No 'pairs' tree in LS file\n";

    std::cout << "[DatasetCreator] Initialised";
    if (_osTree) std::cout << "  OS=" << _osTree->GetEntries();
    if (_lsTree) std::cout << "  LS=" << _lsTree->GetEntries();
    std::cout << std::endl;
}

void DatasetCreator::setMassRange(double m1Min, double m1Max,
                                   double m2Min, double m2Max)
{
    _m1Min = m1Min; _m1Max = m1Max;
    _m2Min = m2Min; _m2Max = m2Max;
}

void DatasetCreator::setPtRange(double ptMin, double ptMax)
{
    _ptMin = ptMin; _ptMax = ptMax;
}

void DatasetCreator::loadEfficiencyMap(const std::string& filename,
                                        const std::string& histName)
{
    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[DatasetCreator] Cannot open efficiency file: " << filename << std::endl;
        return;
    }
    TH2D* h = dynamic_cast<TH2D*>(f->Get(histName.c_str()));
    if (!h) {
        std::cerr << "[DatasetCreator] Histogram '" << histName << "' not found\n";
        f->Close(); delete f; return;
    }
    _efficiencyMap.reset(dynamic_cast<TH2D*>(h->Clone("efficiencyMap")));
    _efficiencyMap->SetDirectory(nullptr);
    f->Close(); delete f;
    std::cout << "[DatasetCreator] Efficiency map: " << histName
              << " (" << _efficiencyMap->GetNbinsX()
              << "x"  << _efficiencyMap->GetNbinsY() << ")\n";
}

double DatasetCreator::calculateWeight(float pt, float y) const
{
    if (!_efficiencyMap) return 1.0;
    if (pt < _ptMin || pt > _ptMax ||
        y  < _efficiencyMap->GetYaxis()->GetXmin() ||
        y  > _efficiencyMap->GetYaxis()->GetXmax()) return 1.0;

    int ptBin = _efficiencyMap->GetXaxis()->FindBin(pt);
    int yBin  = _efficiencyMap->GetYaxis()->FindBin(y);
    if (ptBin < 1 || ptBin > _efficiencyMap->GetXaxis()->GetNbins() ||
        yBin  < 1 || yBin  > _efficiencyMap->GetYaxis()->GetNbins()) return 1.0;

    double eff = _efficiencyMap->GetBinContent(
                     _efficiencyMap->GetBin(ptBin, yBin));
    return (eff > 1e-4) ? 1.0 / eff : 1.0;
}

void DatasetCreator::createVariables()
{
    _m1 = std::make_unique<RooRealVar>("m1","m_{1} [GeV/c^{2}]",_m1Min,_m1Max);
    _m2 = std::make_unique<RooRealVar>("m2","m_{2} [GeV/c^{2}]",_m2Min,_m2Max);
    _m1->setBins(50); _m2->setBins(50);

    _y1   = std::make_unique<RooRealVar>("y1",  "y_{1}",    -2.0, 2.0);
    _y2   = std::make_unique<RooRealVar>("y2",  "y_{2}",    -2.0, 2.0);
    _pt1  = std::make_unique<RooRealVar>("pt1", "p_{T1}",    0.0, 50.0);
    _pt2  = std::make_unique<RooRealVar>("pt2", "p_{T2}",    0.0, 50.0);
    _pairPt   = std::make_unique<RooRealVar>("pairPt",  "p_{T pair}", 0.0, 100.0);
    _pairMass = std::make_unique<RooRealVar>("pairMass","m_{pair}",   0.0, 100.0);
    _phi1 = std::make_unique<RooRealVar>("phi1","#phi_{1}",
                                          -2*TMath::Pi(), 2*TMath::Pi());
    _phi2 = std::make_unique<RooRealVar>("phi2","#phi_{2}",
                                          -2*TMath::Pi(), 2*TMath::Pi());
    _y1->setBins(50);  _y2->setBins(50);
    _phi1->setBins(50); _phi2->setBins(50);
    _pt1->setBins(100); _pt2->setBins(100);
    _pairPt->setBins(100); _pairMass->setBins(100);
}

// ── private helper ──────────────────────────────────────────────────────────
static std::unique_ptr<RooDataSet>
fillDataset(TTree* tree, const std::string& dsName, bool isOS,
            RooRealVar& m1, RooRealVar& m2,
            RooRealVar& y1, RooRealVar& y2,
            RooRealVar& phi1, RooRealVar& phi2,
            RooRealVar& pt1, RooRealVar& pt2,
            RooRealVar& pairPt, RooRealVar& pairMass,
            std::function<double(float,float)> weightFn,
            bool useWeights)
{
    RooArgSet vars(m1,m2,y1,y2,phi1,phi2,pt1,pt2,pairPt,pairMass);

    std::unique_ptr<RooDataSet> ds;
    if (useWeights) {
        // Store each candidate's individual weight as a column so that
        // downstream code can inspect them if needed.  The combined weight
        // w1*w2 is passed to RooDataSet::add() and stored internally.
        RooRealVar wCol("_w","event weight",1.,0.,1e9);
        RooArgSet wsVars(vars); wsVars.add(wCol);
        ds = std::make_unique<RooDataSet>(dsName.c_str(),
                                          (dsName+" dataset").c_str(),
                                          wsVars,
                                          RooFit::WeightVar("_w"));
    } else {
        ds = std::make_unique<RooDataSet>(dsName.c_str(),
                                          (dsName+" dataset").c_str(),
                                          vars);
    }

    float m1v,m2v,mbar1v,mbar2v;
    float pt1v,pt2v,y1v,y2v,phi1v,phi2v;
    double pairPtV,pairMassV;
    UChar_t pairTypeV;

    tree->SetBranchAddress("m1",      &m1v);
    tree->SetBranchAddress("m2",      &m2v);
    tree->SetBranchAddress("mbar1",   &mbar1v);
    tree->SetBranchAddress("mbar2",   &mbar2v);
    tree->SetBranchAddress("pairType",&pairTypeV);
    tree->SetBranchAddress("y1",      &y1v);
    tree->SetBranchAddress("y2",      &y2v);
    tree->SetBranchAddress("phi1",    &phi1v);
    tree->SetBranchAddress("phi2",    &phi2v);
    tree->SetBranchAddress("pt1",     &pt1v);
    tree->SetBranchAddress("pt2",     &pt2v);
    tree->SetBranchAddress("pairPt",  &pairPtV);
    tree->SetBranchAddress("pairMass",&pairMassV);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        double mass1, mass2;
        if (isOS) {
            if      (TESTBIT(pairTypeV,DDbar))    { mass1=m1v;    mass2=mbar2v; }
            else if (TESTBIT(pairTypeV,DbarD))    { mass1=mbar1v; mass2=m2v;   }
            else continue;
        } else {
            if      (TESTBIT(pairTypeV,DD))       { mass1=m1v;    mass2=m2v;   }
            else if (TESTBIT(pairTypeV,DbarDbar)) { mass1=mbar1v; mass2=mbar2v;}
            else continue;
        }

        m1.setVal(mass1);  m2.setVal(mass2);
        y1.setVal(y1v);    y2.setVal(y2v);
        phi1.setVal(phi1v); phi2.setVal(phi2v);
        pt1.setVal(pt1v);  pt2.setVal(pt2v);
        pairPt.setVal(pairPtV); pairMass.setVal(pairMassV);

        if (useWeights) {
            double w = weightFn(pt1v,y1v) * weightFn(pt2v,y2v);
            ds->add(vars, w);
        } else {
            ds->add(vars);
        }
    }
    return ds;
}

void DatasetCreator::createDatasets()
{
    if (!_m1) createVariables();
    bool uw = (_efficiencyMap != nullptr);
    auto wfn = [this](float pt, float y){ return calculateWeight(pt,y); };

    if (_osTree && _osTree->GetEntries() > 0) {
        _osData = fillDataset(_osTree,"osData",true,
                              *_m1,*_m2,*_y1,*_y2,*_phi1,*_phi2,
                              *_pt1,*_pt2,*_pairPt,*_pairMass, wfn, uw);
        std::cout << "[DatasetCreator] OS: " << _osData->numEntries() << " entries";
        if (uw) std::cout << "  sumW=" << _osData->sumEntries();
        std::cout << std::endl;
    }
    if (_lsTree && _lsTree->GetEntries() > 0) {
        _lsData = fillDataset(_lsTree,"lsData",false,
                              *_m1,*_m2,*_y1,*_y2,*_phi1,*_phi2,
                              *_pt1,*_pt2,*_pairPt,*_pairMass, wfn, uw);
        std::cout << "[DatasetCreator] LS: " << _lsData->numEntries() << " entries";
        if (uw) std::cout << "  sumW=" << _lsData->sumEntries();
        std::cout << std::endl;
    }
}

void DatasetCreator::saveDatasets(const std::string& outputFilename) const
{
    if (!_osData && !_lsData) {
        std::cerr << "[DatasetCreator] Nothing to save – call createDatasets() first\n";
        return;
    }
    TFile out(outputFilename.c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "[DatasetCreator] Cannot create: " << outputFilename << std::endl;
        return;
    }
    if (_osData) _osData->Write("osData");
    if (_lsData) _lsData->Write("lsData");
    if (_efficiencyMap) _efficiencyMap->Write("efficiencyMap");
    out.Close();
    std::cout << "[DatasetCreator] Saved to " << outputFilename << std::endl;
}

