// InvariantMassFitter.h  –  Step 2: load pre-built datasets and fit.
//
// No tree reading here.  Datasets come from DatasetCreator output files.
//
// Usage:
//   InvariantMassFitter fitter;
//   fitter.setMassRange(...);
//   fitter.setFittingFunctions("gaussian","chebychev","gaussian");
//   fitter.loadDatasets("period_A.root");
//   fitter.mergeDatasets("period_B.root");   // optional
//   fitter.buildModel();
//   fitter.performPrefit();
//   fitter.fit();
//   fitter.plotFit("out");
//   fitter.saveResultsToFile();

#ifndef INVARIANTMASSFITTER_H
#define INVARIANTMASSFITTER_H

#include <memory>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "TString.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"

class TH2D;

class InvariantMassFitter {
public:
    InvariantMassFitter();
    ~InvariantMassFitter() = default;

    // Configuration
    void setMassRange(double m1Min, double m1Max, double m2Min, double m2Max);
    void setReflOverSgn(float reflOverSgn);
    void setFittingFunctions(TString const& sgnFunc,
                             TString const& bkgFunc,
                             TString const& reflFunc);
    void setOutputFile(const std::string& filename);

    // Dataset loading
    void loadDatasets (const std::string& filename);  // replaces current data
    void mergeDatasets(const std::string& filename);  // appends to current data

    // Core pipeline
    void createVariables();
    void buildModel(bool doReflections = false);
    void fit(const std::string& dataType = "both");   // "OS", "LS", "both"

    // Prefit
    void performPrefit   (const std::string& dataType = "both",
                          bool useReflections = false);
    void setPrefitParameters(bool doReflections = false);
    void plotPrefitResults  (const std::string& outputPrefix);

    // Output
    void plotFit          (const std::string& outputPrefix, bool doReflections = false);
    void saveWorkspace    (const std::string& filename);
    void saveResultsToFile(bool doReflections = false);
    void printFitResults  ();

    double computeChi2(RooPlot* frame, RooDataSet* data, int nFitParams);

    // Accessors
    RooDataSet* osData() const { return _osData.get(); }
    RooDataSet* lsData() const { return _lsData.get(); }

private:
    RooAbsPdf*  createSignalPdf    (const std::string& name, RooRealVar& mass);
    RooAbsPdf*  createBackgroundPdf(const std::string& name, RooRealVar& mass);
    RooAbsPdf*  createReflectedPdf (const std::string& name, RooRealVar& mass);
    void        plotModelOnFrame   (RooPlot* frame, RooAbsPdf* pdf, bool doReflections);
    RooRealVar* getYieldInRange    (RooFitResult* fitResult, RooProdPdf* function,
                                    RooFormulaVar* nCands, TString range);

    RooFitResult* performPrefit1D(RooDataSet* data, RooRealVar& mass,
                                  const std::string& name, bool useReflections);
    void extractPrefitParameters (const std::string& prefitName,
                                  const std::string& targetName);

    /// Internal: read datasets from an open TFile; optionally append.
    void readDatasetsFromFile(TFile* f, bool append);

    double _m1Min{1.7}, _m1Max{2.1};
    double _m2Min{1.7}, _m2Max{2.1};
    float  _reflOverSgn{0.2f};
    bool   _datasetsAreWeighted{false};

    std::unique_ptr<RooRealVar> _m1, _m2;
    std::unique_ptr<RooRealVar> _y1, _y2, _phi1, _phi2;
    std::unique_ptr<RooRealVar> _pt1, _pt2, _pairPt, _pairMass;

    std::unique_ptr<RooFitResult> _prefitResult_os_m1, _prefitResult_os_m2;
    std::unique_ptr<RooFitResult> _prefitResult_ls_m1, _prefitResult_ls_m2;

    std::unique_ptr<RooDataSet>  _osData, _lsData;
    std::unique_ptr<RooWorkspace> _workspace;

    std::string _outputRootFile;
    TString _sgnFuncOption{"gaussian"};
    TString _bkgFuncOption{"chebychev"};
    TString _reflFuncOption{"gaussian"};
};

#endif // INVARIANTMASSFITTER_H