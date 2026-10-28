// InvariantMassFitter.h - Updated with efficiency correction and output file support
#ifndef INVARIANTMASSFITTER_H
#define INVARIANTMASSFITTER_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "TString.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"

class TH2D;

class InvariantMassFitter {
public:
    // Public (not private) so that the ROOT/ACLiC dictionary can see it: mean pair weights per class
    // (order SgnSgn, SgnBkg, BkgSgn, BkgBkg) for the hybrid correction
    struct HybridWeights { std::vector<double> w; std::vector<double> err; bool ok = false; };

    InvariantMassFitter(const std::string& osFilename, const std::string& lsFilename);
    ~InvariantMassFitter();

    // True if both input files and both trees were opened correctly
    bool isValid() const { return _valid; }

    // Configuration methods
    void setMassRange(double m1Min, double m1Max, double m2Min, double m2Max);
    void setPtRange(double ptMin, double ptMax);
    void setReflOverSgn(float const& relfOverSgn);
    void setFittingFunctions(TString const &sgnFunc, TString const &bkgFunc, TString const &reflFunc);
    void setOutputFile(const std::string& filename);

    // Efficiency correction (returns false if the map could not be loaded)
    bool loadEfficiencyMap(const std::string& filename, const std::string& histName);

    // Hybrid mode: the fit is done on the RAW (unweighted) data, and the 3-sigma yields are then
    // multiplied by per-class mean pair weights <1/(eff1*eff2)> estimated from the data by a
    // linear fit over mass-sideband/window regions. Needs loadEfficiencyMap(). No reflections.
    void enableHybrid() { _hybrid = true; }

    // Candidates whose map efficiency is below this value are dropped (default 1e-4).
    // Raising it removes the huge weights 1/eff that dominate the weighted fit.
    void setMinEfficiency(double e) { _minEff = e; }

    // Fix the signal shape parameters (mean, sigma, ...) to the values of the raw ("basic") fit
    // stored in `basicFile` (OS values for the OS fit, LS values for the LS fit).
    // The background parameters stay free. Returns false if the file/snapshots cannot be read.
    bool fixSignalShapeFromFile(const std::string& basicFile);

    // Fix mean and sigma of the weighted LS distribution to the OS values (they have better stats)
    void setLSShapeFromOS(bool b) { _lsShapeFromOS = b; }

    // Core functionality
    void createVariables();
    void loadData();
    void buildModel(bool const doReflections = false);
    void fit(const std::string& dataType = "both"); // "OS", "LS", or "both"

    // Prefit
    void performPrefit(const std::string& dataType = "both", bool useReflections = false);
    void setPrefitParameters(bool const doReflections = false);
    void plotPrefitResults(const std::string& outputPrefix);

    // Output and visualization
    void plotFit(const std::string& outputPrefix, bool const doReflections = false);
    void saveWorkspace(const std::string& filename);
    void saveResultsToFile(bool const doReflections = false);
    void printFitResults();

    // Chi2 computation (works for weighted and unweighted data)
    double computeChi2(RooPlot* frame, RooDataSet* data, int nFitParams);

private:
    // Helper methods
    RooAbsPdf* createSignalPdf(const std::string& name, RooRealVar& mass);
    RooAbsPdf* createBackgroundPdf(const std::string& name, RooRealVar& mass);
    RooAbsPdf* createReflectedPdf(const std::string& name, RooRealVar& mass);
    // Returns 1/efficiency, or -1 if the candidate cannot be corrected
    double calculateWeights(float const &pt, float const &y);
    void plotModelOnFrame(RooPlot* frame, RooAbsPdf* pdf, bool const doReflections,
                          const std::vector<std::string>& textLines = {}, bool showLegend = true);

    // 3-sigma yields
    // value/error: reported yield (hybrid-corrected in hybrid mode); rawValue/rawError: from the fit
    struct YieldEntry { std::string label; double value; double error; bool ok; double rawValue; double rawError; };

    // Hybrid correction: mean pair weights per class (order SgnSgn, SgnBkg, BkgSgn, BkgBkg)
    bool _hybrid = false;
    bool weighted() const { return _efficiencyMap != nullptr && !_hybrid; }
    std::map<std::string, HybridWeights> _hybridW; //! keyed by "os"/"ls" (transient: not for ROOT I/O)
    bool computeHybridWeights(const std::string& tag);
    bool setSigmaWindow();
    std::vector<YieldEntry> computeYields3Sigma(const std::string& snapshot, RooFitResult* res, bool doReflections);
    std::vector<std::string> yieldTextLines(const std::string& label, const std::string& snapshot,
                                            const std::string& fitResultName, bool doReflections);
    RooRealVar *getYieldInRange(RooFitResult *fitResult, RooAbsPdf* function, RooAbsReal* nCands, TString range);

    // Input files and trees
    std::unique_ptr<TFile> _osFile;
    std::unique_ptr<TFile> _lsFile;
    TTree* _osTree;
    TTree* _lsTree;
    bool _valid;

    // Mass variables
    std::unique_ptr<RooRealVar> _m1;
    std::unique_ptr<RooRealVar> _m2;

    // Kinematic variables
    std::unique_ptr<RooRealVar> _y1, _y2;
    std::unique_ptr<RooRealVar> _eta1, _eta2;
    std::unique_ptr<RooRealVar> _phi1, _phi2;
    std::unique_ptr<RooRealVar> _pt1, _pt2, _pairPt, _pairMass;

    double _m1Min, _m1Max;
    double _m2Min, _m2Max;
    double _ptMin, _ptMax;
    float _reflOverSgn;

    // Prefit
    std::unique_ptr<RooFitResult> _prefitResult_os_m1;
    std::unique_ptr<RooFitResult> _prefitResult_os_m2;
    std::unique_ptr<RooFitResult> _prefitResult_ls_m1;
    std::unique_ptr<RooFitResult> _prefitResult_ls_m2;

    RooFitResult* performPrefit1D(RooDataSet* data, RooRealVar& mass,
                                  const std::string& name, bool useReflections);
    void extractPrefitParameters(const std::string& prefitName, const std::string& targetName);

    // Datasets
    std::unique_ptr<RooDataSet> _osData;
    std::unique_ptr<RooDataSet> _lsData;

    // Workspace
    std::unique_ptr<RooWorkspace> _workspace;

    // Efficiency correction
    std::unique_ptr<TH2D> _efficiencyMap;
    double _minEff = 1e-4;

    // Signal shape taken from the raw fit: "os"/"ls" -> parameter name -> value
    std::map<std::string, std::map<std::string, double>> _fixedShape;
    void applyFixedShape(const std::string& tag);

    bool _lsShapeFromOS = false;
    std::vector<std::string> signalShapeNames() const;

    // Output
    std::string _outputRootFile;

    // Function options
    TString _sgnFuncOption;
    TString _bkgFuncOption;
    TString _reflFuncOption;
};

#endif // INVARIANTMASSFITTER_H