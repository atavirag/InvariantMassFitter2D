// DatasetCreator.h
// Step 1: read OS/LS trees, apply efficiency weights, save RooDataSets.

#ifndef DATASETCREATOR_H
#define DATASETCREATOR_H

#include <memory>
#include <string>
#include <functional>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

class DatasetCreator {
public:
    DatasetCreator(const std::string& osFilename, const std::string& lsFilename);
    ~DatasetCreator() = default;

    // Configuration
    void setMassRange(double m1Min, double m1Max, double m2Min, double m2Max);
    void setPtRange(double ptMin, double ptMax);

    /// Load a (pt, y) efficiency TH2D used to compute 1/ε weights.
    /// If never called, all weights default to 1.
    void loadEfficiencyMap(const std::string& filename, const std::string& histName);

    // Core
    void createVariables();   // called automatically by createDatasets()
    void createDatasets();
    void saveDatasets(const std::string& outputFilename) const;

    // In-process access (optional)
    RooDataSet* osData() const { return _osData.get(); }
    RooDataSet* lsData() const { return _lsData.get(); }

private:
    double calculateWeight(float pt, float y) const;

    std::unique_ptr<TFile> _osFile, _lsFile;
    TTree* _osTree{nullptr};
    TTree* _lsTree{nullptr};

    double _m1Min{1.7}, _m1Max{2.1};
    double _m2Min{1.7}, _m2Max{2.1};
    double _ptMin{1.0}, _ptMax{24.0};

    std::unique_ptr<RooRealVar> _m1, _m2;
    std::unique_ptr<RooRealVar> _y1, _y2, _phi1, _phi2;
    std::unique_ptr<RooRealVar> _pt1, _pt2, _pairPt, _pairMass;

    std::unique_ptr<RooDataSet> _osData, _lsData;
    std::unique_ptr<TH2D>       _efficiencyMap;
};

#endif // DATASETCREATOR_H