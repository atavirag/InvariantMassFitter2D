// Mixed-event counterpart of SPlotter: sWeighted ME histograms from the 1D ME prefits.
#pragma once
#include "MESWeights.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <memory>
#include <string>
#include <vector>

class MESPlotter {
public:
    explicit MESPlotter(const std::string& meTreeFile, const std::string& treeName = "pairs");

    // One call per candidate-pT bin: loads the m1 and m2 prefit workspaces
    bool addBin(double ptLo, double ptHi, const std::string& fitFile) { return _weights.addBin(ptLo, ptHi, fitFile); }

    void addVariable(const std::string& name, int nBins, double lo, double hi);
    void addVariable2D(const std::string& nx, int nbx, double lox, double hix,
                       const std::string& ny, int nby, double loy, double hiy);

    // Single pass over the ME tree; fills h_<tag>_<var>_<species> and h2_<tag>_<x>_<y>_<species>
    bool run(const std::string& tag = "me");

    // Clone of a 2D histogram scaled to 1 on the plateau around (0, 0); stored as <name>_norm
    TH2D* addNormalized2D(const std::string& name, double plateauX, double plateauY);
    TH2D* find2D(const std::string& name) const;

    void save(const std::string& outFile);

    struct Var   { std::string name; int nBins; double min, max; };
    struct Var2D { std::string nameX; int nBinsX; double minX, maxX;
                   std::string nameY; int nBinsY; double minY, maxY; };

private:

    std::unique_ptr<TFile> _file;
    TTree* _tree = nullptr;
    MESWeights _weights;
    std::vector<Var> _vars;
    std::vector<Var2D> _vars2D;
    std::vector<std::unique_ptr<TH1D>> _histos;
    std::vector<std::unique_ptr<TH2D>> _histos2D;
};