#ifndef SPLOTTER_H
#define SPLOTTER_H

#include <memory>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "RooWorkspace.h"

struct SPlotVar { std::string name; int nBins; double min, max; };
struct SPlotVar2D {
    std::string nameX; int nBinsX; double minX, maxX;
    std::string nameY; int nBinsY; double minY, maxY;
};

class SPlotter {
public:
    explicit SPlotter(const std::string& fitResultsFile);
    void addVariable(const std::string& name, int nBins, double min, double max);
    void addVariable2D(const std::string& nx, int nbx, double lox, double hix,
                       const std::string& ny, int nby, double loy, double hiy);
    bool run(const std::string& tag);            // "os" or "ls"
    void save(const std::string& outFile);
private:
    std::unique_ptr<TFile> _file;
    RooWorkspace* _ws = nullptr;
    std::vector<SPlotVar> _vars;
    std::vector<std::unique_ptr<TH1D>> _histos;
    std::vector<SPlotVar2D> _vars2D;
    std::vector<std::unique_ptr<TH2D>> _histos2D;
};
#endif