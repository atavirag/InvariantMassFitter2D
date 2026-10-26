#ifndef SPLOTTER_H
#define SPLOTTER_H

#include <memory>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "RooWorkspace.h"

struct SPlotVar { std::string name; int nBins; double min, max; };

class SPlotter {
public:
    explicit SPlotter(const std::string& fitResultsFile);
    void addVariable(const std::string& name, int nBins, double min, double max);
    bool run(const std::string& tag);            // "os" or "ls"
    void save(const std::string& outFile);
private:
    std::unique_ptr<TFile> _file;
    RooWorkspace* _ws = nullptr;
    std::vector<SPlotVar> _vars;
    std::vector<std::unique_ptr<TH1D>> _histos;
};
#endif