// Class to extract the s-plot distributions from the fit results of the previous step
// First version: using raw distributions without reflections

#include "SPlotter.h"
#include <cmath>
#include <iostream>
#include <functional>
#include <map>

#include "TMath.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooStats/SPlot.h"

namespace {
const std::vector<std::pair<std::string, std::string>> kSpecies = {
    {"nSgnSgn", "SgnSgn"}, {"nSgnBkg", "SgnBkg"},
    {"nBkgSgn", "BkgSgn"}, {"nBkgBkg", "BkgBkg"}};

double foldedDeltaPhi(double a, double b) {          // |Δφ| folded into [0, π]
    double d = std::fabs(a - b);
    while (d > TMath::TwoPi()) d -= TMath::TwoPi();
    return d > TMath::Pi() ? TMath::TwoPi() - d : d;
}

double deltaPhiCorrelations(double a, double b) {              // (a - b) wrapped into [-π/2, 3π/2)
    const double lo = -TMath::PiOver2();
    double d = std::fmod(a - b - lo, TMath::TwoPi());
    if (d < 0) d += TMath::TwoPi();
    return d + lo;
}

using Fn = std::function<double(const RooArgSet&)>;
const std::map<std::string, Fn>& derived() {
    static const std::map<std::string, Fn> m = {
        {"deltaY",               [](const RooArgSet& r) { return r.getRealValue("y1") - r.getRealValue("y2"); }},
        {"deltaEta",             [](const RooArgSet& r) { return r.getRealValue("eta1") - r.getRealValue("eta2"); }},
        {"deltaPhi",             [](const RooArgSet& r) { return foldedDeltaPhi(r.getRealValue("phi1"), r.getRealValue("phi2")); }},
        {"deltaPhiCorrelations", [](const RooArgSet& r) { return deltaPhiCorrelations(r.getRealValue("phi1"), r.getRealValue("phi2")); }},
        {"xJ",                   [](const RooArgSet& r) {
            double a = r.getRealValue("pt1"), b = r.getRealValue("pt2");
            return std::min(a, b) / std::max(a, b); }},
    };
    return m;
}
double valueOf(const std::string& name, const RooArgSet& r) {
    auto it = derived().find(name);
    return it != derived().end() ? it->second(r) : r.getRealValue(name.c_str());
}
bool isKnown(const std::string& name, const RooArgSet& r) {
    return derived().count(name) || r.find(name.c_str());
}
} // namespace

SPlotter::SPlotter(const std::string& f) {
    _file.reset(TFile::Open(f.c_str(), "READ"));
    if (!_file || _file->IsZombie()) { std::cerr << "Cannot open " << f << "\n"; return; }
    _ws = dynamic_cast<RooWorkspace*>(_file->Get("w"));
    if (!_ws) std::cerr << "No workspace 'w' in " << f << "\n";
}

void SPlotter::addVariable(const std::string& n, int nb, double lo, double hi) {
    _vars.push_back({n, nb, lo, hi});
}

void SPlotter::addVariable2D(const std::string& nx, int nbx, double lox, double hix,
                             const std::string& ny, int nby, double loy, double hiy) {
    _vars2D.push_back({nx, nbx, lox, hix, ny, nby, loy, hiy});
}

bool SPlotter::run(const std::string& tag) {
    // load the workspace
    if (!_ws) return false;
    _ws->data("osData")->get()->Print();   // eta1, eta2 should be listed
    auto* data = dynamic_cast<RooDataSet*>(_ws->data((tag + "Data").c_str()));
    auto* pdf  = dynamic_cast<RooAddPdf*>(_ws->pdf("model"));
    if (!data || !pdf) { std::cerr << "Missing " << tag << "Data or model\n"; return false; }
    if (data->isWeighted()) { std::cerr << "Dataset is weighted: rerun the fit without efficiency map\n"; return false; }
    if (pdf->pdfList().getSize() != (int)kSpecies.size()) {
        std::cerr << "Model has " << pdf->pdfList().getSize()
                  << " components, expected 4: rerun the fit without reflections\n";
        return false;
    }
    if (!_ws->loadSnapshot((tag + "Params").c_str())) { std::cerr << "No snapshot for " << tag << "\n"; return false; }

    RooArgList yields;
    for (const auto& s : kSpecies) yields.add(*_ws->var(s.first.c_str()));

    // Freeze everything except the four yields
    std::unique_ptr<RooArgSet> pars(pdf->getParameters(*data));
    for (RooAbsArg* a : *pars)
        if (auto* v = dynamic_cast<RooRealVar*>(a))
            if (!yields.find(v->GetName())) v->setConstant(true);

    RooStats::SPlot splot("splot", "sPlot", *data, pdf, yields);
    RooDataSet* sd = splot.GetSDataSet();

    // Sanity checks
    std::vector<double> sumSw(kSpecies.size(), 0.);
    double maxDev = 0.;
    for (int i = 0; i < sd->numEntries(); ++i) {
        const RooArgSet* r = sd->get(i);
        double tot = 0.;
        for (size_t k = 0; k < kSpecies.size(); ++k) {
            double sw = r->getRealValue((kSpecies[k].first + "_sw").c_str());
            sumSw[k] += sw; tot += sw;
        }
        maxDev = std::max(maxDev, std::fabs(tot - 1.));
    }
    std::cout << "\n[" << tag << "] sPlot checks (max |sum_species sw - 1| = " << maxDev << ")\n";
    for (size_t k = 0; k < kSpecies.size(); ++k)
        std::cout << "  " << kSpecies[k].second << ": sum sw = " << sumSw[k]
                  << "  fitted yield = " << _ws->var(kSpecies[k].first.c_str())->getVal() << "\n";

    // Fill raw sWeighted histograms
    for (const SPlotVar& v : _vars) {
        if (!isKnown(v.name, *sd->get(0))) { std::cerr << "Variable " << v.name << " not available, skipped\n"; continue; }
        std::vector<TH1D*> hs;
        for (const auto& s : kSpecies) {
            auto h = std::make_unique<TH1D>(("h_" + tag + "_" + v.name + "_" + s.second).c_str(),
                                            (tag + " " + s.second + ";" + v.name + ";sWeighted entries").c_str(),
                                            v.nBins, v.min, v.max);
            h->Sumw2();
            h->SetDirectory(nullptr);
            hs.push_back(h.get());
            _histos.push_back(std::move(h));
        }
        for (int i = 0; i < sd->numEntries(); ++i) {
            const RooArgSet* r = sd->get(i);
            double x = valueOf(v.name, *r);
            for (size_t k = 0; k < kSpecies.size(); ++k)
                hs[k]->Fill(x, r->getRealValue((kSpecies[k].first + "_sw").c_str()));
        }
    }
        // Fill raw sWeighted 2D histograms
    for (const SPlotVar2D& v : _vars2D) {
        const RooArgSet& first = *sd->get(0);
        if (!isKnown(v.nameX, first) || !isKnown(v.nameY, first)) {
            std::cerr << "Variables " << v.nameX << " / " << v.nameY << " not available, skipped\n";
            continue;
        }
        std::vector<TH2D*> hs;
        for (const auto& s : kSpecies) {
            auto h = std::make_unique<TH2D>(
                ("h2_" + tag + "_" + v.nameX + "_" + v.nameY + "_" + s.second).c_str(),
                (tag + " " + s.second + ";" + v.nameX + ";" + v.nameY + ";sWeighted entries").c_str(),
                v.nBinsX, v.minX, v.maxX, v.nBinsY, v.minY, v.maxY);
            h->Sumw2();
            h->SetDirectory(nullptr);
            hs.push_back(h.get());
            _histos2D.push_back(std::move(h));
        }
        for (int i = 0; i < sd->numEntries(); ++i) {
            const RooArgSet* r = sd->get(i);
            double x = valueOf(v.nameX, *r);
            double y = valueOf(v.nameY, *r);
            for (size_t k = 0; k < kSpecies.size(); ++k)
                hs[k]->Fill(x, y, r->getRealValue((kSpecies[k].first + "_sw").c_str()));
        }
    }
    return true;
}

void SPlotter::save(const std::string& outFile) {
    std::unique_ptr<TFile> out(TFile::Open(outFile.c_str(), "RECREATE"));
    for (auto& h : _histos) h->Write();
    for (auto& h : _histos2D) h->Write();
}