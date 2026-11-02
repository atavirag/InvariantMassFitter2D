#include "MESPlotter.h"

#include <TLeaf.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>

namespace {

struct Species { const char* name; bool sig1; bool sig2; };
const Species kSpecies[4] = {{"SgnSgn", true, true}, {"SgnBkg", true, false},
                             {"BkgSgn", false, true}, {"BkgBkg", false, false}};
constexpr int kNW = 5;   // four species + "All" (unweighted)
const char* speciesName(int k) { return k < 4 ? kSpecies[k].name : "All"; }

// Same definitions as in SPlotter.cxx
double foldedDeltaPhi(double a, double b) {   // |dphi| folded into [0, pi]
    double d = std::fabs(a - b);
    while (d > TMath::TwoPi()) d -= TMath::TwoPi();
    return d > TMath::Pi() ? TMath::TwoPi() - d : d;
}
double deltaPhiCorrelations(double a, double b) {   // (a - b) wrapped into [-pi/2, 3pi/2)
    const double lo = -TMath::PiOver2();
    double d = std::fmod(a - b - lo, TMath::TwoPi());
    if (d < 0) d += TMath::TwoPi();
    return d + lo;
}
double xJ(double a, double b) { return std::min(a, b) / std::max(a, b); }

using Val = std::function<double()>;

// Returns a getter for a tree variable (derived or plain branch), or nullptr if unavailable
Val resolve(TTree* t, const std::string& name) {
    auto combine = [&](const char* a, const char* b, double (*f)(double, double)) -> Val {
        TLeaf *la = t->GetLeaf(a), *lb = t->GetLeaf(b);
        if (!la || !lb) return nullptr;
        return [la, lb, f] { return f(la->GetValue(), lb->GetValue()); };
    };
    if (name == "deltaY")   return combine("y1", "y2", [](double a, double b) { return a - b; });
    if (name == "deltaEta") return combine("eta1", "eta2", [](double a, double b) { return a - b; });
    if (name == "deltaPhi") return combine("phi1", "phi2", foldedDeltaPhi);
    if (name == "deltaPhiCorrelations") return combine("phi1", "phi2", deltaPhiCorrelations);
    if (name == "xJ")       return combine("pt1", "pt2", xJ);
    if (TLeaf* l = t->GetLeaf(name.c_str())) return [l] { return l->GetValue(); };
    return nullptr;
}

}  // namespace

MESPlotter::MESPlotter(const std::string& f, const std::string& treeName) {
    _file.reset(TFile::Open(f.c_str(), "READ"));
    if (!_file || _file->IsZombie()) { std::cerr << "Cannot open " << f << "\n"; return; }
    _tree = dynamic_cast<TTree*>(_file->Get(treeName.c_str()));
    if (!_tree) std::cerr << "No tree '" << treeName << "' in " << f << "\n";
}

void MESPlotter::addVariable(const std::string& n, int nb, double lo, double hi) {
    _vars.push_back({n, nb, lo, hi});
}

void MESPlotter::addVariable2D(const std::string& nx, int nbx, double lox, double hix,
                               const std::string& ny, int nby, double loy, double hiy) {
    _vars2D.push_back({nx, nbx, lox, hix, ny, nby, loy, hiy});
}

bool MESPlotter::run(const std::string& tag) {
    if (!_tree) return false;
    if (_weights.nBins() == 0) { std::cerr << "No pT bins loaded: call addBin first\n"; return false; }

    TLeaf *lm1 = _tree->GetLeaf("m1"), *lpt1 = _tree->GetLeaf("pt1");
    TLeaf *lm2 = _tree->GetLeaf("m2"), *lpt2 = _tree->GetLeaf("pt2");
    if (!lm1 || !lpt1 || !lm2 || !lpt2) { std::cerr << "Branches m1/pt1/m2/pt2 not found\n"; return false; }

    struct H1 { Val x; TH1D* h[kNW]; };
    struct H2 { Val x, y; TH2D* h[kNW]; };
    std::vector<H1> h1;
    std::vector<H2> h2;

    for (const Var& v : _vars) {
        Val x = resolve(_tree, v.name);
        if (!x) { std::cerr << "Variable " << v.name << " not available, skipped\n"; continue; }
        H1 e; e.x = x;
        for (int k = 0; k < kNW; k++) {
            auto h = std::make_unique<TH1D>(("h_" + tag + "_" + v.name + "_" + speciesName(k)).c_str(),
                                            (tag + " " + speciesName(k) + ";" + v.name + ";sWeighted entries").c_str(),
                                            v.nBins, v.min, v.max);
            h->Sumw2();
            h->SetDirectory(nullptr);
            e.h[k] = h.get();
            _histos.push_back(std::move(h));
        }
        h1.push_back(e);
    }
    for (const Var2D& v : _vars2D) {
        Val x = resolve(_tree, v.nameX), y = resolve(_tree, v.nameY);
        if (!x || !y) { std::cerr << "Variables " << v.nameX << " / " << v.nameY << " not available, skipped\n"; continue; }
        H2 e; e.x = x; e.y = y;
        for (int k = 0; k < kNW; k++) {
            auto h = std::make_unique<TH2D>(("h2_" + tag + "_" + v.nameX + "_" + v.nameY + "_" + speciesName(k)).c_str(),
                                            (tag + " " + speciesName(k) + ";" + v.nameX + ";" + v.nameY + ";sWeighted entries").c_str(),
                                            v.nBinsX, v.minX, v.maxX, v.nBinsY, v.minY, v.maxY);
            h->Sumw2();
            h->SetDirectory(nullptr);
            e.h[k] = h.get();
            _histos2D.push_back(std::move(h));
        }
        h2.push_back(e);
    }

    const int nb = _weights.nBins();
    std::vector<std::vector<double>> nIn(nb, std::vector<double>(nb, 0.));
    std::vector<double> sum(kNW, 0.);
    long nSkipped = 0;

    const Long64_t n = _tree->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        _tree->GetEntry(i);
        const double m1 = lm1->GetValue(), pt1 = lpt1->GetValue();
        const double m2 = lm2->GetValue(), pt2 = lpt2->GetValue();
        const int a = _weights.binOf(pt1), b = _weights.binOf(pt2);

        // Pairs outside the fitted pT bins or mass windows have no defined weight: skip them for all species
        if (a < 0 || b < 0 || !_weights.inRange1(a, m1) || !_weights.inRange2(b, m2)) { nSkipped++; continue; }

        const double s1 = _weights.weight1(m1, pt1), s2 = _weights.weight2(m2, pt2);
        const double b1 = 1. - s1, b2 = 1. - s2;      // w_sig + w_bkg = 1 per candidate
        const double w[kNW] = {s1 * s2, s1 * b2, b1 * s2, b1 * b2, 1.};

        nIn[a][b] += 1.;
        for (int k = 0; k < kNW; k++) sum[k] += w[k];
        for (auto& e : h1) { const double x = e.x(); for (int k = 0; k < kNW; k++) e.h[k]->Fill(x, w[k]); }
        for (auto& e : h2) { const double x = e.x(), y = e.y(); for (int k = 0; k < kNW; k++) e.h[k]->Fill(x, y, w[k]); }
    }

    // Sanity checks: expected sum of each species weight from the fitted signal fractions
    std::cout << "\n[" << tag << "] sPlot checks (" << n << " pairs, " << nSkipped << " skipped)\n";
    for (int k = 0; k < 4; k++) {
        double expected = 0.;
        for (int a = 0; a < nb; a++)
            for (int b = 0; b < nb; b++) {
                const double p1 = _weights.signalFraction1(a), p2 = _weights.signalFraction2(b);
                expected += nIn[a][b] * (kSpecies[k].sig1 ? p1 : 1. - p1) * (kSpecies[k].sig2 ? p2 : 1. - p2);
            }
        std::cout << "  " << kSpecies[k].name << ": sum sw = " << sum[k] << "  expected = " << expected << "\n";
    }
    std::cout << "  All: " << sum[4] << " pairs\n";
    return true;
}

TH2D* MESPlotter::find2D(const std::string& name) const {
    for (const auto& h : _histos2D) if (name == h->GetName()) return h.get();
    return nullptr;
}

TH2D* MESPlotter::addNormalized2D(const std::string& name, double plateauX, double plateauY) {
    TH2D* src = find2D(name);
    if (!src) { std::cerr << "Histogram " << name << " not found\n"; return nullptr; }

    double sum = 0.; int nPlateau = 0;
    for (int ix = 1; ix <= src->GetNbinsX(); ix++)
        for (int iy = 1; iy <= src->GetNbinsY(); iy++)
            if (std::fabs(src->GetXaxis()->GetBinCenter(ix)) < plateauX &&
                std::fabs(src->GetYaxis()->GetBinCenter(iy)) < plateauY) {
                sum += src->GetBinContent(ix, iy);
                nPlateau++;
            }
    if (nPlateau == 0 || sum <= 0.) {
        std::cerr << "Non-positive plateau for " << name << " (" << sum << "), not normalized\n";
        return nullptr;
    }
    std::unique_ptr<TH2D> h(static_cast<TH2D*>(src->Clone((name + "_norm").c_str())));
    h->SetDirectory(nullptr);
    h->Scale(nPlateau / sum);
    std::cout << name << ": plateau mean = " << sum / nPlateau << " (" << nPlateau << " bins)\n";
    _histos2D.push_back(std::move(h));
    return _histos2D.back().get();
}

void MESPlotter::save(const std::string& outFile) {
    std::unique_ptr<TFile> out(TFile::Open(outFile.c_str(), "RECREATE"));
    for (auto& h : _histos) h->Write();
    for (auto& h : _histos2D) h->Write();
}