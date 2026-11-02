#include "MESWeights.h"

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <memory>

double MESWeights::WeightFunction::eval(double m) const {
    if (!ok || m < mLo || m > mHi) return 0.;
    const double x = (m - mLo) / step;
    const int i = (int)x;
    if (i >= (int)w.size() - 1) return w.back();
    const double t = x - i;
    return w[i] * (1. - t) + w[i + 1] * t;
}

double MESWeights::signalFraction1(int b) const { const auto& f = _bins[b].f1; return f.nSig / (f.nSig + f.nBkg); }
double MESWeights::signalFraction2(int b) const { const auto& f = _bins[b].f2; return f.nSig / (f.nSig + f.nBkg); }
bool MESWeights::inRange1(int b, double m) const { const auto& f = _bins[b].f1; return m >= f.mLo && m <= f.mHi; }
bool MESWeights::inRange2(int b, double m) const { const auto& f = _bins[b].f2; return m >= f.mLo && m <= f.mHi; }

// tag = "m1" or "m2". Uses the model "prefit_me_<tag>_model" (signal + background).
bool MESWeights::build(RooWorkspace* ws, const std::string& tag, WeightFunction& out) {
    RooRealVar* mass = ws->var(tag.c_str());
    auto* model = dynamic_cast<RooAddPdf*>(ws->pdf(("prefit_me_" + tag + "_model").c_str()));
    if (!mass || !model) {
        std::cerr << "Error: variable '" << tag << "' or model prefit_me_" << tag << "_model not in workspace\n";
        return false;
    }
    ws->loadSnapshot(("meParams_" + tag).c_str());   // fitted values (also stored in the workspace itself)

    const RooArgList& pdfs  = model->pdfList();
    const RooArgList& coefs = model->coefList();
    if (pdfs.size() != 2 || coefs.size() != 2) {
        std::cerr << "Error: " << tag << " model has " << pdfs.size()
                  << " components; only signal + background is supported (no reflections)\n";
        return false;
    }
    auto* fS = dynamic_cast<RooAbsPdf*>(pdfs.at(0));
    auto* fB = dynamic_cast<RooAbsPdf*>(pdfs.at(1));
    auto* yS = dynamic_cast<RooAbsReal*>(coefs.at(0));
    auto* yB = dynamic_cast<RooAbsReal*>(coefs.at(1));
    if (!fS || !fB || !yS || !yB) return false;

    const double nS = yS->getVal(), nB = yB->getVal();
    if (nS <= 0. || nB <= 0.) {
        std::cerr << "Error: non-positive yields in " << tag << " fit (nSig=" << nS << ", nBkg=" << nB << ")\n";
        return false;
    }

    // Normalized component pdfs on a fine grid over the fit range
    const int N = 4000;
    const double lo = mass->getMin(), hi = mass->getMax(), h = (hi - lo) / N;
    RooArgSet nset(*mass);
    const double saved = mass->getVal();
    std::vector<double> fs(N + 1), fb(N + 1), D(N + 1);
    for (int i = 0; i <= N; i++) {
        mass->setVal(lo + i * h);
        fs[i] = fS->getVal(&nset);
        fb[i] = fB->getVal(&nset);
        D[i]  = nS * fs[i] + nB * fb[i];
    }
    mass->setVal(saved);

    auto integrate = [&](auto&& g) {   // trapezoid rule
        double s = 0.;
        for (int i = 0; i <= N; i++) s += ((i == 0 || i == N) ? 0.5 : 1.) * g(i);
        return s * h;
    };
    const double W00 = integrate([&](int i) { return D[i] > 0. ? fs[i] * fs[i] / D[i] : 0.; });
    const double W01 = integrate([&](int i) { return D[i] > 0. ? fs[i] * fb[i] / D[i] : 0.; });
    const double W11 = integrate([&](int i) { return D[i] > 0. ? fb[i] * fb[i] / D[i] : 0.; });
    const double det = W00 * W11 - W01 * W01;
    if (!(det > 0.)) {
        std::cerr << "Error: singular sPlot matrix for " << tag << " (signal and background pdfs too similar?)\n";
        return false;
    }
    const double V00 = W11 / det, V01 = -W01 / det;

    out.w.assign(N + 1, 0.);
    for (int i = 0; i <= N; i++) out.w[i] = D[i] > 0. ? (V00 * fs[i] + V01 * fb[i]) / D[i] : 0.;
    out.mLo = lo; out.mHi = hi; out.step = h;
    out.nSig = nS; out.nBkg = nB; out.ok = true;

    // Closure checks: expect nSig, 1 and 0
    const double cWD  = integrate([&](int i) { return out.w[i] * D[i]; });
    const double cWfs = integrate([&](int i) { return out.w[i] * fs[i]; });
    const double cWfb = integrate([&](int i) { return out.w[i] * fb[i]; });
    std::cout << "  " << tag << ": nSig=" << nS << " nBkg=" << nB
              << " | int(w*D)=" << cWD << " (nSig), int(w*fs)=" << cWfs << " (1), int(w*fb)=" << cWfb << " (0)\n";
    return true;
}

bool MESWeights::addBin(double ptLo, double ptHi, const std::string& fitFile) {
    std::unique_ptr<TFile> f(TFile::Open(fitFile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open fit file " << fitFile << std::endl;
        return false;
    }
    auto* ws = dynamic_cast<RooWorkspace*>(f->Get("w"));
    if (!ws) {
        std::cerr << "Error: no workspace 'w' in " << fitFile << std::endl;
        return false;
    }
    std::cout << "pT bin [" << ptLo << ", " << ptHi << "): " << fitFile << std::endl;
    Bin b;
    b.ptLo = ptLo; b.ptHi = ptHi;
    if (!build(ws, "m1", b.f1) || !build(ws, "m2", b.f2)) return false;
    _bins.push_back(std::move(b));
    return true;
}

int MESWeights::findBin(double pt) const {
    for (size_t i = 0; i < _bins.size(); i++) {
        const bool last = (i + 1 == _bins.size());
        if (pt >= _bins[i].ptLo && (pt < _bins[i].ptHi || (last && pt <= _bins[i].ptHi))) return (int)i;
    }
    return -1;
}

double MESWeights::weight1(double m, double pt) const {
    const int b = findBin(pt);
    return b < 0 ? 0. : _bins[b].f1.eval(m);
}
double MESWeights::weight2(double m, double pt) const {
    const int b = findBin(pt);
    return b < 0 ? 0. : _bins[b].f2.eval(m);
}

bool MESWeights::writeWeightedTree(const std::string& inFile, const std::string& treeName,
                                   const std::string& outFile, const CorrConfig& cc,
                                   const std::string& plotPrefix) const {
    if (_bins.empty()) { std::cerr << "Error: no pT bins loaded\n"; return false; }

    std::unique_ptr<TFile> fin(TFile::Open(inFile.c_str(), "READ"));
    if (!fin || fin->IsZombie()) { std::cerr << "Error: cannot open " << inFile << std::endl; return false; }
    TTree* tin = dynamic_cast<TTree*>(fin->Get(treeName.c_str()));
    if (!tin) { std::cerr << "Error: tree '" << treeName << "' not found in " << inFile << std::endl; return false; }

    float m1 = 0, pt1 = 0, m2 = 0, pt2 = 0;
    if (tin->SetBranchAddress("m1", &m1) < 0 || tin->SetBranchAddress("pt1", &pt1) < 0 ||
        tin->SetBranchAddress("m2", &m2) < 0 || tin->SetBranchAddress("pt2", &pt2) < 0) {
        std::cerr << "Error: branches m1/pt1/m2/pt2 not found (expected float)\n";
        return false;
    }
    // deltaEta/deltaPhi are read through their leaves, so float or double both work
    TLeaf* lDeta = tin->GetLeaf("deltaEta");
    TLeaf* lDphi = tin->GetLeaf("deltaPhi");
    if (!lDeta || !lDphi) { std::cerr << "Error: branches deltaEta/deltaPhi not found\n"; return false; }

    std::unique_ptr<TFile> fout(TFile::Open(outFile.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) { std::cerr << "Error: cannot create " << outFile << std::endl; return false; }
    fout->cd();

    const double pi = TMath::Pi();
    auto* hRaw = new TH2D("hME_raw", "ME pairs, unweighted;#Delta#eta;#Delta#phi (rad)",
                          cc.nDeta, -cc.detaMax, cc.detaMax, cc.nDphi, -pi / 2, 3 * pi / 2);
    auto* hSw  = new TH2D("hME_sw", "ME pairs, sWeighted;#Delta#eta;#Delta#phi (rad)",
                          cc.nDeta, -cc.detaMax, cc.detaMax, cc.nDphi, -pi / 2, 3 * pi / 2);
    hRaw->Sumw2();
    hSw->Sumw2();

    TTree* tout = tin->CloneTree(0);          // same branches, empty
    double sw1 = 0, sw2 = 0, swPair = 0;
    tout->Branch("sw1", &sw1, "sw1/D");
    tout->Branch("sw2", &sw2, "sw2/D");
    tout->Branch("swPair", &swPair, "swPair/D");

    const int nb = (int)_bins.size();
    std::vector<double> sumW1(nb, 0.), sumW2(nb, 0.);
    std::vector<std::vector<double>> nIn(nb, std::vector<double>(nb, 0.));
    double sumPair = 0.;
    long nNoBin = 0;

    const Long64_t n = tin->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        tin->GetEntry(i);
        const int a = findBin(pt1), b = findBin(pt2);
        sw1 = (a >= 0) ? _bins[a].f1.eval(m1) : 0.;
        sw2 = (b >= 0) ? _bins[b].f2.eval(m2) : 0.;
        swPair = sw1 * sw2;

        if (a < 0 || b < 0) nNoBin++;
        if (a >= 0) sumW1[a] += sw1;
        if (b >= 0) sumW2[b] += sw2;
        sumPair += swPair;
        if (a >= 0 && b >= 0) {
            const auto& fa = _bins[a].f1; const auto& fb = _bins[b].f2;
            if (m1 >= fa.mLo && m1 <= fa.mHi && m2 >= fb.mLo && m2 <= fb.mHi) nIn[a][b] += 1.;
        }

        const double deta = lDeta->GetValue(), dphi = lDphi->GetValue();
        hRaw->Fill(deta, dphi);
        if (swPair != 0.) hSw->Fill(deta, dphi, swPair);

        tout->Fill();
    }

    // Normalize the weighted ME to 1 on the plateau around (0, 0)
    auto* hNorm = (TH2D*)hSw->Clone("hME_sw_norm");
    hNorm->SetTitle("ME pairs, sWeighted and normalized;#Delta#eta;#Delta#phi (rad)");
    double plateau = 0.; int nPlateau = 0;
    for (int ix = 1; ix <= hSw->GetNbinsX(); ix++)
        for (int iy = 1; iy <= hSw->GetNbinsY(); iy++) {
            if (std::fabs(hSw->GetXaxis()->GetBinCenter(ix)) < cc.normDeta &&
                std::fabs(hSw->GetYaxis()->GetBinCenter(iy)) < cc.normDphi) {
                plateau += hSw->GetBinContent(ix, iy);
                nPlateau++;
            }
        }
    if (nPlateau > 0 && plateau > 0.) {
        hNorm->Scale(nPlateau / plateau);
        std::cout << "\nME normalization: mean plateau content = " << plateau / nPlateau
                  << " (" << nPlateau << " bins)\n";
    } else {
        std::cerr << "Warning: non-positive plateau (" << plateau << "), normalized ME not scaled\n";
    }

    // Sum-rule checks
    std::cout << "\nChecks (" << n << " pairs):\n";
    double expected = 0.;
    for (int a = 0; a < nb; a++) {
        std::cout << "  bin [" << _bins[a].ptLo << ", " << _bins[a].ptHi << "): sum sw1 = " << sumW1[a]
                  << " (fit nSig1 = " << _bins[a].f1.nSig << "), sum sw2 = " << sumW2[a]
                  << " (fit nSig2 = " << _bins[a].f2.nSig << ")\n";
        for (int b = 0; b < nb; b++) {
            const double pi1 = _bins[a].f1.nSig / (_bins[a].f1.nSig + _bins[a].f1.nBkg);
            const double pi2 = _bins[b].f2.nSig / (_bins[b].f2.nSig + _bins[b].f2.nBkg);
            expected += nIn[a][b] * pi1 * pi2;
        }
    }
    std::cout << "  sum swPair = " << sumPair << " (expected ~ " << expected << ")\n";
    if (nNoBin > 0) std::cout << "  " << nNoBin << " pairs have a candidate outside all pT bins (weight 0)\n";

    // Summary plot
    if (!plotPrefix.empty()) {
        gStyle->SetOptStat(0);
        TCanvas c("c_me_deltaEtaPhi", "ME correlation", 1800, 1000);
        c.Divide(3, 2);
        c.cd(1); hRaw->Draw("COLZ");
        c.cd(2); hSw->Draw("COLZ");
        c.cd(3); hNorm->Draw("COLZ");
        c.cd(4);
        TH1D* pX = hSw->ProjectionX("pME_deta");
        pX->SetTitle("sWeighted ME, #Delta#eta projection;#Delta#eta;weighted pairs");
        pX->SetLineWidth(2); pX->Draw("E");
        c.cd(5);
        TH1D* pY = hSw->ProjectionY("pME_dphi");
        pY->SetTitle("sWeighted ME, #Delta#phi projection;#Delta#phi (rad);weighted pairs");
        pY->SetLineWidth(2); pY->Draw("E");
        c.SaveAs((plotPrefix + "_deltaEtaPhi.pdf").c_str());
        c.SaveAs((plotPrefix + "_deltaEtaPhi.png").c_str());
    }

    fout->cd();
    tout->Write("", TObject::kOverwrite);
    hRaw->Write();
    hSw->Write();
    hNorm->Write();
    fout->Close();
    std::cout << "Weighted ME tree and histograms written to " << outFile << std::endl;
    return true;
}

void MESWeights::plotWeightFunctions(const std::string& outputPrefix) const {
    TCanvas c("c_sw_me", "ME sWeights", 1600, 700);
    c.Divide(2, 1);
    for (int cand = 1; cand <= 2; cand++) {
        c.cd(cand);
        auto* mg = new TMultiGraph();
        auto* leg = new TLegend(0.15, 0.55, 0.45, 0.88);
        for (size_t i = 0; i < _bins.size(); i++) {
            const WeightFunction& f = (cand == 1) ? _bins[i].f1 : _bins[i].f2;
            if (!f.ok) continue;
            auto* g = new TGraph();
            for (size_t k = 0; k < f.w.size(); k += 10) g->SetPoint(g->GetN(), f.mLo + k * f.step, f.w[k]);
            g->SetLineColor(1 + (int)(i % 9));
            g->SetLineWidth(2);
            mg->Add(g, "L");
            leg->AddEntry(g, Form("p_{T} [%g, %g)", _bins[i].ptLo, _bins[i].ptHi), "l");
        }
        mg->SetTitle(Form("signal sWeight, candidate %d;m (GeV/c^{2});w_{s}", cand));
        mg->Draw("A");
        leg->Draw();
    }
    c.SaveAs((outputPrefix + "_sweights.pdf").c_str());
    c.SaveAs((outputPrefix + "_sweights.png").c_str());
}