/**
 * project_sparse.C
 *
 * Step 1 of 2 in the D0 significance optimisation pipeline.
 *
 * Reads the THnSparse from AnalysisResults.root, applies cuts on pT and
 * BDT background score (plus a rapidity fiducial cut), and saves one TH1D
 * per (pT bin, BDT cut) combination into projections.root.
 *
 * The Python script fit_projections.py then reads that file with uproot and
 * runs flarefly fits — no PyROOT required on the Python side.
 *
 * Sparse axis map:
 *   0  - Invariant mass
 *   1  - pT
 *   2  - BDT score background   ← cut axis
 *   3  - BDT score non-prompt
 *   4  - BDT score prompt
 *   5  - y
 *   6  - D0 type
 *   7  - Gap type
 *   8  - FT0-A amplitude
 *   9  - FT0-C amplitude
 *  10  - FV0-A amplitude
 *  11  - FDD-A amplitude
 *  12  - FDD-C amplitude
 *  13  - Number of primary vertex contributors
 *  14  - ZNA energy
 *  15  - ZNC energy
 *  16  - ZNA time
 *  17  - ZNC time
 *
 * Usage:
 *   root -l -b -q project_sparse.C
 */

#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// ─── Axis indices ─────────────────────────────────────────────────────────────
// Axes of hMassVsPtVsPtBVsYVsOriginVsD0Type (hf-task-d0)
namespace Axis {
  constexpr int kMass       =  0;
  constexpr int kPt         =  1;
  constexpr int kY          =  2;
  constexpr int kD0Type     =  3;
}

// Axes of hnDMeson (hf-correlator-d-meson-pairs)
/* namespace Axis {
  constexpr int kMass       =  0;
  constexpr int kPt         =  1;
  constexpr int kY          =  2;
  constexpr int kPVContr    =  3;
  constexpr int kOrigin     =  4;
  constexpr int kD0Type     =  5;
} */

// ─── Configuration ────────────────────────────────────────────────────────────
struct Config {
  std::string inputFile;
  std::string outputFile;
  std::string directory;
  std::string sparseName;
};

Config LoadConfig(const std::string& configPath = "config.json")
{
  std::ifstream cfgFile(configPath);
  if (!cfgFile.is_open()) {
    std::cerr << "[ERROR] Cannot open config file: " << configPath << "\n";
    throw std::runtime_error("Config file not found");
  }

  json cfg;
  cfgFile >> cfg;
  cfgFile.close();

  Config config;
  config.inputFile  = cfg.value("inputFile", "data/AnalysisResults_sigurd.root");
  config.outputFile = cfg.value("outputFile", "tests/projections_d0.root");
  config.directory  = cfg.value("directory", "hf-task-d0");
  config.sparseName = cfg.value("sparseName", "hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type");

  return config;
}

// pT bin edges [GeV/c]
const std::vector<double> kPtBins = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12., 16.};
//const std::vector<double> kPtBins = {0., 16.};

// Rapidity fiducial window
const double kYAbsMax = 0.5;

// ─── Helpers ──────────────────────────────────────────────────────────────────
void SetAxisRange(THnSparse* hs, int axis, double lo, double hi)
{
  TAxis* ax = hs->GetAxis(axis);
  int binLo = std::max(ax->FindBin(lo), 1);
  int binHi = std::min(ax->FindBin(hi), ax->GetNbins());
  ax->SetRange(binLo, binHi);
}

void ResetAllAxes(THnSparse* hs)
{
  for (int i = 0; i < hs->GetNdimensions(); ++i)
    hs->GetAxis(i)->SetRange(0, 0);
}

// ─── Main ─────────────────────────────────────────────────────────────────────
void project_sparse()
{
  // Load configuration from JSON
  Config config;
  try {
    config = LoadConfig("config.json");
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << "\n";
    return;
  }

  // ── Open input ──────────────────────────────────────────────────────────────
  TFile* fIn = TFile::Open(config.inputFile.c_str(), "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open " << config.inputFile << "\n";
    return;
  }

  TDirectory* dir = fIn->GetDirectory(config.directory.c_str());
  if (!dir) {
    std::cerr << "[ERROR] Directory '" << config.directory << "' not found.\n";
    fIn->Close();
    return;
  }

  THnSparse* hs = (THnSparse*)dir->Get(config.sparseName.c_str());
  if (!hs) {
    std::cerr << "[ERROR] THnSparse '" << config.sparseName << "' not found.\n";
    fIn->Close();
    return;
  }

  std::cout << "[INFO] Loaded THnSparse — "
            << hs->GetNdimensions() << " axes, "
            << hs->GetNbins() << " filled bins\n\n";

  // Print axis summary so you can cross-check
  std::cout << "── Axis summary ──────────────────────────────────────────\n";
  for (int i = 0; i < hs->GetNdimensions(); ++i) {
    TAxis* ax = hs->GetAxis(i);
    std::cout << "  Axis " << i << "  \"" << ax->GetTitle() << "\""
              << "  bins=" << ax->GetNbins()
              << "  [" << ax->GetXmin() << ", " << ax->GetXmax() << "]\n";
  }
  std::cout << "──────────────────────────────────────────────────────────\n\n";

  // ── Open output ─────────────────────────────────────────────────────────────
  TFile* fOut = TFile::Open(config.outputFile.c_str(), "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output files\n";
    fIn->Close();
    return;
  }

  const int nPt  = (int)kPtBins.size() - 1;
  int nSaved = 0;
  
  

  for (int iPt = 0; iPt < nPt; ++iPt) {
    double ptLo = kPtBins[iPt];
    double ptHi = kPtBins[iPt + 1];

    // Create a subdirectory per pT bin for tidy navigation from Python
    TDirectory* dPt = fOut->mkdir(
      Form("pt_%04.0f_%04.0f", ptLo * 100, ptHi * 100)  // e.g. "pt_0100_0200"
    );

    ResetAllAxes(hs);
    SetAxisRange(hs, Axis::kPt,     ptLo,     ptHi);
    SetAxisRange(hs, Axis::kY,     -kYAbsMax, kYAbsMax);

    TH1D* hMass = (TH1D*)hs->Projection(Axis::kMass);

    dPt->cd();
    hMass->Write();

    delete hMass;
    ++nSaved;

    std::cout << "  Saved  pT [" << ptLo << "," << ptHi << "]" << "\n";
  }

  ResetAllAxes(hs);
  fOut->Close();
  fIn->Close();

  std::cout << "\n[INFO] Done. Saved " << nSaved
            << " histograms to " << config.outputFile << "\n";
}