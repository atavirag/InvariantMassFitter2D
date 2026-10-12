#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TBox.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

bool isStatLinear = true;

void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels,
                  std::string outname, TFile* fout, bool isOS);
void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst);

// Format a double the same way runFitter does for filenames
std::string fmtMass(double v) {
    std::string s = std::to_string(v);
    return s.substr(0, 4);  // e.g. "1.71", "2.05"
}
std::string fmtPt(double v) {
    return std::to_string((int)v);  // e.g. "0", "2", "14"
}
std::string fmtEta(double v) {
    return std::to_string(v).substr(0, 3);  // e.g. "0.7"
}

void computeRawYieldSys() {

    // --- Configuration (mirrors fitter_config.json) ---
    std::string output_dir = "~/MyMacros/InvariantMassFitter/data/y_0.7_test/raw_yield_sys";
    std::string outfile     = "~/MyMacros/InvariantMassFitter/data/y_0.7_test/rawYieldSystematics.root";
    double eta = 0.7;

    //std::vector<double> pair_pt_bins = {0., 100.};
    std::vector<double> pair_pt_bins = {0., 2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 14.0, 50.};

    std::vector<double> mass_mins = {1.72, 1.73, 1.74, 1.75};
    std::vector<double> mass_maxs = {2.02, 2.03, 2.04, 2.05};

    std::vector<std::string> bkgFuncs = {"exponential",
                                         "chebychev",
                                         "polynomial",
                                         "poly3"};

    // Bad fits to skip (same logic as before)
    auto isBadFit = [](const std::string& bkg, double mmin, double mmax) -> bool {
        //if (bkg == "polynomial" && mmin == 1.75 && mmax == 2.02) return true;
        //if (bkg == "polynomial" && mmin == 1.72 && mmax == 2.04) return true;
        // Add more known-bad combinations here
        return false;
    };

    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");

    // Loop over pt bins
    for (size_t iPt = 0; iPt + 1 < pair_pt_bins.size(); iPt++) {
        double pt_min = pair_pt_bins[iPt];
        double pt_max = pair_pt_bins[iPt + 1];
        std::string pt_label = fmtPt(pt_min) + "_" + fmtPt(pt_max);

        std::vector<std::string> sysLabelsOS, sysLabelsLS;
        std::vector<std::string> filenamesOS, filenamesLS;

        for (const auto& bkg : bkgFuncs) {
            for (double mmin : mass_mins) {
                for (double mmax : mass_maxs) {
                    if (isBadFit(bkg, mmin, mmax)) continue;

                    std::string mass_label = fmtMass(mmin) + "_" + fmtMass(mmax);
                    std::string eta_label  = fmtEta(eta);
                    std::string bin_label  = "mass_" + mass_label + "_eta_" + eta_label + "_pairPt_" + pt_label;

                    // Corrected (efficiency) fit results
                    std::string fnOS = output_dir + "/fit_results_" + bkg + "_corrected_" + bin_label + ".root";
                    std::string fnLS = output_dir + "/fit_results_" + bkg + "_corrected_" + bin_label + ".root"; // same file, OS/LS inside

                    std::string label = bkg + " (" + fmtMass(mmin) + ", " + fmtMass(mmax) + ")";

                    filenamesOS.push_back(fnOS);
                    filenamesLS.push_back(fnLS);
                    sysLabelsOS.push_back(label);
                    sysLabelsLS.push_back(label);
                }
            }
        }

        std::string outnameOS = "RawYieldSys_OS_" + pt_label;
        std::string outnameLS = "RawYieldSys_LS_" + pt_label;

        calculateSys(filenamesOS, sysLabelsOS, outnameOS, fout, true);
        calculateSys(filenamesLS, sysLabelsLS, outnameLS, fout, false);
    }

    fout->Close();
}

void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels,
                  std::string outname, TFile* fout, bool isOS) {

    double average = 0., averageStat = 0., averageSyst = 0.;

    std::vector<double> rawYields, rawYieldsErr, chi2_ndfs, widths, widthsErr;

    std::string yieldHname  = isOS ? "hYields_os"      : "hYields_ls";
    std::string chi2Hname   = isOS ? "hOS_chi2ndf_m1"  : "hLS_chi2ndf_m1";
    std::string sigmaHname  = isOS ? "hOS_m1Sigma"     : "hLS_m1Sigma";

    for (const auto& filename : filenames) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << filename << std::endl;
            continue;
        }

        TH1D* hYields = dynamic_cast<TH1D*>(file->Get(yieldHname.c_str()));
        if (!hYields) {
            std::cerr << "Error: '" << yieldHname << "' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1D* hChi2 = dynamic_cast<TH1D*>(file->Get(chi2Hname.c_str()));
        if (!hChi2) {
            std::cerr << "Error: '" << chi2Hname << "' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1D* hSigma = dynamic_cast<TH1D*>(file->Get(sigmaHname.c_str()));
        if (!hSigma) {
            std::cerr << "Error: '" << sigmaHname << "' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        int binIndex = 1; // SgnSgn yield
        double yield = hYields->GetBinContent(binIndex);
        double yieldErr = hYields->GetBinError(binIndex);
        double chi2 = hChi2->GetBinContent(1);
        double sigma = hSigma->GetBinContent(1);
        double sigmaErr = hSigma->GetBinError(1);

        // --- Quality cuts ---
        bool passQuality = true;
        std::string skipReason = "";

        // Skip if relative error too large
        if (yield > 0 && yieldErr/yield > 0.5) {
            passQuality = false;
            skipReason = "large stat error";
        }
        // Skip if chi2/ndf too bad
        if (chi2 > 3.0 || chi2 < 0.3) {
            passQuality = false;
            skipReason = "bad chi2/ndf";
        }
        // Skip if yield is negative or zero
        if (yield <= 0) {
            passQuality = false;
            skipReason = "non-positive yield";
        }

        if (!passQuality) {
            std::cout << "  SKIPPED (" << skipReason << "): yield=" << yield
                      << " ± " << yieldErr << ", chi2/ndf=" << chi2 << std::endl;
            file->Close();
            continue;
        }

        rawYields.push_back(yield);
        rawYieldsErr.push_back(yieldErr);
        chi2_ndfs.push_back(chi2);
        widths.push_back(sigma);
        widthsErr.push_back(sigmaErr);

        std::cout << "  ACCEPTED: yield=" << yield << " ± " << yieldErr
                  << " (" << 100.*yieldErr/yield << "%)"
                  << ", chi2/ndf=" << chi2
                  << ", sigma=" << sigma << std::endl;

        file->Close();
    }

    if (rawYields.empty()) {
        std::cerr << "No valid fits found for " << outname << std::endl;
        return;
    }

    // --- Outlier rejection: remove yields far from median ---
    std::vector<double> sortedYields = rawYields;
    std::sort(sortedYields.begin(), sortedYields.end());
    double median = (sortedYields.size() % 2 == 0) 
                    ? (sortedYields[sortedYields.size()/2 - 1] + sortedYields[sortedYields.size()/2]) / 2.0
                    : sortedYields[sortedYields.size()/2];

    // Compute MAD (median absolute deviation)
    std::vector<double> absDeviations;
    for (double y : rawYields) {
        absDeviations.push_back(std::abs(y - median));
    }
    std::sort(absDeviations.begin(), absDeviations.end());
    double mad = (absDeviations.size() % 2 == 0)
                 ? (absDeviations[absDeviations.size()/2 - 1] + absDeviations[absDeviations.size()/2]) / 2.0
                 : absDeviations[absDeviations.size()/2];

    // Reject outliers: |yield - median| > 3 * MAD
    // (MAD is more robust than RMS for outlier detection)
    double outlierThreshold = 3.0; // adjust as needed (3.0 = conservative, 2.5 = stricter)
    
    std::vector<double> filteredYields, filteredYieldsErr, filteredChi2, filteredWidths, filteredWidthsErr;
    std::vector<std::string> filteredLabels;

    for (size_t i = 0; i < rawYields.size(); i++) {
        double deviation = std::abs(rawYields[i] - median);
        bool isOutlier = (deviation > outlierThreshold * mad);

        if (isOutlier) {
            std::cout << "  OUTLIER REJECTED: " << sysLabels[i]
                      << " yield=" << rawYields[i]
                      << " (|dev from median| = " << deviation
                      << " > " << outlierThreshold << "*MAD = " << outlierThreshold*mad << ")" << std::endl;
            continue;
        }

        filteredYields.push_back(rawYields[i]);
        filteredYieldsErr.push_back(rawYieldsErr[i]);
        filteredChi2.push_back(chi2_ndfs[i]);
        filteredWidths.push_back(widths[i]);
        filteredWidthsErr.push_back(widthsErr[i]);
        filteredLabels.push_back(sysLabels[i]);
    }

    std::cout << "Kept " << filteredYields.size() << "/" << rawYields.size() 
              << " fits after outlier rejection (median=" << median 
              << ", MAD=" << mad << ")" << std::endl;

    if (filteredYields.empty()) {
        std::cerr << "No fits survived quality cuts for " << outname << std::endl;
        return;
    }

    // Use filtered vectors from here on
    rawYields = filteredYields;
    rawYieldsErr = filteredYieldsErr;
    chi2_ndfs = filteredChi2;
    widths = filteredWidths;
    widthsErr = filteredWidthsErr;
    sysLabels = filteredLabels;

    // --- Fill histograms ---
    int n = rawYields.size();
    TH1D* hRawYieldSys = new TH1D("hRawYieldSys", "Raw Yield Systematics",    n, 0, n);
    TH1D* hWidthSys    = new TH1D("hWidthSys",    "Widths",                   n, 0, n);
    TH1D* hChi2Sys     = new TH1D("hChi2Sys",     "#chi^{2}/ndf",             n, 0, n);

    for (int i = 0; i < n; i++) {
        hRawYieldSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        hWidthSys   ->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        hChi2Sys    ->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());

        if (rawYieldsErr[i] == 0) continue;
        hRawYieldSys->SetBinContent(i+1, rawYields[i]);
        hRawYieldSys->SetBinError  (i+1, rawYieldsErr[i]);
        hWidthSys   ->SetBinContent(i+1, widths[i]);
        hWidthSys   ->SetBinError  (i+1, widthsErr[i]);
        hChi2Sys    ->SetBinContent(i+1, chi2_ndfs[i]);
    }

    GetStatistics(hRawYieldSys, average, averageStat, averageSyst);
    std::cout << outname << " | Average: " << average
              << "  stat: " << averageStat
              << "  syst: " << averageSyst << std::endl;

    // --- Styling ---
    hRawYieldSys->SetMarkerStyle(21);
    hRawYieldSys->SetMarkerColor(kBlue-4);
    hRawYieldSys->SetMarkerSize(1.2);
    hRawYieldSys->SetTitle("");
    hRawYieldSys->SetLineColor(kBlue-4);
    hRawYieldSys->GetYaxis()->SetTitle("Efficiency-corrected Yield");
    hRawYieldSys->GetYaxis()->SetTitleSize(0.045);
    hRawYieldSys->GetYaxis()->SetLabelSize(0.04);
    hRawYieldSys->GetYaxis()->SetRangeUser(average/3., average*5./3.);

    TPaveText* pave = new TPaveText(0.15, 0.2, 0.8, 0.45, "NDC");
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    pave->SetFillStyle(0);
    pave->SetBorderSize(0);
    pave->SetTextAlign(11);
    pave->AddText(Form("#it{N}_{DD}: %.0f #pm %.0f (stat., %.2f%%) #pm %.0f (syst., %.2f%%)",
                       average, averageStat, (averageStat/average)*100,
                       averageSyst, (averageSyst/average)*100));

    TPaveText* pt = new TPaveText(0.15, 0.72, 0.4, 0.88, "brNDC");
    pt->SetTextFont(42); pt->SetTextSize(0.04);
    pt->SetBorderSize(0); pt->SetFillStyle(0); pt->SetTextAlign(11);
    pt->AddText("#scale[1.35]{ALICE Work in progress}");
    pt->AddText("#scale[1.05]{pp #sqrt{s} = 13.6 TeV}");

    TPaveText* pt_left = new TPaveText(0.55, 0.72, 0.87, 0.88, "brNDC");
    pt_left->SetTextFont(42); pt_left->SetTextSize(0.04);
    pt_left->SetBorderSize(0); pt_left->SetFillStyle(0); pt_left->SetTextAlign(11);
    pt_left->AddText("#scale[1.05]{D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.}");
    if (isOS) pt_left->AddText("#scale[1.05]{D^{0}#bar{D}^{0} + #bar{D}^{0}D^{0} pairs}");
    else       pt_left->AddText("#scale[1.05]{D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0} pairs}");

    TLine* line       = new TLine(0, average,            n, average);
    TLine* statsUpper = new TLine(0, average+averageStat, n, average+averageStat);
    TLine* statsLower = new TLine(0, average-averageStat, n, average-averageStat);
    line->SetLineWidth(2); line->SetLineColor(kRed);
    statsUpper->SetLineWidth(2); statsUpper->SetLineColor(kGray+1); statsUpper->SetLineStyle(kDashed);
    statsLower->SetLineWidth(2); statsLower->SetLineColor(kGray+1); statsLower->SetLineStyle(kDashed);

    TBox* colorBox = new TBox(0.1, average-averageSyst, n-0.1, average+averageSyst);
    colorBox->SetFillColor(kRed-10);

    // --- Yield canvas ---
    TCanvas* c = new TCanvas(("c_" + outname).c_str(), "Yield Systematics", 1200, 600);
    c->SetTopMargin(0.04); c->SetLeftMargin(0.1);
    c->SetRightMargin(0.03); c->SetBottomMargin(0.25);
    c->SetTickx(); c->SetTicky();
    gStyle->SetOptStat(0);

    hRawYieldSys->LabelsOption("v", "X");
    hRawYieldSys->Draw("PE");
    colorBox->Draw("same");
    line->Draw("same");
    statsLower->Draw("same");
    statsUpper->Draw("same");
    pave->Draw("same");
    hRawYieldSys->Draw("samePE");
    pt->Draw("same");
    pt_left->Draw("same");

    // --- Chi2 canvas ---
    TCanvas* c_chi2 = new TCanvas(("c_chi2_" + outname).c_str(), "#chi^{2}/ndf", 1200, 400);
    hChi2Sys->SetMarkerStyle(21); hChi2Sys->SetMarkerSize(1.5);
    hChi2Sys->GetYaxis()->SetRangeUser(hChi2Sys->GetMinimum()-0.6, hChi2Sys->GetMaximum()+0.6);
    hChi2Sys->Draw("P");

    // --- Width canvas ---
    TCanvas* c_width = new TCanvas(("c_width_" + outname).c_str(), "Widths", 1200, 600);
    c_width->SetTopMargin(0.04); c_width->SetLeftMargin(0.1);
    c_width->SetRightMargin(0.03); c_width->SetBottomMargin(0.25);
    c_width->SetTickx(); c_width->SetTicky();

    hWidthSys->SetMarkerStyle(21); hWidthSys->SetMarkerSize(1.3);
    hWidthSys->SetTitle("");
    hWidthSys->GetYaxis()->SetTitle("Width (GeV/#it{c})");
    hWidthSys->GetYaxis()->SetTitleSize(0.045); hWidthSys->GetYaxis()->SetLabelSize(0.04);
    hWidthSys->GetXaxis()->SetTitleSize(0.045); hWidthSys->GetXaxis()->SetLabelSize(0.04);
    hWidthSys->GetYaxis()->SetRangeUser(hWidthSys->GetBinContent(1)-0.003, hWidthSys->GetBinContent(1)+0.003);
    hWidthSys->SetLineColor(kBlack);
    hWidthSys->LabelsOption("v", "X");
    hWidthSys->Draw("PE");
    pt->Draw("same");
    pt_left->Draw("same");

    // --- Save to file ---
    fout->cd();
    c->Write();
    c_chi2->Write();
    c_width->Write();

    c->Close();
    c_chi2->Close();
    c_width->Close();

    delete hRawYieldSys;
    delete hWidthSys;
    delete hChi2Sys;
}

void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst) {
    double av = 0., avSt = 0., rms = 0.;
    int nFits = 0;

    for (int i = 1; i <= histo->GetNbinsX(); i++) {
        double value = histo->GetBinContent(i);
        double unc   = histo->GetBinError(i);
        if (value > 0) {
            av += value;
            avSt += isStatLinear ? unc : unc*unc;
            nFits++;
        }
    }
    av   /= nFits;
    avSt  = isStatLinear ? avSt/nFits : sqrt(avSt)/nFits;

    for (int i = 1; i <= histo->GetNbinsX(); i++) {
        double value = histo->GetBinContent(i);
        if (value > 0) rms += (value - av)*(value - av);
    }
    rms = sqrt(rms / (nFits - 1));

    average     = av;
    averageStat = avSt;
    averageSyst = rms;
}