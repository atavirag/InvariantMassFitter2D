/// Quick macro to calculate the luminosity of the triggered dataset

#include <iostream>
#include <functional>

#include <cmath>
#include <iostream>
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"


using std::cout;
using std::endl;

bool verbose = true;
bool debug = false;

struct HistogramValues {
    double analysedTriggersValue;
    double inspectedTVXValue;
    double scalerValue;
};

void computeTriggerEfficiency() {

    TFile *fileMassMB = TFile::Open("~/MyMacros/Double-Dzero-correlations_copy/trigger_efficiency/rawYields_MB.root", "read");
    TFile *fileMassTrigger = TFile::Open("~/MyMacros/Double-Dzero-correlations_copy/trigger_efficiency/rawYields_singleCharm.root", "read");
    TFile *fileAnaResMB = TFile::Open("~/MyMacros/Double-Dzero-correlations_copy/trigger_efficiency/AnalysisResults/AnalysisResults_2024_sampled_triggered.root", "read");
    TFile *fileAnaResTrigger = TFile::Open("~/MyMacros/Double-Dzero-correlations_copy/trigger_efficiency/AnalysisResults/AnalysisResults_singleCharm_triggered.root", "read");

    TFile *fout = TFile::Open("~/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/trigger_efficiency/triggerEfficiency_singleCharm.root", "RECREATE");

    if (!fileMassMB || !fileMassTrigger) {
        cout << ">> ERROR Files with inv. mass fits not well readout" << endl;
        return;
    }
    if (!fileAnaResMB || !fileAnaResTrigger) {
        cout << ">> ERROR AnalysisResults files not well readout" << endl;
        return;
    }

    TDirectory *dir = dynamic_cast<TDirectory*>(fileAnaResTrigger->Get("hf-candidate-creator-2prong/Zorro"));
    // Map to store histogram values for each directory
    std::map<std::string, HistogramValues> directoryValues;

    // Recursive function to traverse directories and process histograms
    std::function<void(TDirectory*, const std::string&)> processDir =
        [&](TDirectory* dir, const std::string& path) {
            TIter nextKey(dir->GetListOfKeys());
            TKey* key;

            HistogramValues values = {-1, -1, -1}; // Default invalid values

            while ((key = (TKey*)nextKey())) {
                TObject* obj = key->ReadObj();

                if (obj->InheritsFrom("TDirectory")) {
                    // If it's a directory, recurse into it
                    TDirectory* subDir = (TDirectory*)obj;
                    std::string newPath = path + "/" + subDir->GetName();
                    processDir(subDir, newPath);
                } else if (obj->InheritsFrom("TH1")) {
                    // If it's a histogram, process it
                    TH1* hist = (TH1*)obj;

                    if (std::string(hist->GetName()) == "AnalysedTriggers") {
                        int binAnalysedTriggers = hist->GetXaxis()->FindBin("fHfSingleCharm2P");
                        values.analysedTriggersValue = hist->GetBinContent(binAnalysedTriggers);
                    } else if (std::string(hist->GetName()) == "InspectedTVX") {
                        values.inspectedTVXValue = hist->Integral();
                    } else if (std::string(hist->GetName()) == "Scalers") {
                        //int bin = hist->GetNbinsX(); // the last bin contains the triggered events
                        int bin = hist->GetXaxis()->FindBin("fHfSingleCharm2P"); // the last bin contains the triggered events
                        values.scalerValue = hist->GetBinContent(bin);
                    }
                }
            }

            // If any values were found, store them
            if (values.analysedTriggersValue != -1 || values.inspectedTVXValue != -1 || values.scalerValue != -1) {
                directoryValues[path] = values;
            }
        };

    // Start processing from the root directory
    processDir(dir, "");

    // Get normalisation of triggered data
    double totAnalysedTriggers = 0.;
    double totScalers = 0.;
    double totInspectedTVX = 0.;
    double totNevTrigger = 0.;
    for (const auto& [dirPath, values] : directoryValues) {

        double nEvents = (values.inspectedTVXValue * values.analysedTriggersValue) / values.scalerValue;
        //totNev += nEvents;

        totAnalysedTriggers += values.analysedTriggersValue;
        totScalers += values.scalerValue;
        totInspectedTVX += values.inspectedTVXValue;

        if (debug) {
            std::cout << "Directory: " << dirPath << std::endl;
            std::cout << "  AnalysedTriggersOfInterest: " << values.analysedTriggersValue << std::endl;
            std::cout << "  InspectedTVX: " << values.inspectedTVXValue << std::endl;
            std::cout << "  Scalers: " << values.scalerValue << std::endl;
            std::cout << "  Number of events: " << nEvents << std::endl;
        }
    }

    totNevTrigger = totInspectedTVX * totAnalysedTriggers / totScalers;

    if (verbose) {
        std::cout << "  NORMALISATION OF TRIGGERED DATA" << std::endl;
        std::cout << "  Total AnalysedTriggersOfInterest: " << totAnalysedTriggers << std::endl;
        std::cout << "  Total InspectedTVX: " << totInspectedTVX << std::endl;
        std::cout << "  Total Scalers: " << totScalers << std::endl;
        std::cout << "  Number of events: " << totNevTrigger << std::endl;
    }

    // Normalisation of MB data
    TDirectory *dirMB = dynamic_cast<TDirectory*>(fileAnaResMB->Get("bc-selection-task"));
    TH1F *hCounterTVX = (TH1F *)dirMB->Get("hCounterTVX");
    double countsTVX = hCounterTVX->GetEntries();

    if (verbose) {
        std::cout << "\n  NORMALISATION OF MB DATA" << std::endl;
        std::cout << "  Total TVX entries: " << countsTVX << std::endl;
    }
    double ratioTriggerMB = totNevTrigger / countsTVX;
    std::cout << "\n trigger / MB = " << totNevTrigger / countsTVX << endl; // Fabrizio dijo que el singleCharm tenía un downscale muy grande

    // Divide normalised raw yields
    TH1F *hRawYieldMB = (TH1F *)fileMassMB->Get("hRawYields");
    TH1F *hRawYieldTrigger = (TH1F *)fileMassTrigger->Get("hRawYields");

    if (!hRawYieldMB || !hRawYieldTrigger) {
        cerr << "ERROR: raw yields not found" << endl;
        return;
    }

    TH1F *hTriggerEfficiency = (TH1F *)hRawYieldMB->Clone();

    for (int i = 1; i <= hRawYieldMB->GetNbinsX(); i++) {
        double triggerEff = (hRawYieldTrigger->GetBinContent(i) / hRawYieldMB->GetBinContent(i)) / ratioTriggerMB;
        double triggerEffError = triggerEff * sqrt(pow(hRawYieldTrigger->GetBinError(i) / hRawYieldTrigger->GetBinContent(i), 2) +
                                      pow(hRawYieldMB->GetBinError(i)/ hRawYieldMB->GetBinContent(i), 2));
        cout << triggerEff<< endl;
        hTriggerEfficiency->SetBinContent(i, triggerEff);
        hTriggerEfficiency->SetBinError(i, triggerEffError);
    }

    // Fit with a pol0 (constant function)
    TH1F *hTriggerEfficiency_clone = (TH1F *)hTriggerEfficiency->Clone();
    TF1 *fitPol0 = new TF1("fitPol0", "pol0", hTriggerEfficiency->GetXaxis()->GetXmin(), hTriggerEfficiency->GetXaxis()->GetXmax());
    hTriggerEfficiency_clone->Fit(fitPol0, "Q");  // "Q" suppresses printout, remove if you want to see fit info

    // Optional: Print fit result
    std::cout << "Fit result (constant value): " << fitPol0->GetParameter(0) << " ± " << fitPol0->GetParError(0) << std::endl;

    fout->cd();
    TCanvas *canvas = new TCanvas("c", "c", 800, 600);
    canvas->SetGrid();  // Add grid
    // Set margins
    canvas->SetTopMargin(0.04);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.03);
    canvas->SetBottomMargin(0.12);

    // Remove stats box globally
    gStyle->SetOptStat(0);

    // Add ticks on all sides
    canvas->SetTickx();
    canvas->SetTicky();

    // Enable grid lines
    canvas->SetGridx();
    canvas->SetGridy();

    // Grid style customization
    gStyle->SetGridColor(kGray + 2);   // Light gray
    gStyle->SetGridStyle(3);           // Dotted lines
    gStyle->SetGridWidth(1);           // Thin lines

    // Apply style settings immediately (optional but recommended)
    gROOT->ForceStyle();
    hTriggerEfficiency->SetTitle("");

    hTriggerEfficiency->SetLineColor(kRed-4);
    hTriggerEfficiency->SetLineWidth(2);
    hTriggerEfficiency->SetMarkerStyle(20);
    hTriggerEfficiency->SetMarkerSize(1);
    hTriggerEfficiency->SetMarkerColor(kRed-4);

    hTriggerEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hTriggerEfficiency->GetYaxis()->SetTitle("#it{N}^{trigger}_{corr}/#it{N}^{MB}_{corr}");
    hTriggerEfficiency->GetXaxis()->SetTitleSize(0.045);
    hTriggerEfficiency->GetYaxis()->SetTitleSize(0.045);
    hTriggerEfficiency->GetXaxis()->SetLabelSize(0.04);
    hTriggerEfficiency->GetYaxis()->SetLabelSize(0.04);


    hTriggerEfficiency->GetYaxis()->SetRangeUser(0.6, 1.4); // Adjust Y range
    hTriggerEfficiency->GetXaxis()->SetRangeUser(1.0, 16.0); // Adjust Y range

    // Draw a red line at y=1
    //TLine *line = new TLine(hTriggerEfficiency->GetXaxis()->GetXmin(), 1, hTriggerEfficiency->GetXaxis()->GetXmax(), 1);
    //line->SetLineColor(kRed);
    //line->SetLineWidth(2);
    //line->SetLineStyle(2);  // Dashed line
    fitPol0->SetLineColor(kAzure-3);
    fitPol0->SetLineWidth(2);
    fitPol0->SetLineStyle(9);

    hTriggerEfficiency->Draw();
    fitPol0->Draw("same");

    // Get fit parameters
    double constant = fitPol0->GetParameter(0);
    double constantError = fitPol0->GetParError(0);

    TLegend* legend = new TLegend(0.45, 0.2, 0.9, 0.45); // Adjust coordinates as needed

    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hTriggerEfficiency, "Data", "ep");
    legend->AddEntry(fitPol0, Form("Fit: pol0 = %.4f #pm %.4f", constant, constantError), "l");
    legend->Draw();

    TPaveText *pt = new TPaveText(0.15,0.72,0.4,0.85,"brNDC");
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);

    pt->AddText("#scale[1.35]{This Thesis}");
    pt->AddText("#scale[1.05]{pp #sqrt{s} = 13.6 TeV}");
    pt->AddText("#scale[1.05]{D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.}");
    pt->Draw("same");

    canvas->SaveAs("~/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/trigger_efficiency/trigger_eff.pdf");
    hTriggerEfficiency->Write();
    canvas->Write();
}
