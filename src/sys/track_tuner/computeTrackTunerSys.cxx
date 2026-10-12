#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <string>

void computeTrackTunerSys() {
    TFile* fmain = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/track_tuner/AccEffD0ToKPi.root", "READ");
    if (!fmain || fmain->IsZombie()) {
        std::cerr << "Error: Could not open fmain" << std::endl;
        return;
    }


    TFile* fworse1p = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/track_tuner/AccEffD0ToKPi_worse_reso.root", "READ");
    if (!fworse1p || fworse1p->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p" << std::endl;
        return;
    }
    TFile* fworse1p5 = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/track_tuner/AccEffD0ToKPi_worse_reso_15.root", "READ");
    if (!fworse1p5 || fworse1p5->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p5" << std::endl;
        return;
    }


    TH1F* hEff_main = dynamic_cast<TH1F*>(fmain->Get("hAccEffD0ToKPiAll"));
    if (!hEff_main) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguousAll' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }

    TH1F* hEff_1p = dynamic_cast<TH1F*>(fworse1p->Get("hAccEffD0ToKPi_worse_resoAll"));
    if (!hEff_1p) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_worse_1pAll' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    TH1F* hEff_1p5 = dynamic_cast<TH1F*>(fworse1p5->Get("hAccEffD0ToKPi_worse_reso_15All"));
    if (!hEff_1p5) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_reso1p5All' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }

    // Create ratio histograms
    TH1F *hRatio_1p = (TH1F*)hEff_1p->Clone("hRatio_1p");
    hRatio_1p->Divide(hEff_main);

    TH1F *hRatio_1p5 = (TH1F*)hEff_1p5->Clone("hRatio_1p5");
    hRatio_1p5->Divide(hEff_main);

    // Create a canvas
    TCanvas *c = new TCanvas("c1", "Efficiency Ratios", 800, 400);
    c->SetTopMargin(0.04);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.03);
    c->SetBottomMargin(0.15);

    // Add ticks on all sides
    c->SetTickx();
    c->SetTicky();

    gStyle->SetOptStat(0);

    // Add a grid
    c->SetGridx();
    c->SetGridy();

    gStyle->SetGridColor(kGray + 2);
    gStyle->SetGridStyle(3);  // Dotted
    gStyle->SetGridWidth(1);

    // Customize the ratio histograms
    hRatio_1p->SetLineColor(kRed);
    hRatio_1p->SetMarkerColor(kRed);
    hRatio_1p->SetMarkerStyle(20);
    hRatio_1p->SetLineWidth(3);
    hRatio_1p->SetMarkerSize(1.2);
    hRatio_1p->SetTitle("; #it{p}_{T} (GeV/#it{c}); Efficiency ratios");
    hRatio_1p->GetYaxis()->SetTitleSize(0.05);
    hRatio_1p->GetYaxis()->SetLabelSize(0.045);
    hRatio_1p->GetXaxis()->SetTitleSize(0.05);
    hRatio_1p->GetXaxis()->SetLabelSize(0.05);
    hRatio_1p->GetYaxis()->SetTitleOffset(1);
    hRatio_1p->GetYaxis()->SetRangeUser(0.95, 1.05);  // Adjust range if needed

    hRatio_1p5->SetLineColor(kBlue-4);
    hRatio_1p5->SetMarkerColor(kBlue-4);
    hRatio_1p5->SetLineWidth(3);
    hRatio_1p5->SetMarkerStyle(21);
    hRatio_1p5->SetMarkerSize(1.2);

    // Draw histograms
    hRatio_1p->Draw("E P");
    

    TPaveText *pt = new TPaveText(0.15,0.68,0.5,0.88,"brNDC");
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);

    pt->AddText("#scale[1.35]{This Thesis}");
    pt->AddText("#scale[1.05]{pp #sqrt{s} = 13.6 TeV}");
    pt->AddText("#scale[1.05]{D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.}");

    pt->Draw("same");

    // Draw a horizontal line at 1 for reference
    TLine* line = new TLine(hRatio_1p->GetXaxis()->GetXmin(), 1., hRatio_1p->GetXaxis()->GetXmax(), 1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9);
    line->SetLineWidth(3);


    line->Draw("same");

    hRatio_1p->Draw("same E P");
    hRatio_1p5->Draw("E P SAME");

    // Add legend
    TLegend *leg = new TLegend(0.6, 0.25, 0.9, 0.38);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->AddEntry(hRatio_1p, "10% smearing", "lp");
    leg->AddEntry(hRatio_1p5, "15% smearing", "lp");

    leg->Draw();

    // Save the plot
    c->SaveAs("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/track_tuner/efficiency_ratios.pdf");

    c->Close();
    // Close files
    fmain->Close();
    fworse1p->Close();
    fworse1p5->Close();
}

