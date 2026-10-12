#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <iostream>
#include <vector>
#include <string>

bool isStatLinear=true; 
void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname, TFile* fout);
void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst);

vector<double> defaultCuts = {0.50, 0.60, 0.55, 0.55, 0.4, 0.25};

void computeBDTSys() {

    std::vector<std::string> filenames = {};
    std::vector<std::string> sysLabels = {};

    int const nPos = 20;
    
    for (int pos = 0; pos < nPos; pos++) {
        //std::string filename = Form("systematics/bdt_sys/rawYield_bdtSys_pos%d.root", pos);
        std::string filename = Form("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/raw_yields/rawYield_bdtSys_pos%d.root", pos);
        std::string sysLabel = Form("%d", pos);

        filenames.push_back(filename);
        sysLabels.push_back(sysLabel);
    }
    TFile* fout = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/bdt_systematics_results.root", "RECREATE");
    calculateSys(filenames, sysLabels, "bdtSys", fout);
}

void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname, TFile* fout) {

    TFile* fEffs = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/efficiencies/efficiencies_bdt.root", "READ");

    TGraphErrors *gEff_01 = (TGraphErrors *)fEffs->Get("gEffPerCut_0_1");
    TGraphErrors *gEff_12 = (TGraphErrors *)fEffs->Get("gEffPerCut_1_2");
    TGraphErrors *gEff_23 = (TGraphErrors *)fEffs->Get("gEffPerCut_2_3");
    TGraphErrors *gEff_34 = (TGraphErrors *)fEffs->Get("gEffPerCut_3_4");
    TGraphErrors *gEff_45 = (TGraphErrors *)fEffs->Get("gEffPerCut_4_5");
    TGraphErrors *gEff_56 = (TGraphErrors *)fEffs->Get("gEffPerCut_5_6");

    std::vector<TGraphErrors *> gEffs = {gEff_01, gEff_12, gEff_23, gEff_34, gEff_45, gEff_56};
    std::vector<std::vector<double>> effs(gEffs.size()); // Initialize vector with correct size
    std::vector<std::vector<double>> effsErr(gEffs.size()); // Initialize vector with correct size

    for (size_t i = 0; i < gEffs.size(); ++i) {
        TGraphErrors *g = dynamic_cast<TGraphErrors *>(gEffs[i]); // Ensure it's a TGraphErrors
        if (!g) continue; // Check if the graph exists
    
        int nPoints = g->GetN();
        std::vector<double> cut(nPoints), eff(nPoints), err(nPoints);

        for (int j = 0; j < nPoints; ++j) {
            g->GetPoint(j, cut[j], eff[j]); // Extract each (x, y) point
            err[j] = g->GetErrorY(j);       // Extract error on y-axis
            effs[i].push_back(eff[j]);      // Store efficiency values
            effsErr[i].push_back(err[j]);   // Store efficiency errors
            cout << "Efficiency: " << eff[j] << " +- " << err[j] << endl;
        }
        double x_last, y_last;
gEffs.back()->GetPoint(nPoints - 1, x_last, y_last);
double err_last = gEffs.back()->GetErrorY(nPoints - 1);
std::cout << "Last point: x = " << x_last << ", y = " << y_last << " ± " << err_last << std::endl;
    }

    std::vector<double> averages = {}, averageStats = {}, averageSysts = {};

    std::vector<double> ptLims = {1., 2., 3., 5., 7., 10., 24.};

    std::vector<std::vector<double>> rawYields = {};
    std::vector<std::vector<double>> rawYieldsErr = {};
    std::vector<std::vector<double>> chi2_ndfs = {};
    std::vector<std::vector<double>> widths = {};
    std::vector<std::vector<double>> widthsErr = {};
    for (const auto& filename : filenames) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        cout << "In file " << filename.c_str() << endl;
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << filename << std::endl;
            continue;
        }

        TH1D* hYields = dynamic_cast<TH1D*>(file->Get("hRawYields"));
        if (!hYields) {
            std::cerr << "Error: Histogram 'hYields' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1D* hChi2_ndf = dynamic_cast<TH1D*>(file->Get("hRawYieldsChiSquare"));
        if (!hChi2_ndf) {
            std::cerr << "Error: Histogram 'hRawYieldsChiSquare' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1D* hWidth = dynamic_cast<TH1D*>(file->Get("hRawYieldsSigma"));
        if (!hWidth) {
            std::cerr << "Error: Histogram 'hSigma' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        std::vector<double> rawYieldsPt = {};
        std::vector<double> rawYieldsErrPt = {};
        std::vector<double> chi2_ndfsPt = {};
        std::vector<double> widthsPt = {};
        std::vector<double> widthsErrPt = {};
        for (int ipt = 0; ipt < ptLims.size(); ipt++) {
            double binContent = hYields->GetBinContent(ipt+1);
            rawYieldsPt.push_back(binContent);
            double binError = hYields->GetBinError(ipt+1);
            rawYieldsErrPt.push_back(binError);

            double chi2_ndf = hChi2_ndf->GetBinContent(ipt+1);
            chi2_ndfsPt.push_back(chi2_ndf);

            double width = hWidth->GetBinContent(ipt+1);
            widthsPt.push_back(width);
            double widthError = hWidth->GetBinError(ipt+1);
            widthsErrPt.push_back(widthError);

            cout << "binContent: " << binContent << ", chi2_ndf: " << chi2_ndf << ", width: " << width << endl;
            }
            rawYields.push_back(rawYieldsPt);
            rawYieldsErr.push_back(rawYieldsErrPt);
            chi2_ndfs.push_back(chi2_ndfsPt);
            widths.push_back(widthsPt);
            widthsErr.push_back(widthsErrPt);

            file->Close();
        }

        // plot systematics
        for (int ipt = 0; ipt < ptLims.size()-1; ipt++) {

            double average = 0., averageStat = 0., averageSyst = 0.;


            TH1D *hRawYieldSys = new TH1D(Form("hRawYieldSys_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), Form("BDT selection Systematic, pt (%.1f, %.1f) GeV/#it{c}", ptLims[ipt], ptLims[ipt+1]), rawYields.size(), 1, rawYields.size());
            TH1D *hWidthSys = new TH1D(Form("hWidthSys_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "Widths", rawYields.size(), 0, rawYields.size());
            TH1D *hChi2Sys = new TH1D(Form("hChi2Sys_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "#chi^{2}/ndf", rawYields.size(), 0, rawYields.size());
            double squaredSum = 0.;

            TH1D *hEfficiencies = new TH1D(Form("hEfficiencies_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), Form("BDT efficiencies, pt (%.1f, %.1f) GeV/#it{c}", ptLims[ipt], ptLims[ipt+1]), effs[ipt].size(), 1, effs[ipt].size());

            if (rawYields.size() != effs[ipt].size()) {
                cout << "ERROR: EFFICIENCIES AND YIELDS ARE NOT DONE WITH THE SAME NUMBER OF CUTS" << endl;
                return;
            }

            if (sysLabels.size() != rawYields.size()) {
                std::cerr << "Mismatch between systematics labels and raw yields size!" << std::endl;
            }
            for (int i = 0; i < rawYields.size(); i++) {
                if (rawYieldsErr[i][ipt] == 0) continue;
                hRawYieldSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());

                hWidthSys->GetXaxis()->SetRange(2, rawYields.size());
                hChi2Sys->GetXaxis()->SetRange(2, rawYields.size());
                hWidthSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
                hChi2Sys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());

                hRawYieldSys->SetBinContent(i+1, rawYields[i][ipt]);
                hRawYieldSys->SetBinError(i+1, rawYieldsErr[i][ipt]);

                hWidthSys->SetBinContent(i+1, widths[i][ipt]);
                hWidthSys->SetBinError(i+1, widthsErr[i][ipt]);

                hChi2Sys->SetBinContent(i+1, chi2_ndfs[i][ipt]);

                hEfficiencies->SetBinContent(i+1, effs[ipt][i]);
                hEfficiencies->SetBinError(i+1, effsErr[ipt][i]);
                hEfficiencies->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());

                //cout << "Raw yield error for bin " << i << ": " << rawYieldsErr[i][ipt] << endl;
        }

        // Divide hRawYieldSys by effs
        TCanvas *c_divide = new TCanvas(Form("c_divide_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "N_raw/effs", 800, 600);

        TH1D *hSys = (TH1D *)hRawYieldSys->Clone();
        hSys->SetNameTitle("hSys", "N_raw/efficiency");
        hSys->Sumw2();  // Calculate and track errors for h1
        hEfficiencies->Sumw2();

        // Loop over bins
        for (int bin = 1; bin <= hRawYieldSys->GetNbinsX(); bin++) {
            double yield = hRawYieldSys->GetBinContent(bin);
            double yieldError = hRawYieldSys->GetBinError(bin);

            double efficiency = hEfficiencies->GetBinContent(bin);
            double efficiencyError = hEfficiencies->GetBinError(bin);

            // Ensure no division by zero
            if (efficiency != 0)
            {
                // Manually calculate the result of division
                double result = yield / efficiency;

                // Propagate binomial error for the ratio
                double resultError = sqrt(pow(yieldError * efficiency, 2) + pow(yield * efficiencyError, 2)) / (efficiency * efficiency);

                // Set the result and error in the new histogram
                hSys->SetBinContent(bin, result);
                hSys->SetBinError(bin, resultError);
            } else {
                // If efficiency is 0, set the result to 0 and the error to a large value or 0
                hSys->SetBinContent(bin, 0);
                hSys->SetBinError(bin, 0);  // Or set some large value as error if needed
            }
        }

        // Calculate the ratio of trial/default
        TH1D *hSysRatio = (TH1D *)hSys->Clone();
        hSysRatio->SetNameTitle("hSysRatio", "Trial/Default");
        // Loop over bins
        for (int i = 1; i <= hSysRatio->GetNbinsX(); i++) {
            hSysRatio->SetBinContent(i, hSys->GetBinContent(i)/hSys->GetBinContent(1));
            double error = sqrt( pow(hSys->GetBinError(i) / hSys->GetBinContent(i), 2) + pow(hSys->GetBinError(1) / hSys->GetBinContent(1), 2) );
            hSysRatio->SetBinError(i, error);
        }

        GetStatistics(hSysRatio, average, averageStat, averageSyst);
        std::cout << "Average: " << average << ", average stat: " << averageStat << ", average syst: " << averageSyst << endl;

        averages.push_back(average);
        averageStats.push_back(averageStat);
        averageSysts.push_back(averageSyst);
        cout << "average " << average << " averageSyst " << averageSyst << endl;


        hSysRatio->GetXaxis()->SetRange(2, rawYields.size());
        hSys->GetXaxis()->SetRange(2, rawYields.size());

        TLine *lineRatio = new TLine(2, hSysRatio->GetBinContent(1), rawYields.size(), hSysRatio->GetBinContent(1));
        lineRatio->SetLineWidth(2);
        lineRatio->SetLineColor(kRed);

        TLine *lineSys = new TLine(2, hSys->GetBinContent(1), rawYields.size(), hSys->GetBinContent(1));
        lineSys->SetLineWidth(2);
        lineSys->SetLineColor(kRed);

        TLine *statsUpperSys = new TLine(2, hSys->GetBinContent(1) + hSys->GetBinError(1), rawYields.size(), hSys->GetBinContent(1) + hSys->GetBinError(1));
        statsUpperSys->SetLineWidth(2);
        statsUpperSys->SetLineColor(kGray+1);
        statsUpperSys->SetLineStyle(kDashed);

        TLine *statsLowerSys = new TLine(2, hSys->GetBinContent(1) - hSys->GetBinError(1), rawYields.size(), hSys->GetBinContent(1) - hSys->GetBinError(1));
        statsLowerSys->SetLineWidth(2);
        statsLowerSys->SetLineColor(kGray+1);
        statsLowerSys->SetLineStyle(kDashed);

        //hSys->Divide(hSys, hEfficiencies, 1, 1, "B");

        hSys->SetMarkerStyle(21);
        hSys->SetMarkerSize(1.3);
        hSys->SetLineColor(kBlack);
        hSys->Draw();
        lineSys->Draw("same");
        statsLowerSys->Draw("same");
        statsUpperSys->Draw("same");

        //c_divide->SaveAs(Form("yieldOverEff_bdtSys%.1f_%.1f.png", ptLims[ipt], ptLims[ipt+1]));

        // Make a canvas to plot the efficiencies
        hEfficiencies->GetXaxis()->SetRange(2, rawYields.size());
        hEfficiencies->SetMarkerStyle(21);
        hEfficiencies->SetMarkerSize(1.3);

        TCanvas *cEffs = new TCanvas(Form("c_effs_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "Efficiencies", 800, 600);

        TLine *lineEffs = new TLine(2, hEfficiencies->GetBinContent(1), rawYields.size(), hEfficiencies->GetBinContent(1));
        lineEffs->SetLineWidth(2);
        lineEffs->SetLineColor(kRed);

        TLine *statsUpperEffs = new TLine(2, hEfficiencies->GetBinContent(1) + hEfficiencies->GetBinError(1), rawYields.size(), hEfficiencies->GetBinContent(1) + hEfficiencies->GetBinError(1));
        statsUpperEffs->SetLineWidth(2);
        statsUpperEffs->SetLineColor(kGray+1);
        statsUpperEffs->SetLineStyle(kDashed);

        TLine *statsLowerEffs = new TLine(2, hEfficiencies->GetBinContent(1) - hEfficiencies->GetBinError(1), rawYields.size(), hEfficiencies->GetBinContent(1) - hEfficiencies->GetBinError(1));
        statsLowerEffs->SetLineWidth(2);
        statsLowerEffs->SetLineColor(kGray+1);
        statsLowerEffs->SetLineStyle(kDashed);

        hEfficiencies->Draw();
        lineEffs->Draw("same");
        statsUpperEffs->Draw("same");
        statsLowerEffs->Draw("same");
        //cEffs->SaveAs(Form("efficiencies_bdtSys%.1f_%.1f.png", ptLims[ipt], ptLims[ipt+1]));


        //hRawYieldSys->SetMarkerColor(kBlue);
        hRawYieldSys->SetMarkerStyle(21);
        //hRawYieldSys->SetLineColor(kBlue);
        //hRawYieldSys->GetYaxis()->SetRangeUser(0.6, 1.6);

        // Create TPaveText
        TPaveText *pave = new TPaveText(0.15, 0.7, 0.8, 0.85, "NDC");  // Normalized coordinates
        pave->SetFillColor(0);  // Transparent background
        pave->SetTextAlign(12); // Align left
        pave->SetTextSize(0.03);
        //pave->AddText(Form("#it{N}_{D}: %.0f #pm %.0f (stat., %.2f%%) #pm %.0f (syst., %.2f%%)", average, averageStat, (averageStat/average)*100, averageSyst, (averageSyst/average)*100));


        TLine *line = new TLine(2, rawYields[0][ipt], rawYields.size(), rawYields[0][ipt]);
        line->SetLineWidth(2);
        line->SetLineColor(kRed);

        double statErr = rawYieldsErr[0][ipt]/rawYields[0][ipt];
        //cout << rawYieldsErr[0][ipt]/rawYields[0][ipt] << endl;

        TLine *statsUpper = new TLine(2, rawYields[0][ipt] + rawYieldsErr[0][ipt], rawYields.size(), rawYields[0][ipt] + rawYieldsErr[0][ipt]);
        statsUpper->SetLineWidth(2);
        statsUpper->SetLineColor(kGray+1);
        statsUpper->SetLineStyle(kDashed);

        TLine *statsLower = new TLine(2, rawYields[0][ipt] - rawYieldsErr[0][ipt], rawYields.size(), rawYields[0][ipt] - rawYieldsErr[0][ipt]);
        statsLower->SetLineWidth(2);
        statsLower->SetLineColor(kGray+1);
        statsLower->SetLineStyle(kDashed);

        TCanvas *c = new TCanvas(Form("c_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "Raw Yield Systematics", 800, 600);

        gStyle->SetOptStat(0);

        TBox *colorBox = new TBox(0.01, average-averageSyst, rawYields.size()-1, average + averageSyst);  // (x1, y1, x2, y2)
        colorBox->SetFillColor(kRed-10);               // Set box color

        hRawYieldSys->GetXaxis()->SetRange(2, rawYields.size());
        hRawYieldSys->Draw("PE");
        //colorBox->Draw("same");
        line->Draw("same");
        statsLower->Draw("same");
        //pave->Draw("same");
        statsUpper->Draw("same");
        hRawYieldSys->Draw("samePE");
        //legend->Draw();

        //c->SaveAs(Form("rawYield_bdtSys%.1f_%.1f.png", ptLims[ipt], ptLims[ipt+1]));

        // Add plot with the chi2_ndf
        TCanvas *c_chi2 = new TCanvas(Form("c_chi2_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "#chi^{2}/ndf", 800, 400);

        hChi2Sys->SetMarkerStyle(21);
        hChi2Sys->SetMarkerSize(1.5);
        hChi2Sys->GetYaxis()->SetRangeUser(hChi2Sys->GetBinContent(1)-0.6, hChi2Sys->GetBinContent(1)+0.6);
        hChi2Sys->Draw("P");
        std::string outChiPng = "Chi2_over_ndf_" + outname + ".png";
        //c_chi2->SaveAs(outChiPng.c_str());

        // Add plot with the widths
        TCanvas *c_width = new TCanvas(Form("c_width_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "Widths", 800, 400);
        hWidthSys->SetMarkerStyle(21);
        hWidthSys->SetMarkerSize(1.3);
        hWidthSys->GetYaxis()->SetRangeUser(hWidthSys->GetBinContent(1)-0.003, hWidthSys->GetBinContent(1)+0.003);
        //hWidthSys->SetLineSize(2);
        hWidthSys->SetLineColor(kBlack);

        hWidthSys->Draw("PE");
        std::string outWidthPng = "Width_" + outname + ".png";
        //c_width->SaveAs(outWidthPng.c_str());

        TCanvas *c_ratio = new TCanvas(Form("c_ratio_%.1f_%.1f", ptLims[ipt], ptLims[ipt+1]), "Ratio of N_raw/effs", 800, 600);
        hSysRatio->GetYaxis()->SetRangeUser(0.9, 1.1);
        hSysRatio->SetMarkerStyle(21);
        hSysRatio->SetMarkerSize(1.3);
        hSysRatio->Draw();
        lineRatio->Draw("same");

        //c_ratio->SaveAs(Form("ratioSys_bdtSys%.1f_%.1f.png", ptLims[ipt], ptLims[ipt+1]));

        fout->cd();
        c->Write();
        c_chi2->Write();
        c_ratio->Write();
        c_divide->Write();
        c_width->Write();
        cEffs->Write();

        hRawYieldSys->Write();
        hEfficiencies->Write();
        hSysRatio->Write();
        
        hWidthSys->Write();
        hChi2Sys->Write();

        c->Close();
        c_chi2->Close();
        c_width->Close();

        delete hRawYieldSys;
        delete hWidthSys;
        delete hChi2Sys;
    }

    // Define bin centers and asymmetric errors for bin widths
    std::vector<double> ptCenters, ptErrLow, ptErrHigh, systOverAvg;
    for (size_t i = 0; i < ptLims.size() - 1; i++) {
        double center = 0.5 * (ptLims[i] + ptLims[i + 1]);
        ptCenters.push_back(center);
        ptErrLow.push_back(center - ptLims[i]);       // Distance to lower edge
        ptErrHigh.push_back(ptLims[i + 1] - center);  // Distance to upper edge

        // Compute systematic uncertainty / average
        if (averages[i] != 0) {
            systOverAvg.push_back(averageSysts[i] / averages[i]);
        } else {
            systOverAvg.push_back(0);  // Avoid division by zero
        }
    }

    // Create TGraphAsymmErrors to represent bin ranges
    TGraphAsymmErrors* gSystOverAvg = new TGraphAsymmErrors(ptCenters.size(), 
                                                            ptCenters.data(), 
                                                            systOverAvg.data(), 
                                                            ptErrLow.data(), 
                                                            ptErrHigh.data(), 
                                                            nullptr, nullptr);
    gSystOverAvg->SetTitle("Systematic Uncertainty / Average; p_{T} (GeV/c); Syst / Avg");
    gSystOverAvg->SetMarkerStyle(21);
    gSystOverAvg->SetMarkerSize(1.2);
    gSystOverAvg->SetMarkerColor(kRed);
    gSystOverAvg->SetLineColor(kRed);
    gSystOverAvg->SetLineWidth(2);

    // Draw the graph
    // Set y-axis range (modify these values as needed)
    double yMin = 0.0, yMax = 0.15;  // Example values, adjust based on your data
    gSystOverAvg->GetYaxis()->SetRangeUser(yMin, yMax);

TCanvas* c = new TCanvas("c", "Systematic Uncertainty / Average", 1200, 600);
c->SetTopMargin(0.04);
c->SetLeftMargin(0.1);
c->SetRightMargin(0.03);
c->SetBottomMargin(0.25);

// Add ticks on all sides
c->SetTickx();
c->SetTicky();

gStyle->SetOptStat(0);

gSystOverAvg->Draw("AP");

TPaveText *pt = new TPaveText(0.15,0.72,0.4,0.88,"brNDC");
pt->SetTextFont(42);
pt->SetTextSize(0.04);
pt->SetBorderSize(0);
pt->SetFillStyle(0);
pt->SetTextAlign(11);

pt->AddText("#scale[1.35]{This Thesis}");
pt->AddText("#scale[1.05]{pp #sqrt{s} = 13.6 TeV}");

TPaveText *pt_left = new TPaveText(0.55,0.72,0.87,0.88,"brNDC");
pt_left->SetTextFont(42);
pt_left->SetTextSize(0.04);
pt_left->SetBorderSize(0);
pt_left->SetFillStyle(0);
pt_left->SetTextAlign(11);
pt_left->AddText("#scale[1.05]{D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.}");

pt->Draw("same");
pt_left->Draw("same");

// Save the plot
c->SaveAs("SystOverAvg.png");
gSystOverAvg->Write();

}

void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst) {
    double av=0., avSt=0., rms=0.;
    int nbins = histo->GetNbinsX();
    int nFits=0;

    for (int i=1; i<nbins; i++){
        double value = histo->GetBinContent(i);
        double unc = histo->GetBinError(i);
        if(value>0){
            av += value;
            if(isStatLinear) {
                avSt += unc;
            } else{
                avSt += unc*unc;
            }
            nFits++;
        }
        cout <<" i="<<i<<" val="<<value;
    }
    cout<<endl;
    av /= nFits;
    if(isStatLinear) {
        avSt /= nFits;
    } else  avSt = TMath::Sqrt( avSt )/nFits;

    for (int i=1; i<nbins; i++){
        double value = histo->GetBinContent(i);
        if(value>0) rms += (value - av)*(value - av);
    }

    rms = TMath::Sqrt( rms /(nFits-1) );
    average=av;
    averageStat=avSt;
    averageSyst=rms;
}