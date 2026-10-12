#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TMath.h"
#include "TSystem.h"
#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

bool verbose, debug;
bool isLog;
// pT ranges
int const nPtBins = 6;
double ptBinsEdges[nPtBins+1];

void DrawHist(TH1F *h, TString saveName, TFile *f, TString yAxisLabel, bool isLog);
void DrawEffHists(TH1F *h, TH1F *hPID, TH1F *hCand, TH1F *hTopol, TH1F *hHFFlag, 
                    TString saveName, TFile *f, TString yAxisLabel, bool isLog);
void DrawHist(TH1F *h, TH1F *h2, TString saveName, TFile *f, TString yAxisLabel, bool isLog);

TH1F * GetRawYield(vector<TH1F *> vh, TString const& name, double ptBinsEdges);
TH1F * CalculateEfficiency(TH1F *hYieldGen, TH1F *hYieldReco, TString nameReco);
vector<TH1F *> GetYPerPtRange(int const& nbins, TH2F* hPtVsY, TString nameParticle);
vector<TH1F *> MakeAnalysis(TFile *fRecoPrompt, TFile *fRecoNonPrompt, TFile *fGenPrompt, TFile *fGenNonPrompt, TFile *outputFile,
                    TString const& dirnameReco, TString const& suffix, TString const& hPromptName, TString const& hNonPromptName, bool const& isLog);
void efficiencyFromSparse() {
    // Load JSON config
    ifstream jsonFile("/home/andrea/MyMacros/Double-Dzero-correlations/.config/config-efficiencyFromSparse.json");
    if (!jsonFile.is_open()) {
        cerr << "Error: Could not open JSON file." << endl;
        return;
    }

    // Read the JSON file into a string.
    std::string jsonString;
    jsonFile.seekg(0, std::ios::end);
    jsonString.reserve(jsonFile.tellg());
    jsonFile.seekg(0, std::ios::beg);
    jsonString.assign((std::istreambuf_iterator<char>(jsonFile)), std::istreambuf_iterator<char>());

    // Parse the JSON data.
    json jsonData = json::parse(jsonString);
    verbose = jsonData["verbose"];
    debug = jsonData["debug"];
    isLog = jsonData["isLog_eff"];
    int const nPositions = jsonData["parameters"]["nPositions"];

    vector<double> vPtBinsEdges = jsonData["parameters"]["ptBinsEdges"];
    for (size_t i = 0; i < (nPtBins+1); ++i) {
            ptBinsEdges[i] = vPtBinsEdges[i];
    }

    TString const fileNameRecoPrompt = TString(jsonData["files"]["filenameRecoPrompt"].get<std::string>());
    TString const fileNameRecoNonPrompt = TString(jsonData["files"]["filenameRecoNonPrompt"].get<std::string>());
    TString const fileNameGenPrompt = TString(jsonData["files"]["filenameGenPrompt"].get<std::string>());
    TString const fileNameGenNonPrompt = TString(jsonData["files"]["filenameGenNonPrompt"].get<std::string>());
    std::string const outputDirName = jsonData["outputDir"].get<std::string>();
    TString const hPromptName = TString(jsonData["histos"]["hPromptName"].get<std::string>());
    TString const hNonPromptName = TString(jsonData["histos"]["hNonPromptName"].get<std::string>());

    // Load data
    TFile *fRecoPrompt = TFile::Open(fileNameRecoPrompt);
    TFile *fRecoNonPrompt = TFile::Open(fileNameRecoNonPrompt);
    TFile *fGenPrompt = TFile::Open(fileNameGenPrompt);
    TFile *fGenNonPrompt = TFile::Open(fileNameGenNonPrompt);
    if (!fRecoPrompt || !fRecoNonPrompt) {
        cerr << "ERROR: files reco not found\n";
        return;
    }
    if (!fGenPrompt || !fGenNonPrompt) {
        cerr << "ERROR: files gen not found\n";
        return;
    }

    vector<vector<double>> effsVsPtVsPosPrompt = {};
    vector<vector<double>> effsVsPtVsPosNonPrompt = {};
    vector<vector<double>> effsVsPtVsPosErrPrompt = {};
    vector<vector<double>> effsVsPtVsPosErrNonPrompt = {};
    for (int iPos = 0; iPos < nPositions; iPos++) {
        if (verbose) {
            cout << "----------- POSITION " << iPos << " ----------- \n";
        }
        auto outputName = outputDirName + "Eff_times_Acc_pos" + std::to_string(iPos) + "_bdtSys.root";
        //auto outputName = "Eff_times_Acc_noFDCuts" + std::to_string(iPos) + ".root";
        TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");
        outputFile->cd();
        TString suffixPos = "_pos"+ std::to_string(iPos);
        
        auto dirReco = "pos"+std::to_string(iPos);
        if (verbose) {
            cout << "dirReco: " << dirReco << endl;
        }
        std::vector<TH1F *> vEffs = MakeAnalysis(fRecoPrompt, fRecoNonPrompt, fGenPrompt, fGenNonPrompt, outputFile, dirReco.c_str(), suffixPos, hPromptName, hNonPromptName, isLog);
        TH1F *hEffPrompt = vEffs[0];
        TH1F *hEffNonPrompt = vEffs[1];

        vector<double> effVsPtPrompt = {};
        vector<double> effVsPtNonPrompt = {};
        vector<double> effVsPtErrPrompt = {};
        vector<double> effVsPtErrNonPrompt = {};
        for (int ipt = 0; ipt < nPtBins; ipt++) {
            double iEffPrompt = hEffPrompt->GetBinContent(ipt+1);
            double iEffNonPrompt = hEffNonPrompt->GetBinContent(ipt+1);

            double iEffErrPrompt = hEffPrompt->GetBinError(ipt+1);
            double iEffErrNonPrompt = hEffNonPrompt->GetBinError(ipt+1);

            effVsPtPrompt.push_back(iEffPrompt);
            effVsPtNonPrompt.push_back(iEffNonPrompt);

            effVsPtErrPrompt.push_back(iEffErrPrompt);
            effVsPtErrNonPrompt.push_back(iEffErrNonPrompt);
        }
        effsVsPtVsPosPrompt.push_back(effVsPtPrompt);
        effsVsPtVsPosNonPrompt.push_back(effVsPtNonPrompt);

        effsVsPtVsPosErrPrompt.push_back(effVsPtErrPrompt);
        effsVsPtVsPosErrNonPrompt.push_back(effVsPtErrNonPrompt);
    }
    TCanvas *cEffVsPos = new TCanvas("cEffVsPos", "Efficiency per Cut", 1200, 800);
    cEffVsPos->Divide(ceil(nPtBins/2.0), 2); // Arrange canvases

    std::vector<TGraphErrors *> graphsPrompt;
    std::vector<TGraphErrors *> graphsNonPrompt;

    TFile *fout = TFile::Open("/home/andrea/MyMacros/Double-Dzero-correlations/res/results/thesis_august_2025/systematics/bdt_selections/efficiencies/efficiencies_bdt.root", "RECREATE");

    for (int ipt = 0; ipt < nPtBins; ipt++) {
        std::vector<double> posIndices;
        std::vector<double> efficienciesPrompt;
        std::vector<double> efficienciesNonPrompt;
        std::vector<double> efficienciesErrPrompt;
        std::vector<double> efficienciesErrNonPrompt;

        double sumEffPrompt = 0, sumEffNonPrompt = 0;
        double sumEffErrPrompt = 0, sumEffErrNonPrompt = 0;

        for (int iPos = 0; iPos < nPositions; iPos++) {
            posIndices.push_back(iPos);
            efficienciesPrompt.push_back(effsVsPtVsPosPrompt[iPos][ipt]);
            efficienciesNonPrompt.push_back(effsVsPtVsPosNonPrompt[iPos][ipt]);

            double errPrompt = effsVsPtVsPosErrPrompt[iPos][ipt];
            double errNonPrompt = effsVsPtVsPosErrNonPrompt[iPos][ipt];

            efficienciesErrPrompt.push_back(errPrompt);
            efficienciesErrNonPrompt.push_back(errNonPrompt);

            sumEffPrompt += effsVsPtVsPosPrompt[iPos][ipt];
            sumEffNonPrompt += effsVsPtVsPosNonPrompt[iPos][ipt];

            sumEffErrPrompt += errPrompt * errPrompt;  // Sum of squares for error propagation
            sumEffErrNonPrompt += errNonPrompt * errNonPrompt;
        }

        sumEffErrPrompt = std::sqrt(sumEffErrPrompt);
        sumEffErrNonPrompt = std::sqrt(sumEffErrNonPrompt);

        double meanEffPrompt = sumEffPrompt / nPositions;
        double meanEffNonPrompt = sumEffNonPrompt / nPositions;

        // Compute the 20% upper and lower limits
        double upperBoundPrompt = meanEffPrompt * 1.2;
        double lowerBoundPrompt = meanEffPrompt * 0.8;

        graphsPrompt.push_back(new TGraphErrors(
            nPositions, 
            posIndices.data(), 
            efficienciesPrompt.data(), 
            nullptr,  // No x-axis errors
            efficienciesErrPrompt.data()  // y-axis errors
        ));
        
        graphsNonPrompt.push_back(new TGraphErrors(
            nPositions, 
            posIndices.data(), 
            efficienciesNonPrompt.data(), 
            nullptr,  // No x-axis errors
            efficienciesErrNonPrompt.data()  // y-axis errors
        ));

        cEffVsPos->cd(ipt + 1);
        graphsPrompt[ipt]->SetMarkerStyle(20);
        graphsPrompt[ipt]->SetMarkerColor(kBlack);
        graphsPrompt[ipt]->GetYaxis()->SetRangeUser(lowerBoundPrompt*0.8, upperBoundPrompt*1.2);

        graphsPrompt[ipt]->SetTitle(Form("Efficiency vs Position for pT bin %d", ipt));
        graphsPrompt[ipt]->GetXaxis()->SetTitle("Position index");
        graphsPrompt[ipt]->GetYaxis()->SetTitle("Efficiency");
        graphsPrompt[ipt]->Draw("AP");

        graphsPrompt[ipt]->SetName(Form("gEffPerCut_%d_%d", ipt, ipt+1));
        fout->cd();
        graphsPrompt[ipt]->Write();



        // Draw horizontal lines at mean efficiencies
        double xmin = posIndices.front();
        double xmax = posIndices.back();
    
        TLine *linePrompt = new TLine(xmin, meanEffPrompt, xmax+1, meanEffPrompt);
        linePrompt->SetLineColor(kRed);
        linePrompt->SetLineStyle(2); // Dashed line
        linePrompt->Draw("same");

        // Add 20% upper and lower bounds for prompt
        TLine *lineUpperPrompt = new TLine(xmin, upperBoundPrompt, xmax+1, upperBoundPrompt);
        lineUpperPrompt->SetLineColor(kRed);
        lineUpperPrompt->SetLineStyle(3); // Dotted line
        lineUpperPrompt->Draw("same");

        TLine *lineLowerPrompt = new TLine(xmin, lowerBoundPrompt, xmax+1, lowerBoundPrompt);
        lineLowerPrompt->SetLineColor(kRed);
        lineLowerPrompt->SetLineStyle(3); // Dotted line
        lineLowerPrompt->Draw("same");

        auto legend = new TLegend(0.6, 0.2, 0.8, 0.4);
        legend->AddEntry(graphsPrompt[ipt], "Prompt", "p");
        legend->Draw();
    }

    cEffVsPos->SaveAs("Efficiency_vs_Position_per_pT_bin.pdf");
    fout->Close();

}

// Function to plot and save a single TH1F
void DrawHist(TH1F *h, TString saveName, TFile *f, TString yAxisLabel, bool isLog) {
    TCanvas *c = new TCanvas();
    c->SetTickx();
    c->SetTicky();
    if (isLog) {
        c->SetLogy();
    }

    h->GetYaxis()->SetTitle(yAxisLabel);
    h->GetXaxis()->SetTitle("p_{T}");
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
    h->SetMarkerStyle(21);
    h->Draw("EP");

    //c->SaveAs(saveName);

    f->cd();
    if (isLog) {
        c->Write();
    }
    h->Write();
    if (c) {
        c->Close();
        gSystem->ProcessEvents();
        c->Clear();
    }
}

// Function to plot and save two TH1Fs together
void DrawHist(TH1F *h, TH1F *h2, TString saveName, TFile *f, TString yAxisLabel, bool isLog) {
    TCanvas *c = new TCanvas();
    c->SetTickx();
    c->SetTicky();
    if (isLog) {
        c->SetLogy();
    }

    h->GetYaxis()->SetTitle(yAxisLabel);
    h->GetXaxis()->SetTitle("p_{T}");
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
    h->SetMarkerStyle(21);
    h->Draw("EP");

    h2->GetYaxis()->SetTitle(yAxisLabel);
    h2->GetXaxis()->SetTitle("p_{T}");
    h2->SetMarkerColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->SetMarkerStyle(22);
    h2->Draw("sameEP");

    TLegend* legend = new TLegend(0.5, 0.2, 0.85, 0.5); // Adjust coordinates as needed
    legend->AddEntry(h, "Prompt", "lep");
    legend->AddEntry(h2, "Non-prompt", "lep");
    legend->Draw();

    //c->SaveAs(saveName);


    f->cd();
    if (isLog) {
        c->Write();
    }
    h->Write();
    if (c) {
        c->Close();
        gSystem->ProcessEvents();
        c->Clear();
    }
}

// Function to calculate the raw-yield
TH1F * GetRawYield(std::vector<TH1F *> vh, TString const& name) {
    auto title = "Raw Yield"+name;
    TH1F *hYield = new TH1F(name, title, nPtBins, ptBinsEdges);
    if (verbose) {std::cout << "------- " << name << " ------- \n";}
    for (int i = 0; i < nPtBins; i++) {
        if (vh[i]->Integral() != 0) {
            Double_t yieldErr = 0.;
            hYield->SetBinContent(i+1, vh[i]->IntegralAndError(0, vh[i]->GetXaxis()->GetNbins(), yieldErr));
            hYield->SetBinError(i+1, yieldErr);
            if (verbose) {
                std::cout << hYield->GetBinContent(i+1) << " +- " << hYield->GetBinError(i+1) << std::endl;
            }
        }
    }
    return hYield;
}

// Calculate efficiencies as yield_reco / yield_gen
TH1F * CalculateEfficiency(TH1F *hYieldGen, TH1F *hYieldReco, TString nameReco) {
    auto name = "hEff" + nameReco;
    auto title = "Efficiency " + nameReco;
    if (verbose) {
        std::cout << "------- Eficiencias " << nameReco << " -------- \n";
    }
    TH1F *hEff = new TH1F(name, title, nPtBins, ptBinsEdges);
    for (int i = 0; i < nPtBins; i++) {
        hEff->SetBinContent(i+1, hYieldReco->GetBinContent(i+1)/hYieldGen->GetBinContent(i+1));
        hEff->SetBinError(i+1, hYieldReco->GetBinContent(i+1)/hYieldGen->GetBinContent(i+1) *
                          (hYieldReco->GetBinError(i+1)/hYieldReco->GetBinContent(i+1) + hYieldGen->GetBinError(i+1)/hYieldGen->GetBinContent(i+1)));
        if (verbose) {
            std::cout << hEff->GetBinContent(i+1) << " +- " << hEff->GetBinError(i+1) << std::endl;
        }
    }
    return hEff;
}

// Get vector of THF1s containing Y per selected pT range 
std::vector<TH1F *> GetYPerPtRange(int const& nbins, TH2F* hPtVsY, TString nameParticle) {
    auto xProjName = "hPt"+nameParticle;
    TH1F *hPt = (TH1F *)hPtVsY->ProjectionX(xProjName, 0, nbins);
    std::vector<TH1F *> vHYPt = {};
    if (verbose) {
        std::cout << "pT projections in y" << std::endl;
    }
    for (int i=0; i < nPtBins; i++) {
        TString yProjName = Form("pT%d", i)+nameParticle;

        vHYPt.push_back((TH1F *)hPtVsY->ProjectionY(yProjName, hPt->FindBin(ptBinsEdges[i]), hPt->FindBin(ptBinsEdges[i+1])-1));
        if (debug) {
            std::cout << "                                   " << hPt->FindBin(ptBinsEdges[i]) << " " << hPt->FindBin(ptBinsEdges[i+1])-1 << std::endl;
        }
        if (verbose) {
            std::cout << hPt->GetBinContent(i+1) << " +- " << hPt->GetBinError(i+1) << std::endl;
        }
    }
    return vHYPt;
}

std::vector<TH1F*> MakeAnalysis(TFile *fRecoPrompt, TFile *fRecoNonPrompt, TFile *fGenPrompt, TFile *fGenNonPrompt, TFile *outputFile,
                    TString const& dirnameReco, TString const& suffix, TString const& hPromptName, TString const& hNonPromptName, bool const& isLog) {
    TString dirnameGen = "hf-correlator-d-meson-pairs";
    TDirectory *dirGenPrompt = (TDirectory *)fGenPrompt->Get(dirnameGen);
    TDirectory *dirGenNonPrompt = (TDirectory *)fGenNonPrompt->Get(dirnameGen);
    TDirectory *dirRecoPrompt = (TDirectory *)fRecoPrompt->Get(dirnameReco);
    TDirectory *dirRecoNonPrompt = (TDirectory *)fRecoNonPrompt->Get(dirnameReco);

    // Load MC histos
    TH3F *hPtVsYVsNContribGenPrompt = (TH3F *)dirGenPrompt->Get("hPtVsYVsNContribMcGenPrompt");
    TH3F *hPtVsYVsNContribGenNonPrompt = (TH3F *)dirGenNonPrompt->Get("hPtVsYVsNContribMcGenNonPrompt");

    TH2F *hPtVsYGenPrompt = (TH2F *)hPtVsYVsNContribGenPrompt->Project3D("yx");
    TH2F *hPtVsYGenNonPrompt = (TH2F *)hPtVsYVsNContribGenNonPrompt->Project3D("yx");

    TString recoHistPromptName =  hPromptName;
    TString recoHistNonPromptName =  hNonPromptName;
    TH2F *hPtVsYRecoPrompt = (TH2F *)dirRecoPrompt->Get(recoHistPromptName);
    TH2F *hPtVsYRecoNonPrompt = (TH2F *)dirRecoNonPrompt->Get(recoHistNonPromptName);
    if (verbose) {
        cout << "            Data loaded \n";
    }

    // Project hPtVsY as a function of pT
    int nbins = hPtVsYGenPrompt->GetXaxis()->GetNbins();

    // Project hPtVsY as a function of y
    int const nbinsYGen = hPtVsYGenPrompt->GetYaxis()->GetNbins();
    int const nbinsYRecoPrompt = hPtVsYRecoPrompt->GetYaxis()->GetNbins();
    int const nbinsYRecoNonPrompt = hPtVsYRecoNonPrompt->GetYaxis()->GetNbins();

    if (debug) {
        std::cout << "Number of pT bins " << suffix << ": " << nbins << std::endl;
        std::cout << "Min pT bin: " << hPtVsYGenPrompt->GetXaxis()->GetBinLowEdge(0)
        << "     ; Max pT bin: " << hPtVsYGenPrompt->GetXaxis()->GetBinUpEdge(nbins) << std::endl;
        std::cout << "Number of bins Y gen " << suffix << ": " << nbinsYGen << std::endl;
        std::cout << "Number of bins Y reco prompt" << suffix << ": " << nbinsYRecoPrompt << std::endl;
        std::cout << "Number of bins Y reco non-prompt" << suffix << ": " << nbinsYRecoNonPrompt << std::endl;
    }

    TString genPromptName = "GenPrompt"+suffix, genNonPromptName = "GenNonPrompt"+suffix;
    TString recoPromptName = "RecoPrompt"+suffix, recoNonPromptName = "RecoNonPrompt"+suffix;
    std::vector<TH1F *> hYGenPromptPt = GetYPerPtRange(nbinsYGen, hPtVsYGenPrompt, genPromptName);
    std::vector<TH1F *> hYGenNonPromptPt = GetYPerPtRange(nbinsYGen, hPtVsYGenNonPrompt, genNonPromptName);
    std::vector<TH1F *> hYRecoPtPrompt = GetYPerPtRange(nbinsYRecoPrompt, hPtVsYRecoPrompt, recoPromptName);
    std::vector<TH1F *> hYRecoPtNonPrompt = GetYPerPtRange(nbinsYRecoNonPrompt, hPtVsYRecoNonPrompt, recoNonPromptName);

    // Calculate integrals, which will give us the raw yields
    // Gen
    TH1F *hYieldGenPrompt = GetRawYield(hYGenPromptPt, "yieldGenPrompt"+suffix);
    TH1F *hYieldGenNonPrompt = GetRawYield(hYGenNonPromptPt, "yieldGenNonPrompt"+suffix);
    // Reco
    TH1F *hYieldRecoPrompt = GetRawYield(hYRecoPtPrompt, "yieldRecoPrompt"+suffix);
    TH1F *hYieldRecoNonPrompt = GetRawYield(hYRecoPtNonPrompt, "yieldRecoNonPrompt"+suffix);

    // Calculate the efficiencies
    TH1F *hEffNonPrompt = CalculateEfficiency(hYieldGenNonPrompt, hYieldRecoNonPrompt, "NonPrompt");
    DrawHist(hEffNonPrompt, "efficiencyNonPrompt"+suffix+".pdf", outputFile, "Efficiency", isLog);
    TH1F *hEffPrompt = CalculateEfficiency(hYieldGenPrompt, hYieldRecoPrompt, "Prompt");
    DrawHist(hEffPrompt, "efficiencyPrompt"+suffix+".pdf", outputFile, "Efficiency"+suffix, isLog);
    DrawHist(hEffPrompt, hEffNonPrompt, "efficiencyPrompt"+suffix+".pdf", outputFile, "Efficiency"+suffix, isLog);

    std::vector<TH1F *> vEffs = {hEffPrompt, hEffNonPrompt};
    return vEffs;
}
