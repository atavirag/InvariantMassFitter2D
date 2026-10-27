#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "Riostream.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

bool verbose, debug;
void computeEfficiencyMap()
{
    // Load JSON config
    ifstream jsonFile("efficiency_map_config.json");
    if (!jsonFile.is_open())
    {
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

// Use .get<std::string>() before assigning to TString
TString const fileNameGen  = TString(jsonData["files"]["filenameGen"].get<std::string>());
TString const fileNameReco = TString(jsonData["files"]["filenameReco"].get<std::string>());
TString const fileNameData = TString(jsonData["files"]["filenameData"].get<std::string>());

std::string const outputDirName = jsonData["outputDir"].get<std::string>();

TString const hNameReco = TString(jsonData["histos"]["hNameReco"].get<std::string>());
TString const hNameGen  = TString(jsonData["histos"]["hNameGen"].get<std::string>());
TString const hNameData = TString(jsonData["histos"]["hNameData"].get<std::string>());

TString const hNameRecoPrompt = TString(jsonData["histos"]["hNameRecoPrompt"].get<std::string>());
TString const hNameGenPrompt  = TString(jsonData["histos"]["hNameGenPrompt"].get<std::string>());

TString const hNameRecoNonPrompt = TString(jsonData["histos"]["hNameRecoNonPrompt"].get<std::string>());
TString const hNameGenNonPrompt  = TString(jsonData["histos"]["hNameGenNonPrompt"].get<std::string>());

    // Load data
    TFile *fGen = TFile::Open(fileNameGen);
    TFile *fReco = TFile::Open(fileNameReco);
    TFile *fData = TFile::Open(fileNameData);

    if (!fReco || !fGen || !fData)
    {
        cerr << "ERROR: files not found\n";
        return;
    }

    std::string const outputFileName = jsonData["outputFile"].get<std::string>();
    auto outputName = outputDirName + outputFileName;
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");
    outputFile->cd();

    //TString dirnameGen = "hf-task-d0";
    TString dirname = "hf-correlator-d-meson-pairs";

    TDirectory *dirReco = (TDirectory *)fReco->Get(dirname);
    TDirectory *dirGen = (TDirectory *)fGen->Get(dirname);
    TDirectory *dirData = (TDirectory *)fData->Get(dirname);

    // Load MC histos
    TH3F *hPtVsYReco = (TH3F *)dirReco->Get(hNameReco);
    TH3F *hPtVsYGen = (TH3F *)dirGen->Get(hNameGen);

    TH3F *hPtVsYRecoPrompt = (TH3F *)dirReco->Get(hNameRecoPrompt);
    TH3F *hPtVsYGenPrompt = (TH3F *)dirGen->Get(hNameGenPrompt);

    TH3F *hPtVsYRecoNonPrompt = (TH3F *)dirReco->Get(hNameRecoNonPrompt);
    TH3F *hPtVsYGenNonPrompt = (TH3F *)dirGen->Get(hNameGenNonPrompt);
    // Load data histo
    TH1F *hNContrib = (TH1F *)dirData->Get(hNameData);


    double dataIntegral = hNContrib->Integral();
    double recoIntegral = hPtVsYReco->ProjectionZ()->Integral();
    double genIntegral = hPtVsYGen->ProjectionZ()->Integral();

    TH1F *hNContribReco = (TH1F *)hPtVsYReco->ProjectionZ();
    TH1F *hNContribGen = (TH1F *)hPtVsYGen->ProjectionZ();

    if (verbose)
    {
        cout << "            Data loaded \n";
    }

    TH2F *hEfficiencyMap = (TH2F *)hPtVsYReco->Project3D("yx");
    hEfficiencyMap->Reset();
    hEfficiencyMap->SetNameTitle("hEfficiencyMap", "Efficiency map");

    TH2F *hPtVsYRecoWeighted = (TH2F *)hPtVsYReco->Project3D("yx");
    hPtVsYRecoWeighted->Reset();
    hPtVsYRecoWeighted->SetNameTitle("hPtVsYRecoWeighted", "hPtVsYRecoWeighted");

    TH2F *hPtVsYGenWeighted = (TH2F *)hPtVsYGen->Project3D("yx");
    hPtVsYGenWeighted->Reset();
    hPtVsYGenWeighted->SetNameTitle("hPtVsYGenWeighted", "hPtVsYGenWeighted");

    TH2F *hEfficiencyMapPrompt = (TH2F *)hPtVsYRecoPrompt->Project3D("yx");
    hEfficiencyMapPrompt->Reset();
    hEfficiencyMapPrompt->SetNameTitle("hEfficiencyMapPrompt", "Efficiency map prompt");

    TH2F *hPtVsYRecoWeightedPrompt = (TH2F *)hPtVsYRecoPrompt->Project3D("yx");
    hPtVsYRecoWeightedPrompt->Reset();
    hPtVsYRecoWeightedPrompt->SetNameTitle("hPtVsYRecoWeightedPrompt", "hPtVsYRecoWeightedPrompt");

    TH2F *hPtVsYGenWeightedPrompt = (TH2F *)hPtVsYGenPrompt->Project3D("yx");
    hPtVsYGenWeightedPrompt->Reset();
    hPtVsYGenWeightedPrompt->SetNameTitle("hPtVsYGenWeightedPrompt", "hPtVsYGenWeightedPrompt");

    TH2F *hEfficiencyMapNonPrompt = (TH2F *)hPtVsYRecoNonPrompt->Project3D("yx");
    hEfficiencyMapNonPrompt->Reset();
    hEfficiencyMapNonPrompt->SetNameTitle("hEfficiencyMapNonPrompt", "Efficiency map NonPrompt");

    TH2F *hPtVsYRecoWeightedNonPrompt = (TH2F *)hPtVsYRecoNonPrompt->Project3D("yx");
    hPtVsYRecoWeightedNonPrompt->Reset();
    hPtVsYRecoWeightedNonPrompt->SetNameTitle("hPtVsYRecoWeightedNonPrompt", "hPtVsYRecoWeightedNonPrompt");

    TH2F *hPtVsYGenWeightedNonPrompt = (TH2F *)hPtVsYGenNonPrompt->Project3D("yx");
    hPtVsYGenWeightedNonPrompt->Reset();
    hPtVsYGenWeightedNonPrompt->SetNameTitle("hPtVsYGenWeightedNonPrompt", "hPtVsYGenWeightedNonPrompt");

    // Get number of bins
    int const nbinsPt = hPtVsYGen->GetXaxis()->GetNbins();
    int const nbinsY = hPtVsYGen->GetYaxis()->GetNbins();
    int const nbinsNContrib = hPtVsYGen->GetZaxis()->GetNbins();

    int y_min = hPtVsYGen->GetYaxis()->FindBin(-0.5);
    int y_max = hPtVsYGen->GetYaxis()->FindBin(0.5);

    if (debug)
    {
        std::cout << "Number of pT bins: " << nbinsPt << std::endl;
        std::cout << "Min pT bin: " << hPtVsYGen->GetXaxis()->GetBinLowEdge(0)
                  << "     ; Max pT bin: " << hPtVsYGen->GetXaxis()->GetBinUpEdge(nbinsPt) << std::endl;
        std::cout << "Number of bins Y: " << nbinsY << std::endl;
        std::cout << "Number of PV contrib bins: " << nbinsNContrib << std::endl;
    }

    // Divide histos point by point
    for (int ipt = 1; ipt <= nbinsPt; ipt++)
    {
        for (int iy = y_min; iy <= y_max; iy++)
        {
            // Weight the distributions using nContrib info
            double weightSumReco = 0.0, weightSumGen = 0.0;

            double weightSumRecoPrompt = 0.0, weightSumGenPrompt = 0.0;
            double weightSumRecoNonPrompt = 0.0, weightSumGenNonPrompt = 0.0;

            for (int in = 1; in <= nbinsNContrib; in++)
            {
                double binContentReco = hPtVsYReco->GetBinContent(ipt, iy, in);
                double binContentGen = hPtVsYGen->GetBinContent(ipt, iy, in);

                double binContentRecoPrompt = hPtVsYRecoPrompt->GetBinContent(ipt, iy, in);
                double binContentGenPrompt = hPtVsYGenPrompt->GetBinContent(ipt, iy, in);

                double binContentRecoNonPrompt = hPtVsYRecoNonPrompt->GetBinContent(ipt, iy, in);
                double binContentGenNonPrompt = hPtVsYGenNonPrompt->GetBinContent(ipt, iy, in);

                double nContribValueData = hNContrib->GetBinContent(in);
                double nContribValueMCReco = hNContribReco->GetBinContent(in); // MC shape in NContrib
                double nContribValueMCGen = hNContribGen->GetBinContent(in); // MC shape in NContrib

                double dataRecoWeight = 1.0;
                double dataGenWeight = 1.0;
                if (nContribValueMCReco > 0 && nContribValueMCGen > 0 && 
                    nContribValueData > 0 && recoIntegral > 0 && genIntegral > 0) {
                    dataRecoWeight = (nContribValueData/dataIntegral) / (nContribValueMCReco / recoIntegral);
                    dataGenWeight = (nContribValueData/dataIntegral) / (nContribValueMCGen / genIntegral);
                }

                weightSumReco += (dataRecoWeight * binContentReco);
                weightSumGen += ( dataGenWeight * binContentGen);

                weightSumRecoPrompt += (dataRecoWeight * binContentRecoPrompt);
                weightSumGenPrompt += ( dataGenWeight * binContentGenPrompt);

                weightSumRecoNonPrompt += (dataRecoWeight * binContentRecoNonPrompt);
                weightSumGenNonPrompt += ( dataGenWeight * binContentGenNonPrompt);

            }
            double tmpEff = (weightSumGen > 0) ? weightSumReco / weightSumGen : 0;
            cout << "tmpEff: " << tmpEff << endl;

            double tmpEffPrompt = (weightSumGen > 0) ? weightSumRecoPrompt / weightSumGenPrompt : 0;
            cout << "tmpEffPrompt: " << tmpEffPrompt << endl;

            double tmpEffNonPrompt = (weightSumGen > 0) ? weightSumRecoNonPrompt / weightSumGenNonPrompt : 0;
            cout << "tmpEffNonPrompt: " << tmpEffNonPrompt << endl;

            // Compute the weighted average (or leave as summed weights)
            if (weightSumReco > 0)
            {
                hPtVsYRecoWeighted->SetBinContent(ipt, iy, weightSumReco);
                //cout << "weightSumReco " << weightSumReco << endl;
            }
            if (weightSumGen > 0)
            {
                hPtVsYGenWeighted->SetBinContent(ipt, iy, weightSumGen);
            }

            if (weightSumRecoPrompt > 0)
            {
                hPtVsYRecoWeightedPrompt->SetBinContent(ipt, iy, weightSumRecoPrompt);
                //cout << "weightSumReco " << weightSumReco << endl;
            }
            if (weightSumGenPrompt > 0)
            {
                hPtVsYGenWeightedPrompt->SetBinContent(ipt, iy, weightSumGenPrompt);
            }

            if (weightSumRecoNonPrompt > 0)
            {
                hPtVsYRecoWeightedNonPrompt->SetBinContent(ipt, iy, weightSumRecoNonPrompt);
                //cout << "weightSumReco " << weightSumReco << endl;
            }
            if (weightSumGenNonPrompt > 0)
            {
                hPtVsYGenWeightedNonPrompt->SetBinContent(ipt, iy, weightSumGenNonPrompt);
            }
        }
    }

    hPtVsYRecoWeighted->RebinX(5);
    hPtVsYGenWeighted->RebinX(5);
    hEfficiencyMap->RebinX(5);

    hPtVsYRecoWeightedPrompt->RebinX(5);
    hPtVsYGenWeightedPrompt->RebinX(5);
    hEfficiencyMapPrompt->RebinX(5);

    hPtVsYRecoWeightedNonPrompt->RebinX(5);
    hPtVsYGenWeightedNonPrompt->RebinX(5);
    hEfficiencyMapNonPrompt->RebinX(5);


    // Divide histos point by point
    for (int ipt = 1; ipt <= hPtVsYRecoWeighted->GetXaxis()->GetNbins(); ipt++)
    {
        for (int iy = y_min; iy < y_max; iy++) {
            auto eff = 0.;
            auto effError = 0.;

            auto effPrompt = 0.;
            auto effErrorPrompt = 0.;

            auto effNonPrompt = 0.;
            auto effErrorNonPrompt = 0.;

            if (hPtVsYRecoWeighted->GetBinContent(ipt, iy) != 0 && hPtVsYGenWeighted->GetBinContent(ipt, iy) != 0)
            {
                eff = hPtVsYRecoWeighted->GetBinContent(ipt, iy) / hPtVsYGenWeighted->GetBinContent(ipt, iy);
                effError = eff * sqrt(pow(hPtVsYRecoWeighted->GetBinError(ipt, iy) / hPtVsYRecoWeighted->GetBinContent(ipt, iy), 2) + pow(hPtVsYGenWeighted->GetBinError(ipt, iy) / hPtVsYGenWeighted->GetBinContent(ipt, iy), 2));
                hEfficiencyMap->SetBinContent(ipt, iy, eff);
                hEfficiencyMap->SetBinError(ipt, iy, effError);
            }

            if (hPtVsYRecoWeightedPrompt->GetBinContent(ipt, iy) != 0 && hPtVsYGenWeightedPrompt->GetBinContent(ipt, iy) != 0)
            {
                effPrompt = hPtVsYRecoWeightedPrompt->GetBinContent(ipt, iy) / hPtVsYGenWeightedPrompt->GetBinContent(ipt, iy);
                effErrorPrompt = effPrompt * sqrt(pow(hPtVsYRecoWeightedPrompt->GetBinError(ipt, iy) / hPtVsYRecoWeightedPrompt->GetBinContent(ipt, iy), 2) + pow(hPtVsYGenWeightedPrompt->GetBinError(ipt, iy) / hPtVsYGenWeightedPrompt->GetBinContent(ipt, iy), 2));
                hEfficiencyMapPrompt->SetBinContent(ipt, iy, effPrompt);
                hEfficiencyMapPrompt->SetBinError(ipt, iy, effErrorPrompt);
            }

            if (hPtVsYRecoWeightedNonPrompt->GetBinContent(ipt, iy) != 0 && hPtVsYGenWeightedNonPrompt->GetBinContent(ipt, iy) != 0)
            {
                effNonPrompt = hPtVsYRecoWeightedNonPrompt->GetBinContent(ipt, iy) / hPtVsYGenWeightedNonPrompt->GetBinContent(ipt, iy);
                effErrorNonPrompt = effNonPrompt * sqrt(pow(hPtVsYRecoWeightedNonPrompt->GetBinError(ipt, iy) / hPtVsYRecoWeightedNonPrompt->GetBinContent(ipt, iy), 2)
                                                         + pow(hPtVsYGenWeightedNonPrompt->GetBinError(ipt, iy) / hPtVsYGenWeightedNonPrompt->GetBinContent(ipt, iy), 2));
                hEfficiencyMapNonPrompt->SetBinContent(ipt, iy, effNonPrompt);
                hEfficiencyMapNonPrompt->SetBinError(ipt, iy, effErrorNonPrompt);
            }
        }
    }

    hEfficiencyMap->GetZaxis()->SetTitle("(Acc x #epsilon)");   // Set Z-axis title
    hEfficiencyMap->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");   // Set Z-axis title
    hEfficiencyMap->GetZaxis()->SetTitleSize(0.05);       // Adjust title size
    hEfficiencyMap->GetZaxis()->SetTitleOffset(1.0);      // Adjust title position
    hEfficiencyMap->SetStats(0);

    hEfficiencyMapPrompt->GetZaxis()->SetTitle("(Acc x #epsilon)");   // Set Z-axis title
    hEfficiencyMapPrompt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");   // Set Z-axis title
    hEfficiencyMapPrompt->GetZaxis()->SetTitleSize(0.05);       // Adjust title size
    hEfficiencyMapPrompt->GetZaxis()->SetTitleOffset(1.0);      // Adjust title position
    hEfficiencyMapPrompt->SetStats(0);

    hEfficiencyMapNonPrompt->GetZaxis()->SetTitle("(Acc x #epsilon)");   // Set Z-axis title
    hEfficiencyMapNonPrompt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");   // Set Z-axis title
    hEfficiencyMapNonPrompt->GetZaxis()->SetTitleSize(0.05);       // Adjust title size
    hEfficiencyMapNonPrompt->GetZaxis()->SetTitleOffset(1.0);      // Adjust title position
    hEfficiencyMapNonPrompt->SetStats(0);

    TPaveText *ptTitle = new TPaveText(0.12, 0.8, 0.4, 0.94, "brNDC");
    ptTitle->SetTextFont(42);
    ptTitle->SetTextSize(0.04);
    ptTitle->SetBorderSize(0);
    ptTitle->SetFillStyle(0);
    ptTitle->SetTextAlign(11);

    // Add text lines
    ptTitle->AddText("#scale[1.35]{This Thesis}");

    TPaveText *pt = new TPaveText(0.12, 0.8, 0.4, 0.94, "brNDC");
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);
    pt->AddText("#scale[1.2]{pp, #sqrt{#it{s}} = 13.6 TeV}");
    pt->AddText("D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.");

    TCanvas *c2d = new TCanvas("cEfficiencyMap", "cEfficiencyMap", 1000, 700);
    hEfficiencyMap->Draw("colz");
    ptTitle->Draw("same");
    pt->Draw("same");

    TCanvas *cPV = new TCanvas("cPVContributors", "cPVContributors", 800, 600);
    // Set margins
cPV->SetTopMargin(0.04);
cPV->SetLeftMargin(0.15);
cPV->SetRightMargin(0.03);
cPV->SetBottomMargin(0.12);

// Remove stats box globally
gStyle->SetOptStat(0);

// Add ticks on all sides
cPV->SetTickx();
cPV->SetTicky();

// Enable grid lines
cPV->SetGridx();
cPV->SetGridy();

// Grid style customization
gStyle->SetGridColor(kGray + 2);   // Light gray
gStyle->SetGridStyle(3);           // Dotted lines
gStyle->SetGridWidth(1);           // Thin lines

// Apply style settings immediately (optional but recommended)
gROOT->ForceStyle();

    hNContrib->Scale(1.0 / dataIntegral);
    hNContribReco->Scale(1.0 / recoIntegral);

    hNContrib->GetYaxis()->SetRangeUser(0.0, 0.02);
    hNContrib->GetXaxis()->SetTitleSize(0.045);
    hNContrib->GetYaxis()->SetTitleSize(0.045);
    hNContrib->GetXaxis()->SetLabelSize(0.04);
    hNContrib->GetYaxis()->SetLabelSize(0.04);
    hNContrib->SetTitle("");
    hNContrib->GetXaxis()->SetTitle("Number of PV contributors");
    hNContrib->GetYaxis()->SetTitle("Normalised entries");

    hNContrib->SetStats(0);
    hNContrib->SetLineWidth(5);
    hNContrib->SetLineColor(kBlue);

    hNContribReco->SetLineWidth(5);
    hNContribReco->SetLineStyle(kDashed);
    hNContribReco->SetLineColor(kRed);

    hNContrib->Draw();
    hNContribReco->Draw("same");
    ptTitle->Draw("same");
    pt->Draw("same");

    auto legend = new TLegend(0.6, 0.7, 0.88, 0.85);
    legend->SetFillStyle(0);  // 0 = transparent
    legend->SetFillColor(0);  // Just to be sure
    legend->SetLineColor(0);  // Just to be sure
    legend->AddEntry(hNContrib, "Data", "l");
    legend->AddEntry(hNContribReco, "Reco", "l");
    legend->Draw();

    outputFile->cd();
    cPV->Write();

    hEfficiencyMap->Write();
    hEfficiencyMapPrompt->Write();
    hEfficiencyMapNonPrompt->Write();

    c2d->Write();

    hPtVsYRecoWeighted->Write();
    hPtVsYGenWeighted->Write();

    hPtVsYRecoWeightedPrompt->Write();
    hPtVsYGenWeightedPrompt->Write();

    hPtVsYRecoWeightedNonPrompt->Write();
    hPtVsYGenWeightedNonPrompt->Write();
    outputFile->Close();

    cout << "Plots saved. Exit!" << endl;
}