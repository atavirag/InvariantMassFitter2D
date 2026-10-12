// Task para ajustar el pico de masa invariante del archivo de datos del Run3

#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TMath.h"
#include "THnSparse.h"
#include "THStack.h"
#include "TF1.h"
#include "TSystem.h"

#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

void DrawProjection(TH1F* hD, TH1F* hDbar, int const& ptMin, int const& ptMax, TString filename);

void PeakAnalysisMc()
{
    // Cargamos la configuración escrita en el JSON
    ifstream jsonFile("config-PeakAnalysis.json");
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
    bool const verbose = jsonData["verbose"];

    // Cargamos los datos
    TString const fileName = jsonData["filename"];
    TString const dirName = jsonData["dirname"];
    TString const histName = jsonData["histname"];
    TString const outputName = jsonData["outputfile"];
    TString const outputDir = jsonData["outputdir"];
    TFile *f = TFile::Open(fileName);
    TDirectory *dir = (TDirectory *)f->Get(dirName);

    if (!jsonData["configMC"]["isMC"]) {
        cerr << "Error: Enable MC." << endl;
        return;
    }
    TString const histNameD0 = jsonData["configMC"]["histnameMCD0"];
    TString const histNameD0bar = jsonData["configMC"]["histnameMCD0bar"];

    // Cargamos los datos

    TH2F *hMass = (TH2F *)dir->Get(histName);
    TH3F *hMassSigD0 = (TH3F *)dir->Get(histNameD0);
    TH3F *hMassSigD0bar = (TH3F *)dir->Get(histNameD0bar);

    if (verbose) {std::cout << "Datos cargados \n";}

    // Proyectamos hMass en función del pT
    int const nbinsY = hMassSigD0->GetZaxis()->GetNbins();
    int const nbinsPt = hMassSigD0->GetYaxis()->GetNbins();
    if (verbose) {
        std::cout << "Bines y: " << nbinsY << " Bines pT: "<< nbinsPt << std::endl;
        std::cout << "   TAMAÑO BINES   \n";
        for (int i=0; i<nbinsPt; i++) {
            std::cout << i << ".-   Low: " << hMassSigD0->GetYaxis()->GetBinLowEdge(i) << "     Up: " << hMassSigD0->GetYaxis()->GetBinUpEdge(i) << std::endl;
        }
    }

    TFile *outputFile = TFile::Open(outputName, "RECREATE");
    outputFile->mkdir(outputDir);
    outputFile->cd(outputDir);

    // Rangos de pT
    int const nPtBins = jsonData["nPtBins"]; 
    vector<double> ptBinsEdges = jsonData["ptBinsEdges"];
    TH1F *hMassPt = (TH1F *)hMassSigD0->ProjectionX("hMassPt", 0, nbinsPt);
    TH1F *hMassPtbar = (TH1F *)hMassSigD0bar->ProjectionX("hMassPtbar", 0, nbinsPt);
    std::vector<TH1F *> vMassPt = {};
    std::vector<TH1F *> vMassPtbar = {};

    for (int i=0; i < nPtBins; i++) {
        auto nameMass = "pT"+std::to_string(i)+"Mass";
        auto nameMassBar = "pT"+std::to_string(i)+"MassBar";

        vMassPt.push_back((TH1F *)hMassSigD0->ProjectionX(nameMass.c_str(), hMassSigD0->GetYaxis()->FindBin(ptBinsEdges[i]), hMassSigD0->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));
        vMassPtbar.push_back((TH1F *)hMassSigD0bar->ProjectionX(nameMassBar.c_str(), hMassSigD0bar->GetYaxis()->FindBin(ptBinsEdges[i]), hMassSigD0bar->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));
        if (verbose) {
            std::cout << "                                   " << hMassSigD0bar->GetYaxis()->FindBin(ptBinsEdges[i]) << " " << hMassSigD0bar->GetYaxis()->FindBin(ptBinsEdges[i+1])-1 << std::endl;
        }
    }

    // Pintamos las proyecciones
    for (int i=0; i < nPtBins; i++) {
        auto projName = "InvMassPt"+std::to_string(i);
        DrawProjection(vMassPt[i], vMassPtbar[i], ptBinsEdges[i], ptBinsEdges[i+1], projName.c_str());
    }
}

void DrawProjection(TH1F* hD, TH1F* hDbar, int const& ptMin, int const& ptMax, TString filename) {
    TCanvas *c = new TCanvas();
    c->SetTickx();
    c->SetTicky();
    c->SetGrid();

    auto title = TString::Format("Invariant Mass, pT =(%d, %d); Mass (GeV); Entries", ptMin, ptMax);
    hD->Add(hDbar);
    hD->SetMarkerStyle(7);
    hD->SetMarkerColor(kBlue);
    hD->SetTitle(title);
    hD->GetXaxis()->SetRangeUser(1.65, 2.10);
    hD->Draw("PE");

    hD->Write("",TObject::kOverwrite);
    auto savePdf = filename+".pdf";
    //c->SaveAs(savePdf);
    if (c)
    {
        c->Close();
        gSystem->ProcessEvents();
        c->Clear();
    }

}