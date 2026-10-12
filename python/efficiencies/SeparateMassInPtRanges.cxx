// Task para separar el plot de masa invariante procedente de taskD0.cxx en rangos de pT

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

void DrawProjection(TH1F* h, int const& ptMin, int const& ptMax, TString filename);

void SeparateMassInPtRanges()
{
    // Cargamos la configuración escrita en el JSON
    ifstream jsonFile("config-SeparateMassInPtRanges.json");
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
    TH2F *hMass = (TH2F *)dir->Get(histName);
    if (verbose) {std::cout << "Datos cargados \n";}

    // Proyectamos hMass en función del pT
    int nbinsPt = hMass->GetYaxis()->GetNbins();
    if (verbose) {
        std::cout << " Bines pT: " << nbinsPt << std::endl;
        std::cout << "   TAMAÑO BINES   \n";
        for (int i=0; i<nbinsPt; i++) {
            std::cout << i << ".-   Low: " << hMass->GetYaxis()->GetBinLowEdge(i) << "     Up: " << hMass->GetYaxis()->GetBinUpEdge(i) << std::endl;
        }
    }

    TFile *outputFile = TFile::Open(outputDir+outputName, "RECREATE");
    //outputFile->mkdir(outputDir);
    //outputFile->cd(outputDir);

    // Rangos de pT
    int const nPtBins = jsonData["nPtBins"]; 
    vector<double> ptBinsEdges = jsonData["ptBinsEdges"];
    TH1F *hMassPt = (TH1F *)hMass->ProjectionX("hMassPt", 0, nbinsPt);
    std::vector<TH1F *> vMassPt = {};

    for (int i=0; i < nPtBins; i++) {
        auto name = "pT"+std::to_string(i);

        vMassPt.push_back((TH1F *)hMass->ProjectionX(name.c_str(), hMass->GetYaxis()->FindBin(ptBinsEdges[i]), hMass->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));
        if (verbose) {
            std::cout << "                                   " << hMass->GetYaxis()->FindBin(ptBinsEdges[i]) << " " << hMass->GetYaxis()->FindBin(ptBinsEdges[i+1])-1 << std::endl;
        }
    }
    // Pintamos las proyecciones
    for (int i=0; i < nPtBins; i++) {
        auto projName = "InvMassPt"+std::to_string(i);
        DrawProjection(vMassPt[i], ptBinsEdges[i], ptBinsEdges[i+1], projName.c_str());
    }
}

void DrawProjection(TH1F* h, int const& ptMin, int const& ptMax, TString filename) {
    TCanvas *c = new TCanvas();
    c->SetTickx();
    c->SetTicky();
    c->SetGrid();

    h->SetMarkerStyle(7);
    h->SetMarkerColor(kBlue);
    auto title = TString::Format("Invariant Mass, pT =(%d, %d) GeV; Mass (GeV); Entries", ptMin, ptMax);
    h->SetTitle(title);
    h->GetXaxis()->SetRangeUser(1.65, 2.10);
    h->Draw("PE");
    h->Write("",TObject::kOverwrite);

    auto savePdf = filename+".pdf";
    //c->SaveAs(savePdf);
    if (c)
    {
        c->Close();
        gSystem->ProcessEvents();
        c->Clear();
    }
}
