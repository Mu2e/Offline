#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMultigraph.h"
#include "TCanvas.h"
#include <sstream>
#include "dataAnalysis.C"

int plotRejectionRate()
{
	TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

	int numberOfBins = 2000;

	TFile e("electronFit25ns_8DoublePeak6Uniform.root");
	TTree* electronData25ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F electronHist25ns8 = TH1F("electronHist25ns8","electronHist25ns8",numberOfBins,0.0,0.2); 
	electronData25ns8->Draw("fittingFunction->GetParameter(1)/7.0e6>>electronHist25ns8");

	TH1F electronHistMCTrig = TH1F("electronHistMCTrig","electronHistMCTrig",numberOfBins,0.0,0.2);
	electronData25ns8->Draw("qMctrigenergy>>electronHistMCTrig");

	TH1F electronHistSum = TH1F("electronHistSum","electronHistSum",numberOfBins,0.0,10000.0);
	electronData25ns8->Draw("(qAdc[0]+qAdc[1]+qAdc[2]+qAdc[3]+qAdc[4]+qAdc[5]+qAdc[6]+qAdc[7] - 4*(qAdc[0] + qAdc[1]))>>electronHistSum");

	TFile f("electronFit25ns_8.root");
	TTree* electronData25ns81 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F electronHist25ns81 = TH1F("electronHist25ns81","electronHist25ns81",numberOfBins,0.0,0.2); 
	electronData25ns81->Draw("func->GetParameter(1)/7.0e6>>electronHist25ns81");

	TFile g("electronFit25ns_8DoublePeak6Gaussian.root");
	TTree* electronData25ns8Gauss = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F electronHist25ns8Gauss = TH1F("electronHist25ns8Gauss","electronHist25ns8Gauss",numberOfBins,0.0,0.2); 
	electronData25ns8Gauss->Draw("fittingFunction->GetParameter(1)/7.0e6>>electronHist25ns8Gauss");

	TFile k("protonFit25ns_8DoublePeak6Uniform.root");
	TTree* protonData25ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F protonHist25ns8 = TH1F("protonHist25ns8","protonHist25ns8",numberOfBins,0.0,0.2); 
	protonData25ns8->Draw("fittingFunction->GetParameter(1)/7.0e6>>protonHist25ns8");

	TH1F protonHistMCTrig = TH1F("protonHistMCTrig","protonHistMCTrig",numberOfBins,0.0,0.2);
	protonData25ns8->Draw("qMctrigenergy>>protonHistMCTrig");

	TFile l("protonFit25ns_8DoublePeak6Gaussian.root");
	TTree* protonData25ns8Gauss = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F protonHist25ns8Gauss = TH1F("protonHist25ns8Gauss","protonHist25ns8Gauss",numberOfBins,0.0,0.2); 
	protonData25ns8Gauss->Draw("fittingFunction->GetParameter(1)/7.0e6>>protonHist25ns8Gauss");


	TFile j("protonFit25ns_8.root");
	TTree* protonData25ns81 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F protonHist25ns81 = TH1F("protonHist25ns81","protonHist25ns81",numberOfBins,0.0,0.2); 
	protonData25ns81->Draw("func->GetParameter(1)/7.0e6>>protonHist25ns81");


	TH1F protonHistSum = TH1F("protonHistSum","protonHistSum",numberOfBins,0.0,10000.0);
	protonData25ns8->Draw("(qAdc[0]+qAdc[1]+qAdc[2]+qAdc[3]+qAdc[4]+qAdc[5]+qAdc[6]+qAdc[7]- 4.0*(qAdc[0] + qAdc[1]))>>protonHistSum");

	TGraph *gr25ns8 = computeRejectionGraph(electronHist25ns8,protonHist25ns8,numberOfBins);
	gr25ns8->SetMarkerColor(kRed);
	gr25ns8->SetLineColor(kRed);

	TGraph *grSum = computeRejectionGraph(electronHistSum,protonHistSum,numberOfBins);
	grSum->SetMarkerColor(kBlue);
	grSum->SetLineColor(kBlue);

	TGraph *gr25ns81 = computeRejectionGraph(electronHist25ns81,protonHist25ns81,numberOfBins);
	gr25ns81->SetMarkerColor(kGreen);
	gr25ns81->SetLineColor(kGreen);

	TGraph *gr25nsMCTrig = computeRejectionGraph(electronHistMCTrig,protonHistMCTrig,numberOfBins);
	gr25nsMCTrig->SetMarkerColor(kOrange);
	gr25nsMCTrig->SetLineColor(kOrange);

	TGraph *gr25nsGauss = computeRejectionGraph(electronHist25ns8Gauss,protonHist25ns8Gauss,numberOfBins);
	gr25nsGauss->SetMarkerColor(kBlack);
	gr25nsGauss->SetLineColor(kBlack);





	 TMultiGraph *mg = new TMultiGraph();
	 mg->SetTitle("title;xaxis title; yaxis title");
	 mg->Add(gr25ns8);
	 mg->Add(gr25ns81);
	 mg->Add(grSum);
	 mg->Add(gr25nsMCTrig);
	 mg->Add(gr25nsGauss);


   	 mg->Draw("apl");

   	mg->SetTitle("Rejection Rate of Protons vs. Acceptance Rate of Electrons");
   //	mg->CenterTitle();
   	mg->GetXaxis()->SetTitle("Acceptance Rate of Electrons");
   	mg->GetYaxis()->SetTitle("Rejection Rate of Protons");

   	TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
   	legend->SetTextFont(72);
   	legend->SetTextSize(0.04);
   	legend->SetFillColor(0);
   	legend->AddEntry(gr25ns8, "Uniform", "lp");
   	legend->AddEntry(grSum, "Sum", "lp");
   	legend->AddEntry(gr25ns81, "Gauss", "lp");
   	legend->AddEntry(gr25nsMCTrig, "MC Trigenergy", "lp");
   	legend->AddEntry(gr25nsGauss, "Gauss - New", "lp");
   	legend->Draw();


   	return 0;
}

TGraph* computeRejectionGraph(TH1F &electronHist, TH1F &protonHist, const int numberOfBins)
{
  Double_t truncX[numberOfBins], truncY[numberOfBins];

  int protonSum = 0;
  int electronSum = 0;
  for (int i = 1; i <= numberOfBins; ++i)
  {
    // For some reason bin number starts with 1 for TH1F
    protonSum += protonHist.GetBinContent(i);
    electronSum += electronHist.GetBinContent(i);
    // acceptance rate of electrons
    truncX[i - 1] = electronSum / (double) electronHist.GetEntries();
    // 1 - rejection rate
    truncY[i - 1] = 1 - (protonSum / (double) protonHist.GetEntries());
  }

  TGraph * rejectionGraph = new TGraph(numberOfBins,truncX,truncY);
  return rejectionGraph;

}