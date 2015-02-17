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

int plotRejectionRate()
{
	TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

	int numberOfBins = 2000;

	TFile e("electronFit25ns_8.root");
	TTree* electronData25ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist25ns8 = new TH1F("electronHist25ns8","electronHist25ns8",numberOfBins,0.0,0.2); 
	electronData25ns8->Draw("func->GetParameter(1)/7.0e6>>electronHist25ns8");

	TFile f("electronFit25ns_12.root");
	TTree* electronData25ns12 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist25ns12 = new TH1F("electronHist25ns12","electronHist25ns12",numberOfBins,0.0,0.2); 
	electronData25ns12->Draw("func->GetParameter(1)/7.0e6>>electronHist25ns12");

	TFile g("electronFit25ns_16.root");
	TTree* electronData25ns16 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist25ns16 = new TH1F("electronHist25ns16","electronHist25ns16",numberOfBins,0.0,0.2); 
	electronData25ns16->Draw("func->GetParameter(1)/7.0e6>>electronHist25ns16");

	TFile h("electronFit50ns_8.root");
	TTree* electronData50ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist50ns8 = new TH1F("electronHist50ns8","electronHist50ns8",numberOfBins,0.0,0.2); 
	electronData50ns8->Draw("func->GetParameter(1)/7.0e6>>electronHist50ns8");

	TFile i("electronFit50ns_12.root");
	TTree* electronData50ns12 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist50ns12 = new TH1F("electronHist50ns12","electronHist50ns12",numberOfBins,0.0,0.2); 
	electronData50ns12->Draw("func->GetParameter(1)/7.0e6>>electronHist50ns12");

	TFile j("electronFit50ns_16.root");
	TTree* electronData50ns16 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *electronHist50ns16 = new TH1F("electronHist50ns16","electronHist50ns16",numberOfBins,0.0,0.2); 
	electronData50ns16->Draw("func->GetParameter(1)/7.0e6>>electronHist50ns16");

	TFile k("protonFit25ns_8.root");
	TTree* protonData25ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist25ns8 = new TH1F("protonHist25ns8","protonHist25ns8",numberOfBins,0.0,0.2); 
	protonData25ns8->Draw("func->GetParameter(1)/7.0e6>>protonHist25ns8");

	TFile l("protonFit25ns_12.root");
	TTree* protonData25ns12 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist25ns12 = new TH1F("protonHist25ns12","protonHist25ns12",numberOfBins,0.0,0.2); 
	protonData25ns12->Draw("func->GetParameter(1)/7.0e6>>protonHist25ns12");

	TFile m("protonFit25ns_16.root");
	TTree* protonData25ns16 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist25ns16 = new TH1F("protonHist25ns16","protonHist25ns16",numberOfBins,0.0,0.2); 
	protonData25ns16->Draw("func->GetParameter(1)/7.0e6>>protonHist25ns16");

	TFile n("protonFit50ns_8.root");
	TTree* protonData50ns8 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist50ns8 = new TH1F("protonHist50ns8","protonHist50ns8",numberOfBins,0.0,0.2); 
	protonData50ns8->Draw("func->GetParameter(1)/7.0e6>>protonHist50ns8");

	TFile o("protonFit50ns_12.root");
	TTree* protonData50ns12 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist50ns12 = new TH1F("protonHist50ns12","protonHist50ns12",numberOfBins,0.0,0.2); 
	protonData50ns12->Draw("func->GetParameter(1)/7.0e6>>protonHist50ns12");

	TFile p("protonFit50ns_16.root");
	TTree* protonData50ns16 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	TH1F *protonHist50ns16 = new TH1F("protonHist50ns16","protonHist50ns16",numberOfBins,0.0,0.2); 
	protonData50ns16->Draw("func->GetParameter(1)/7.0e6>>protonHist50ns16");

	TGraph *gr25ns8 = computeRejectionGraph(electronHist25ns8,protonHist25ns8,numberOfBins);
	gr25ns8->SetMarkerColor(kRed);
	gr25ns8->SetLineColor(kRed);

	TGraph *gr25ns12 = computeRejectionGraph(electronHist25ns12,protonHist25ns12,numberOfBins);
	gr25ns12->SetMarkerColor(kGreen);
	gr25ns12->SetLineColor(kGreen);

	TGraph *gr25ns16 = computeRejectionGraph(electronHist25ns16,protonHist25ns16,numberOfBins);
	gr25ns16->SetMarkerColor(kViolet);
	gr25ns16->SetLineColor(kViolet);

	TGraph *gr50ns8 = computeRejectionGraph(electronHist50ns8,protonHist50ns8,numberOfBins);
	gr50ns8->SetMarkerColor(kOrange);
	gr50ns8->SetLineColor(kOrange);

	TGraph *gr50ns12 = computeRejectionGraph(electronHist50ns12,protonHist50ns12,numberOfBins);
	gr50ns12->SetMarkerColor(kBlack);
	gr50ns12->SetLineColor(kBlack);

	TGraph *gr50ns16 = computeRejectionGraph(electronHist50ns16,protonHist50ns16,numberOfBins);
	gr50ns16->SetMarkerColor(kBlue);
	gr50ns16->SetLineColor(kBlue);



	 TMultiGraph *mg = new TMultiGraph();
	 mg->SetTitle("title;xaxis title; yaxis title");
	 mg->Add(gr25ns8);
	 mg->Add(gr25ns12);
	 mg->Add(gr25ns16);
	 mg->Add(gr50ns8);
	 mg->Add(gr50ns12);
	 mg->Add(gr50ns16);


   	 mg->Draw("apl");

   	mg->SetTitle("(1 - Rejection) Rate of Protons vs. Acceptance Rate of Electrons");
   //	mg->CenterTitle();
   	mg->GetXaxis()->SetTitle("Acceptance Rate of Electrons");
   	mg->GetYaxis()->SetTitle("(1 - Rejection) Rate of Protons");

   	TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
   	legend->SetTextFont(72);
   	legend->SetTextSize(0.04);
   	legend->SetFillColor(0);
   	legend->AddEntry(gr25ns8, "25 ns - 8", "lp");
   	legend->AddEntry(gr25ns12, "25 ns - 12", "lp");
   	legend->AddEntry(gr25ns16, "25 ns - 16", "lp");
   	legend->AddEntry(gr50ns8, "50 ns - 8", "lp");
   	legend->AddEntry(gr50ns12, "50 ns - 12", "lp");
   	legend->AddEntry(gr50ns16, "50 ns - 16", "lp");
   	legend->Draw();


   	return 0;


}