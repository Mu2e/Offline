#define mmAlg2_cxx
#include "mmAlg2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TProfile.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
using namespace std;

//class constructor

void mmAlg2::Loop()
{
	int N = 100;
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (dem__status > 0 && demmc_proc == 56 && dem__mom > 100.0) {
	// if (Cut(ientry) < 0) continue;
	for (int i = 0; i < demtshmc_ ; ++i){
	  if(demtsh__active[i] && demtshmc__rel[i] == 0){
	    double proc = demtshmc__proc[i];
	    double rdrift = demtsh__rdrift[i];
	    double mom = demtshmc__mom[i];
	    double doca = abs(demtsh__doca[i]);
	    double distVal = rdrift - doca;
	    h1->Fill(doca,distVal);
	    prof1->Fill(doca,distVal);

	  }
	}
      }
    }

    TCanvas*c1 = new TCanvas("c1","c1",1000,1000);
    prof1->SetMarkerSize(0.7);
    h1->Draw("colorZ");
    prof1->Draw("same");
    TCanvas*c4 = new TCanvas("c4","c4",1000,1000);
    prof1->Draw();



    for (int i=1; i<101;++i){
      double valueMean = prof1->GetBinContent(i);
      double valueRms = prof1->GetBinErrorLow(i);
      means.push_back(valueMean);
      RMSs.push_back(valueRms);
      double x1 = (2.5/100)*i;
      double x2 = (2.5/100)*(i+1);
      double x3  = (x1+x2)/2;
      slices.push_back(x3);

    }


    double *slicesArray = &slices[0];
    double *meansArray = &means[0];  //array containing the means of the 2d plot coming from TProfile
    double *RMSarray = &RMSs[0];  //array containing the RMSs of the 2d plot coming from TProfile
    gStyle->SetOptStat(111100);

    TFile* fout2 = new TFile("TProfileMeans.root","RECREATE");
    TCanvas* can1 = new TCanvas("can1","can1",1000,1000); //means of slices
    TGraph* g1 = new TGraph(100,slicesArray,meansArray);
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(.75);
    g1->SetTitle("means of slices from 2d plot");
    g1->GetXaxis()->SetTitle("center of slice, doca (mm)");
    g1->GetYaxis()->SetTitle("mean (mm)");
    g1->Draw("apl");
    g1->Write();
    fout2->Close();

    TFile* fout = new TFile("TProfileRMSs.root","RECREATE");
    TCanvas* can2 = new TCanvas("can2","can2",1000, 1000);  // rms's of slices
    TGraph* g2 = new TGraph(100,slicesArray,RMSarray);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(.75);
    g2->SetTitle("RMSs of slices from 2d plot");
    g2->GetXaxis()->SetTitle("center of slice, doca (mm)");
    g2->GetYaxis()->SetTitle("RMS (mm)");
    g2->Draw("apl");
    g2->Write();
    fout->Close();


}



