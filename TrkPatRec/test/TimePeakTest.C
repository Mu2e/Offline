#include <algorithm>
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void TimePeakTest (TTree* tp, char* page="cleaning" ) {
  TString spage(page);
  gStyle->SetOptStat(0);
  TCut convpeak("ncphits>9");
  TCut convhit("_mcgen==2");
  TCut nonconvhit("_mcgen!=2");
  
  if(spage == "cleaning") {
    TH1F* dtc = new TH1F("dtc","Peak hit #Delta t",100,-40.0,40.0);
    TH1F* dtnc = new TH1F("dtnc","Peak hit #Delta t",100,-40.0,40.0);
    dtc->SetLineColor(kRed);
    dtnc->SetLineColor(kBlue);
    TH1F* dfc = new TH1F("dfc","Peak hit #Delta #phi",100,-4,4);
    TH1F* dfnc = new TH1F("dfnc","Peak hit #Delta #phi",100,-4,4);
    dfc->SetLineColor(kRed);
    dfnc->SetLineColor(kBlue);
    TH1F* rhoc = new TH1F("rhoc","Peak hit #rho",100,380.0,680.0);
    TH1F* rhonc = new TH1F("rhonc","Peak hit #rho",100,380.0,680.0);
    rhoc->SetLineColor(kRed);
    rhonc->SetLineColor(kBlue);
    TH1F* mvac = new TH1F("mvac","Peak hit MVA",100,-0.1,1.1);
    TH1F* mvanc = new TH1F("mvanc","Peak hit MVA",100,-0.1,1.1);
    mvac->SetLineColor(kRed);
    mvanc->SetLineColor(kBlue);

    tp->Project("dtc","_dt",convpeak+convhit);
    tp->Project("dtnc","_dt",convpeak+nonconvhit);
    tp->Project("dfc","_dphi",convpeak+convhit);
    tp->Project("dfnc","_dphi",convpeak+nonconvhit);
    tp->Project("rhoc","_rho",convpeak+convhit);
    tp->Project("rhonc","_rho",convpeak+nonconvhit);
    tp->Project("mvac","_mva",convpeak+convhit);
    tp->Project("mvanc","_mva",convpeak+nonconvhit);

    TCanvas* ccan = new TCanvas("ccan","Peak Cleaning",800,800);
    ccan->Divide(2,2);
    ccan->cd(1);
    dtc->Draw();
    dtnc->Draw("same");
    ccan->cd(2);
    dfc->Draw();
    dfnc->Draw("same");
    ccan->cd(3);
    rhoc->Draw();
    rhonc->Draw("same");
    ccan->cd(4);
    mvac->Draw();
    mvanc->Draw("same");

  }

}


