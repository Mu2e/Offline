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

void TimePeakTest (TTree* tp, char* page="hits" ) {
  TString spage(page);
  TCut convevent("nchits>19");
  TCut anyconvpeak("ncphits>0");
  TCut convpeak("ncphits/nchits>0.5");
  TCut convhit("_mcgen==2");
  TCut nonconvhit("_mcgen!=2");
  TCut hittime("abs(_dt)<30.0");
  
  if(spage == "hits") {
    gStyle->SetOptStat(0);
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

    tp->Project("dtc","_dt",convevent+convpeak+convhit+hittime);
    tp->Project("dtnc","_dt",convevent+convpeak+nonconvhit+hittime);
    tp->Project("dfc","_dphi",convevent+convpeak+convhit+hittime);
    tp->Project("dfnc","_dphi",convevent+convpeak+nonconvhit+hittime);
    tp->Project("rhoc","_rho",convevent+convpeak+convhit+hittime);
    tp->Project("rhonc","_rho",convevent+convpeak+nonconvhit+hittime);
    tp->Project("mvac","_mva",convevent+convpeak+convhit+hittime);
    tp->Project("mvanc","_mva",convevent+convpeak+nonconvhit+hittime);

    TCanvas* ccan = new TCanvas("ccan","Peak Hits",800,800);
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

    TLegend* leg = new TLegend(0.1,0.6,0.6,0.9);
    leg->AddEntry(mvac,"CE hits","L");
    leg->AddEntry(mvanc,"Non-CE hits","L");
    double hitpur = dtc->GetEntries()/(dtc->GetEntries()+dtnc->GetEntries());
    char label[50];
    snprintf(label,50,"Peak CE purity =%3.3f",hitpur);
    leg->AddEntry((TObject*)0,label,"");
    leg->Draw();

  } else if (spage == "peak"){
    gStyle->SetOptStat(111111);
    TH1F* pphires = new TH1F("pphires","Time Peak #phi resolution;#Delta #phi",100,-1.0,1.0);
    pphires->SetStats(1);
    TH1F* ptimeres = new TH1F("ptimeres","Time Peak time resolution;#Delta t (nsec)",100,-20.0,20.0);
    ptimeres->SetStats(1);
    TH2F* pfrac = new TH2F("pfrac","Time Peak Hit Fractions;Peak CE hits/Total CE hits (efficiency);Peak CE hits/Total Peak Hits (purity)",50,0.0,1.05,50,0.0,1.05);
    pfrac->SetStats(0);
    tp->Project("ptimeres","ptime-ctime",convevent+convpeak);
    tp->Project("pphires","pphi-cphi",convevent+convpeak);
    tp->Project("pfrac","ncphits/nphits:ncphits/nchits",convevent+anyconvpeak);
    TCanvas* pcan = new TCanvas("pcan","Peak Properties",800,800);
    pcan->Divide(2,2);
    pcan->cd(1);
    ptimeres->Draw();
    pcan->cd(2);
    pphires->Draw();
    pcan->cd(3);
    pfrac->Draw("colorz");
  }

}


