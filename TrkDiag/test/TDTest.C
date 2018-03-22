#include "TTree.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;

Double_t ALine(Double_t *x, Double_t *par) {
  if(x[0] > par[0])
    return par[1];
  else
    return par[1]+(x[0]-par[0])*par[2];
}
 
void TDTest(TTree* sh,const char* page="edep") {
  std::string spage(page);
  TCut conv("mcproc==56&&mcgen==2");
  TCut proton("mcpdg==2212");
  TCut other("mcpdg==11&&mcgen!=2");
  if(spage == "edep"){
    TH1F* cedep = new TH1F("cedep","Ce Straw Energy Deposit;Edep (KeV)",50,0,8.0);
    TH1F* pedep = new TH1F("pedep","Proton Straw Energy Deposit;Edep (KeV)",50,0,8.0);
    cedep->SetStats(0);
    pedep->SetStats(0);
    cedep->SetLineColor(kRed);
    pedep->SetLineColor(kBlue);
    sh->Project("cedep","edep*1000.0",conv);
    sh->Project("pedep","edep*1000.0",proton);
    cedep->Scale(pedep->GetEntries()/cedep->GetEntries());
    TCanvas* edep = new TCanvas("edep","edep",400,400);
    cedep->Draw("Hist");
    pedep->Draw("same");
    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(cedep,"Conversion electron","l");
    leg->AddEntry(pedep,"Proton","l");
    leg->Draw();

  } else if(spage == "edlen"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH2F* edslen = new TH2F("edslen","Ce Reco Hit Energy vs Straw Half Length;Straw Half Length (mm);EDep (KeV)",40,200.0,600.0,50,0,6.0);
    TProfile* edslenp = new TProfile("edslenp","Ce Reco Hit Energy vs Straw Half Length;Straw Half Length (mm);EDep (KeV)",40,200.0,600.0,0.0,3.0,"S");
    edslen->SetStats(0);
    sh->Project("edslen","1000.0*edep:slen",conv);
    sh->Project("edslenp","1000.0*edep:slen",conv);

    edslen->FitSlicesY(0,0,-1,20);

    TCanvas* ed = new TCanvas("ed","ed",600,600);
    TH1D *edslen_1 = (TH1D*)gDirectory->Get("edslen_1");
    edslen_1->SetLineColor(kRed);
    edslen_1->Fit("pol1");
    edslen_1->SetStats(1);
    edslen->Draw("colorz");
    edslen_1->Draw("sames");
    edslenp->Draw("same");
  } else if(spage == "TDlen"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH2F* ctdlen = new TH2F("ctdlen","Ce TDRes vs reco wire position;reco wire position (mm);reco - mc wire position (mm)",20,0.0,400.0,50,-100.0,100.0);
    TH2F* ptdlen = new TH2F("ptdlen","Proton TDRes vs reco wire position;reco wire position (mm);reco - mc wire position (mm)",20,0.0,400.0,50,-80.0,80.0);
    sh->Project("ctdlen","shlen-mcshlen:abs(shlen)",conv);
    sh->Project("ptdlen","shlen-mcshlen:abs(shlen)",proton);
    ctdlen->FitSlicesY(0,0,-1);
    ptdlen->FitSlicesY(0,0,-1);
    TCanvas* tdlcan = new TCanvas("tdlcan","tdlcan",700,700);
    tdlcan->Divide(2,2);
    tdlcan->cd(1);
    ctdlen->Draw("colorz");
    tdlcan->cd(2);
    ptdlen->Draw("colorz");
    tdlcan->cd(3);
    TH1D *ctdlen_2 = (TH1D*)gDirectory->Get("ctdlen_2");
    ctdlen_2->SetLineColor(kRed);
    ctdlen_2->Fit("pol3");
    ctdlen_2->SetStats(1);
    tdlcan->cd(4);
    TH1D *ptdlen_2 = (TH1D*)gDirectory->Get("ptdlen_2");
    ptdlen_2->SetLineColor(kRed);
    ptdlen_2->Fit("pol3");
    ptdlen_2->SetStats(1);
    TH2F* ctdlenp = new TH2F("ctdlenp","Ce TD pull vs reco wire position;reco wire position (mm);TD pull",15,0.0,450.0,50,-6.0,6.0);
    TH2F* ptdlenp = new TH2F("ptdlenp","Proton TD pull vs reco wire position;reco wire position (mm);TD pull",15,0.0,450.0,50,-6.0,6.0);
    sh->Project("ctdlenp","(shlen-mcshlen)/wres:abs(shlen)",conv);
    sh->Project("ptdlenp","(shlen-mcshlen)/wres:abs(shlen)",proton);
    ctdlenp->FitSlicesY(0,0,-1);
    ptdlenp->FitSlicesY(0,0,-1);
    TCanvas* tdlpcan = new TCanvas("tdlpcan","tdlpcan",700,700);
    tdlpcan->Divide(2,2);
    tdlpcan->cd(1);
    ctdlenp->Draw("colorz");
    tdlpcan->cd(2);
    ptdlenp->Draw("colorz");
    tdlpcan->cd(3);
    TH1D *ctdlenp_2 = (TH1D*)gDirectory->Get("ctdlenp_2");
    ctdlenp_2->SetLineColor(kRed);
    ctdlenp_2->Fit("pol3");
    ctdlenp_2->SetStats(1);
    tdlpcan->cd(4);
    TH1D *ptdlenp_2 = (TH1D*)gDirectory->Get("ptdlenp_2");
    ptdlenp_2->SetLineColor(kRed);
    ptdlenp_2->Fit("pol3");
    ptdlenp_2->SetStats(1);
  } else if(spage == "TDREdep"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    std::vector<TH2F*> tdrvec;
    char cut[50];
    char name[50];
    char title[80];
    unsigned nplots(0);
    float emin(0.0);
    float ebin(0.5);
    float emax = emin+ebin;
    while(emax < 7.0){
      ++nplots;
      snprintf(cut,50,"1000.0*edep>%3.1f&&1000.0*edep<%3.1f",emin,emax);
      snprintf(name,50,"tdr%3.1f",emin);
      snprintf(title,80,"TD pull vs reco wire position, %3.1fKeV < edep < %3.1fKeV;reco wire position (mm);TD pull",emin,emax);
      TH2F* hist = new TH2F(name,title,50,0.0,650.0,50,-6.0,6.0);
      tdrvec.push_back(hist);
      sh->Project(name,"(shlen-mcshlen)/wres:abs(shlen)",cut);
      emin += ebin;
      emax += ebin;
    }
    TCanvas* tdrcan = new TCanvas("tdrcan","tdrcan",700,700);
    unsigned nx = (unsigned)ceil(sqrt(nplots));
    unsigned ny = ceil(nplots/float(nx));
    tdrcan->Divide(nx,ny);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdrcan->cd(iplot+1);
      tdrvec[iplot]->Draw("colorz");
    }
    TCanvas* tdrfcan = new TCanvas("tdrfcan","tdrfcan",700,700);
    tdrfcan->Divide(nx,ny);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdrfcan->cd(iplot+1);
      tdrvec[iplot]->FitSlicesY(0,0,-1);
      std::string fname(tdrvec[iplot]->GetName());
      fname += "_2";
      TH1D *tdr_2 = (TH1D*)gDirectory->Get(fname.c_str());
      std::string title(tdrvec[iplot]->GetTitle());
      title += " Sigma";
      tdr_2->SetTitle(title.c_str());
      tdr_2->SetMinimum(0.8);
      tdr_2->SetMaximum(1.2);
      tdr_2->SetLineColor(kRed);
      tdr_2->Draw();
    }
  } else if(spage == "Proton"){
    TH2F* pcomp = new TH2F("pcomp","Proton TD Reco wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-650,650,50,-650,650);
    sh->Project("pcomp","shlen:mcshlen",proton);
    TH1F* pres = new TH1F("pres","Proton TD Wire Distance Resolution; Resolution (mm)",50,-300,300);
    sh->Project("pres","shlen-mcshlen",proton);
    TH1F* ppull = new TH1F("ppull","Proton TD Pull;Pull",100,-10,10);
    sh->Project("ppull","(shlen-mcshlen)/wres",proton);
    gStyle->SetOptStat(1111);
    TCanvas* pcan = new TCanvas("pcan","pcan",700,700);
    pcan->Divide(2,2);
    pcan->cd(1);
    pcomp->Draw("colorz");
    pcan->cd(2);
    pres->Fit("gaus","","",-2.0*pres->GetRMS(),2.0*pres->GetRMS());
    pcan->cd(3);
    ppull->Fit("gaus","","",-2.0*ppull->GetRMS(),2.0*ppull->GetRMS());
  } else if(spage == "Ce"){
    TH2F* cecomp = new TH2F("cecomp","Ce TD Reco wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-650,650,50,-650,650);
    sh->Project("cecomp","shlen:mcshlen",conv);
    TH1F* ceres = new TH1F("ceres","Ce TD Wire Distance Resolution; Resolution (mm)",50,-300,300);
    sh->Project("ceres","shlen-mcshlen",conv);
    TH1F* cepull = new TH1F("cepull","Ce TD Pull;Pull",100,-10,10);
    sh->Project("cepull","(shlen-mcshlen)/wres",conv);
    gStyle->SetOptStat(1111);
    TCanvas* cecan = new TCanvas("cecan","cecan",700,700);
    cecan->Divide(2,2);
    cecan->cd(1);
    cecomp->Draw("colorz");
    cecan->cd(2);
    ceres->Fit("gaus","","",-2.0*ceres->GetRMS(),2.0*ceres->GetRMS());
    cecan->cd(3);
    cepull->Fit("gaus","","",-2.0*cepull->GetRMS(),2.0*cepull->GetRMS());
  } else if(spage == "Other"){
    TH2F* otcomp = new TH2F("otcomp","Other TD Reco wire distanot vs true;True TD distanot (mm);Reco wire distanot (mm)",50,-650,650,50,-650,650);
    sh->Project("otcomp","shlen:mcshlen",other);
    TH1F* otres = new TH1F("otres","Other TD Wire Distanot Resolution; Resolution (mm)",50,-300,300);
    sh->Project("otres","shlen-mcshlen",other);
    TH1F* otpull = new TH1F("otpull","Other TD Pull;Pull",100,-10,10);
    sh->Project("otpull","(shlen-mcshlen)/wres",other);
    gStyle->SetOptStat(1111);
    TCanvas* otcan = new TCanvas("otcan","otcan",700,700);
    otcan->Divide(2,2);
    otcan->cd(1);
    otcomp->Draw("colorz");
    otcan->cd(2);
    otres->Fit("gaus","","",-2.0*otres->GetRMS(),2.0*otres->GetRMS());
    otcan->cd(3);
    otpull->Fit("gaus","","",-2.0*otpull->GetRMS(),2.0*otpull->GetRMS());
  }
}
