#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"


void dioacc(TTree* dio, TCut acut) {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1F* amomd = new TH1F("amomd","DIO Acceptance vs Momentum;MeV/c",500,95,105);
  TH1F* amom = new TH1F("amom","DIO Acceptance vs Momentum;MeV/c",500,95,105);
  amom->Sumw2();
  amomd->Sumw2();
  dio->Project("amom","mcmom",acut);
  dio->Project("amomd","mcmom");
  amom->Divide(amomd);

//  TH1F* apitch = new TH1F("apitch","DIO acceptance vs tan#lambda;tan#lambda",200,0.5,1.1);
//  TH1F* apitchd = new TH1F("apitchd","DIO acceptance vs tan#lambda;tan#lambda",200,0.5,1.1);
//  apitch->Sumw2();
//  apitchd->Sumw2();
//  dio->Project("apitchd","mcent.td");
//  dio->Project("apitch","mcent.td",acut);
//  apitch->Divide(apitchd);
//
//  TH1F* aradius = new TH1F("aradius","DIO acceptance vs curvature radius;#rho (mm)",200,200,310);
//  TH1F* aradiusd = new TH1F("aradiusd","DIO acceptance vs curvature radius;#rho (mm)",200,200,310);
//  aradius->Sumw2();
//  aradiusd->Sumw2();
//  dio->Project("aradiusd","1.0/mcent.om");
//  dio->Project("aradius","1.0/mcent.om",acut);
//  aradius->Divide(aradiusd);

//  TH1F* ad0 = new TH1F("ad0","DIO acceptance vs d0;d0 (mm)",200,-100,120);
//  TH1F* ad0d = new TH1F("ad0d","DIO acceptance vs d0;d0 (mm)",200,-100,120);
//  ad0->Sumw2();
//  ad0d->Sumw2();
//  dio->Project("ad0d","mcent.d0");
//  dio->Project("ad0","mcent.d0",acut);
//  ad0->Divide(ad0d);

  TH1F* acost = new TH1F("acost","DIO acceptance vs cos(#theta);cos(#theta)",200,-1.0,1.0);
  TH1F* acostd = new TH1F("acostd","DIO acceptance vs cos(#theta);cos(#theta)",200,-1.0,1.0);
  acost->Sumw2();
  acostd->Sumw2();
  dio->Project("acostd","mctd/sqrt(1+mctd^2)");
  dio->Project("acost","mctd/sqrt(1+mctd^2)",acut);
  acost->Divide(acostd);

  TH1F* aphi = new TH1F("aphi","DIO acceptance vs #phi;#phi (mm)",200,-3.1416,3.1416);
  TH1F* aphid = new TH1F("aphid","DIO acceptance vs #phi;#phi (mm)",200,-3.1416,3.1416);
  aphi->Sumw2();
  aphid->Sumw2();
  dio->Project("aphid","mcent.p0");
  dio->Project("aphi","mcent.p0",acut);
  aphi->Divide(aphid);

  TH1F* atargetr = new TH1F("atargetr","DIO acceptance vs target radius;target #rho (mm)",200,0.0,90.0);
  TH1F* atargetrd = new TH1F("atargetrd","DIO acceptance vs target radius;target #rho (mm)",200,0.0,90.0);
  atargetr->Sumw2();
  atargetrd->Sumw2();
  dio->Project("atargetrd","sqrt(mcx^2+mcy^2)");
  dio->Project("atargetr","sqrt(mcx^2+mcy^2)",acut);
  atargetr->Divide(atargetrd);

  TH1F* atargetz = new TH1F("atargetz","DIO acceptance vs target z;target z (mm)",34,-4765,-3905);
  TH1F* atargetzd = new TH1F("atargetzd","DIO acceptance vs target z;target z (mm)",34,-4765,-3905);
  atargetz->Sumw2();
  atargetzd->Sumw2();
  dio->Project("atargetzd","mcz");
  dio->Project("atargetz","mcz",acut);
  atargetz->Divide(atargetzd);

//  TH1F* armax = new TH1F("armax","DIO acceptance vs rmax;rmax (mm)",200,4250,680);
//  TH1F* armaxd = new TH1F("armaxd","DIO acceptance vs rmax;rmax (mm)",200,450,680);
//  armax->Sumw2();
//  armaxd->Sumw2();
//  dio->Project("armaxd","mcent.d0+2.0/mcent.om");
//  dio->Project("armax","mcent.d0+2.0/mcent.om",acut);
//  armax->Divide(armaxd);

  TCanvas* acan = new TCanvas("acan","acceptance",1200,800);
  TCanvas* acan_2 = new TCanvas("acan_2","acceptance",1200,800);
  TCanvas* acan_3 = new TCanvas("acan_3","acceptance",1200,800);
  acan->Divide(1,3);
  acan->cd(1);
  amom->SetMaximum(0.22);
  amom->Fit("pol1");
  acan_2->Divide(2,1);
  acan_2->cd(1);
  acost->Draw();
  acan_2->cd(2);
  aphi->SetMinimum(0.12);
  aphi->Draw();
  acan_3->Divide(2,1);
  acan_3->cd(1);
  atargetr->SetMinimum(0.1);
  atargetr->Draw();
  acan_3->cd(2);
  atargetz->SetMinimum(0.12);
  atargetz->Draw();


}
