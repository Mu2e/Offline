#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

Double_t splitgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Double_t core;
  Double_t tail;
  Float_t xval = x[0];
  if(xval > par[1]) {
    core = exp(-0.5*pow((xval-par[1])/par[2],2))/par[2];
    tail = par[4]*exp(-0.5*pow((xval-par[1])/par[5],2))/par[5];
  } else {
    core = exp(-0.5*pow((xval-par[1])/par[3],2))/par[3];
    tail = (1/par[2]-1/par[3]+par[4]/par[5])*exp(-0.5*pow((xval-par[1])/par[6],2));
  }
  retval = par[0]*0.398942*(core+tail);
// add a tail Gaussian
  return retval;
}

TCut mcsel("mcent.mom>100&&mcent.td<1.0&&mcent.td>0.5774&&mc.ngood>=20");
TCut helix("helixfail==0");
TCut seed("seedfail==0");

void KalTest (TTree* trk) {
  TCanvas* kcan = new TCanvas("kcan","Kalman Fit",1200,800);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TH1F* t0 = new TH1F("t0","t0 resolution;nsec",100,-5,5);
  trk->Project("t0","t0-mcmid.t0","kalfail==0&&nactive>=20");
  TH1F* nh = new TH1F("nh","N hits",66,-0.5,65.5);
  TH1F* na = new TH1F("na","N hits",66,-0.5,65.5);
  TH1F* nd = new TH1F("nd","N hits",66,-0.5,65.5);
  nh->SetLineColor(kBlack);
  na->SetLineColor(kRed);
  nd->SetLineColor(kGreen);
  TH1F* fstat = new TH1F("fstat","fit status",22,-1.5,20.5);
  TH1F* momr = new TH1F("momr","momentum resolution at start of tracker;MeV",100,-2.5,2.5);
  TH2F* mom = new TH2F("mom","momentum at start of tracker;true momentum (MeV);fit momentum (MeV)",
    100,90,107,100,90,107);
  mom->SetStats(0);
  trk->Project("fstat","fit.status","mcent.mom>100");
  trk->Project("nh","nhits","kalfail==0");
  trk->Project("na","nactive","kalfail==0");
  trk->Project("nd","nhits-nactive","kalfail==0");
//  trk->Project("momr","fit.mom-mcmom","kalfail==0");
//  trk->Project("momr","fit.mom-mcmom","kalfail==0&&nactive>=20&&fit.mom>100&&t0err<1&&chisq/ndof<5&&fit.momerr<0.2");
  trk->Project("momr","fit.mom-mcent.mom","mcent.mom>100&&kalfail==0");
  //"&&t0err<0.8&&fit.momerr<0.1&&chisq/ndof<2");
  kcan->Clear();
  kcan->Divide(2,2);
  kcan->cd(1);
  nd->Draw();
  na->Draw("same");
  nh->Draw("same");
  TLegend* leg = new TLegend(0.3,0.7,0.8,0.9);
  leg->AddEntry(nh,"All hits","l");
  leg->AddEntry(na,"Active hits","l");
  leg->AddEntry(nd,"Disabled hits","l");
  leg->Draw();
  kcan->cd(2);
  t0->SetStats(1);
  t0->Fit("gaus");
  kcan->cd(3);
  trk->Draw("fit.mom:mcent.mom>>mom","kalfail==0&&nactive>=20");
  kcan->cd(4);
  momr->SetStats(1);
  double integral = momr->GetEntries()*momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,momr->GetRMS(),momr->GetRMS(),0.01,2*momr->GetRMS(),2*momr->GetRMS());
  sgau->SetParLimits(5,1.0*momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.8);
  momr->Fit("sgau","L");
  
  double keff = 0.5*momr->GetEntries()/fstat->GetEntries();
  TPaveText* text = new TPaveText(0.1,0.7,0.4,0.8,"NDC");  
  char line[40];
  sprintf(line,"P>100MeV/c Eff. = %3.2f",keff);
  text->AddText(line);
  text->Draw();

  }  
  void CircleTest (TTree* trk) {
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111111);

  TCanvas* ccanxy = new TCanvas("ccanxy","Circle Fit",1200,800);
  TH2F* cx = new TH2F("cx","Circle X center;MC true (mm);Fit (mm)",100,-400,400,100,-400,400);
  TH2F* cy = new TH2F("cy","Circle Y center;MC true (mm);Fit (mm)",100,-400,400,100,-400,400);
  TH2F* cr = new TH2F("cr","Circle radius;MC true (mm);Fit (mm)",100,180,350,100,180,350);
  TH1F* cxr = new TH1F("cxr","X center resolution;mm",100,-200,200);
  TH1F* cyr = new TH1F("cyr","Y center resolution;mm",100,-200,200);
  TH1F* crr = new TH1F("crr","Radius resolution;mm",100,-200,200);
  cx->SetStats(0);
  cy->SetStats(0);
  cr->SetStats(0);

  trk->Project("cxr","hcx-mccx",helix+mcsel);
  trk->Project("cyr","hcy-mccy",helix+mcsel);
  trk->Project("crr","hr-mcr",helix+mcsel);
  ccanxy->Clear();
  ccanxy->Divide(3,2);
  ccanxy->cd(1);
  trk->Draw("hcx:mccx>>cx",helix+mcsel);
  ccanxy->cd(2);
  trk->Draw("hcy:mccy>>cy",helix+mcsel);
  ccanxy->cd(3);
  trk->Draw("hr:mcr>>cr",helix+mcsel);
  ccanxy->cd(4);
  cxr->Fit("gaus");
  ccanxy->cd(5);
  cyr->Fit("gaus");
  ccanxy->cd(6);
  crr->Fit("gaus");
  
  TCanvas* ccanrz = new TCanvas("ccanrz","RZ parameters",1200,800);
  
  TH2F* dfdz = new TH2F("dfdz","Helix pitch (d#phi/dZ);MC true (radians/mm);helix fit (radians/mm)",100,0.0035,0.0065,100,0.0035,0.0065);
  TH2F* fz0 = new TH2F("fz0","Heliz #phi intercept;MC true (radians); helix fit (radians)",100,-3.5,3.5,100,-3.5,3.5);
  TH1F* dfdzr = new TH1F("dfdzr","d#phi/dZ resolution;radians/mm",100,-0.0005,0.0005);
  TH1F* fz0r = new TH1F("fz0r","#phi intercept resolution;radians",100,-0.4,0.4);
  trk->Project("dfdzr","hdfdz-mcdfdz",helix+mcsel);
  trk->Project("fz0r","hfz0-mcfz0",helix+mcsel);
  dfdz->SetStats(0);
  fz0->SetStats(0);
  
  ccanrz->Clear();
  ccanrz->Divide(2,2);
  ccanrz->cd(1);
  trk->Draw("hdfdz:mcdfdz>>dfdz",helix+mcsel);
  ccanrz->cd(2);
  trk->Draw("hfz0:mcfz0>>fz0",helix+mcsel);
  ccanrz->cd(3);
  dfdzr->Fit("gaus");
  ccanrz->cd(4);
  fz0r->Fit("gaus");
 } 

  void HelixTest (TTree* trk) {
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111111);

  TCanvas* hpcan = new TCanvas("hpcan","Helix params",1200,800);
  
  TH2F* d0 = new TH2F("d0","d0;MC true (mm);helix fit (mm)",50,-200,200,50,-200,200);
  TH2F* phi0 = new TH2F("phi0","#phi0;MC true (radians);helix fit (radians)",50,-3.5,3.5,50,-3.5,3.5);
  TH2F* om = new TH2F("om","#omega;MC true (1/mm);helix fit (1/mm)",50,0.002,0.005,50,0.002,0.005);
  TH2F* z0 = new TH2F("z0","z0;MC true (mm);helix fit (mm)",50,-1000,1000,50,-1000,1000);
  TH2F* td = new TH2F("td","tan(#lambda);MC true;helix fit",50,0.5,1.2,50,0.5,1.2);
  d0->SetStats(0);
  phi0->SetStats(0);
  om->SetStats(0);
  z0->SetStats(0);
  td->SetStats(0);


  
  hpcan->Clear();
  hpcan->Divide(3,2);
  hpcan->cd(1);
  trk->Draw("hd0:mcmid.d0>>d0",helix+mcsel);
  hpcan->cd(2);
  trk->Draw("hp0:mcmid.p0>>phi0",helix+mcsel);
  hpcan->cd(3);
  trk->Draw("hom:mcmid.om>>om",helix+mcsel);
  hpcan->cd(4);
  trk->Draw("hz0:mcmid.z0>>z0",helix+mcsel);
  hpcan->cd(5);
  trk->Draw("htd:mcmid.td>>td",helix+mcsel);

  
  TCanvas* hprcan = new TCanvas("hprcan","Helix parameter Resolution",1200,800);
  
  TH1F* d0r = new TH1F("d0r","seed d0 resolution;mm",100,-300,300);
  TH1F* phi0r = new TH1F("phi0r","seed #phi0 resolution;radians",100,-0.25,0.25);
  TH1F* omr = new TH1F("omr","seed #omega resolution;1/mm",100,-0.002,0.002);
  TH1F* z0r = new TH1F("z0r","seed z0 resolution;mm",100,-150,150);
  TH1F* tdr = new TH1F("tdr","seed tan(#lambda) resolution",100,-0.3,0.3);

  trk->Project("d0r","hd0-mcmid.d0",helix+mcsel);
  trk->Project("phi0r","hp0-mcmid.p0",helix+mcsel);
  trk->Project("omr","hom-mcmid.om",helix+mcsel);
  trk->Project("z0r","hz0-mcmid.z0",helix+mcsel);
  trk->Project("tdr","htd-mcmid.td",helix+mcsel);
  
  hprcan->Clear();
  hprcan->Divide(3,2);
  hprcan->cd(1);
  d0r->Fit("gaus");
  hprcan->cd(2);
  phi0r->Fit("gaus");
  hprcan->cd(3);
  omr->Fit("gaus");
  hprcan->cd(4);
  z0r->Fit("gaus");
  hprcan->cd(5);
  tdr->Fit("gaus");
 

}

void AntiMomRes(TTree* trk) {
  TCanvas* amcan = new TCanvas("amcan","Momentum",1200,800);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
    TCut tsel = mcsel +TCut("kalfail==0");
// selection cuts

  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trk->Project("effnorm","mcent.mom",mcsel);
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TH1F* momres[4];
// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = "nactive<20";
  ncuts[1] = "nactive<20";
  ncuts[2] = "nactive<25";
  ncuts[3] = "nactive<30";
  t0cuts[0] = "";
  t0cuts[1] = "t0err>1.5";
  t0cuts[2] = "t0err>1.0";
  t0cuts[3] = "t0err>0.9";
  momcuts[0] = "";
  momcuts[1] = "fit.momerr>0.2";
  momcuts[2] = "fit.momerr>0.18";
  momcuts[3] = "fit.momerr>0.15";
  fitcuts[0] = "";
  fitcuts[1] = "fitcon<1e-4";
  fitcuts[2] = "fitcon<1e-3";
  fitcuts[3] = "fitcon<1e-2";

  amcan->Clear();
  amcan->Divide(2,2);
  for(unsigned ires=0;ires<4;ires++){
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",151,-2.5,2.5);
    TCut total = ncuts[ires] || t0cuts[ires] || momcuts[ires] || fitcuts[ires];
    trk->Project(mname,"fit.mom-mcent.mom",total+tsel);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    if(integral > 0){
      sgau->SetParameters(integral,0.0,momres[ires]->GetRMS(),momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
      sgau->SetParLimits(5,1.0*momres[ires]->GetRMS(),1.0);
      sgau->SetParLimits(6,1.0*momres[ires]->GetRMS(),1.0);
      sgau->SetParLimits(4,0.0,0.8);
      amcan->cd(ires+1);
      momres[ires]->Fit("sgau","L");

      double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
      TPaveText* text = new TPaveText(0.1,0.4,0.4,0.8,"NDC");  
      char line[40];
      sprintf(line,"%s",ncuts[ires].GetTitle());
      text->AddText(line);
      sprintf(line,"%s",t0cuts[ires].GetTitle());
      text->AddText(line);
      sprintf(line,"%s",momcuts[ires].GetTitle());
      text->AddText(line);
      sprintf(line,"%s",fitcuts[ires].GetTitle());
      text->AddText(line);
      sprintf(line,"Eff=%3.2f",keff);
      text->AddText(line);

      text->Draw();
    }
    amcan->cd(0);
  }  

}


void MomRes(TTree* trk) {
  TCanvas* mcan = new TCanvas("mcan","Momentum",1200,800);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
// should have pitch angle and generated hit cuts here, FIXME!!!
  TCut tsel = mcsel +TCut("kalfail==0");
// selection cuts
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trk->Project("effnorm","mcent.mom",mcsel);
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TH1F* momres[4];
// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = "nactive>=20";
  ncuts[1] = "nactive>=20";
  ncuts[2] = "nactive>=25";
  ncuts[3] = "nactive>=30";
  t0cuts[0] = "";
  t0cuts[1] = "t0err<1.5";
  t0cuts[2] = "t0err<1.0";
  t0cuts[3] = "t0err<0.9";
  momcuts[0] = "";
  momcuts[1] = "fit.momerr<0.2";
  momcuts[2] = "fit.momerr<0.18";
  momcuts[3] = "fit.momerr<0.15";
  fitcuts[0] = "";
  fitcuts[1] = "fitcon>1e-4";
  fitcuts[2] = "fitcon>1e-3";
  fitcuts[3] = "fitcon>1e-2";

  mcan->Clear();
  mcan->Divide(2,2);
  for(unsigned ires=0;ires<4;ires++){
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",151,-2.5,2.5);
    trk->Project(mname,"fit.mom-mcent.mom",ncuts[ires]+t0cuts[ires]+momcuts[ires]+fitcuts[ires]+tsel);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momres[ires]->GetRMS(),momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
    sgau->SetParLimits(5,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.8);
    mcan->cd(ires+1);
    momres[ires]->Fit("sgau","L");

    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
    TPaveText* text = new TPaveText(0.1,0.4,0.4,0.8,"NDC");  
    char line[40];
    sprintf(line,"%s",ncuts[ires].GetTitle());
    text->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    text->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    text->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    text->AddText(line);
    sprintf(line,"Eff=%3.2f",keff);
    text->AddText(line);
 
    text->Draw();
  } 
  mcan->cd(0); 
}
  
void SeedTest (TTree* trk) {
  
  TCanvas* scan = new TCanvas("srcan","Seed parameter Resolution",1200,800);
  
  TH1F* d0s = new TH1F("d0s","seed d0 resolution;mm",100,-100,100);
  TH1F* phi0s = new TH1F("phi0s","seed #phi0 resolution;radians",100,-0.07,0.07);
  TH1F* oms = new TH1F("oms","seed #omega resolution;1/mm",100,-0.0005,0.0005);
  TH1F* z0s = new TH1F("z0s","seed z0 resolution;mm",100,-30,30);
  TH1F* tds = new TH1F("tds","seed tan(#lambda) resolution",100,-0.07,0.07);


  trk->Project("d0s","sd0-mcmid.d0",seed+mcsel);
  trk->Project("phi0s","sp0-mcmid.p0",seed+mcsel);
  trk->Project("oms","som-mcmid.om",seed+mcsel);
  trk->Project("z0s","sz0-mcmid.z0",seed+mcsel);
  trk->Project("tds","std-mcmid.td",seed+mcsel);
  
  scan->Clear();
  scan->Divide(3,2);
  scan->cd(1);
  d0s->Fit("gaus");
  scan->cd(2);
  phi0s->Fit("gaus");
  scan->cd(3);
  oms->Fit("gaus");
  scan->cd(4);
  z0s->Fit("gaus");
  scan->cd(5);
  tds->Fit("gaus");
  
}
