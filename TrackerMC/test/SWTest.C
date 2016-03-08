#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLine.h"
#include "TProfile.h"
#include "TFitResult.h"


void SWTest(TTree* sw){
  THStack* vstack = new THStack("vm","Maximum Voltage at Discriminator;Vmax (mV)");
  TH1F* vmo = new TH1F("vmo","Vmax (x20 gain);Vmax (mV)",201,-0.5,200.5);
  vmo->SetFillColor(kBlue);
  sw->Project("vmo","vmax","nxing>0&&(mce<100 || npart>1)");
  vstack->Add(vmo);
  TH1F* vmce = new TH1F("vmce","Vmax (x20 gain);Vmax (mV)",201,-0.5,200.5);
  vmce->SetLineColor(kRed);
  sw->Project("vmce","vmax","nxing>0&&(mce>100 && npart==1)");
  vstack->Add(vmce);
  TH1F* npart = new TH1F("npart","N particles",5,-0.5,4.5);
  sw->Project("npart","npart");
  THStack* estack = new THStack("estack","G4 Energy in Straw Hit;Energy (KeV)");
  TH1F* g4e0x = new TH1F("g4e0x","G4 Energy in Straw Hit;Energy (KeV)",150,0,4.0);
  g4e0x->SetFillColor(kRed);
  TH1F* g4e1x = new TH1F("g4e1x","G4 Energy in Straw Hit;Energy (KeV)",150,0,4.0);
  g4e1x->SetFillColor(kCyan);
  TH1F* g4e2x = new TH1F("g4e2x","G4 Energy in Straw Hit;Energy (KeV)",150,0,4.0);
  g4e2x->SetFillColor(kBlue);
  sw->Project("g4e0x","sesum*1000.0","tmin>300&&nxing==0");
  estack->Add(g4e0x);
  sw->Project("g4e1x","sesum*1000.0","tmin>300&&nxing==1");
  estack->Add(g4e1x);
  sw->Project("g4e2x","sesum*1000.0","tmin>300&&nxing>=2");
  estack->Add(g4e2x);
//
  THStack* dstack = new THStack("dstack","First Cluster Drift Distance;D (mm)");
  TH1F* dw0x = new TH1F("dw0x","First Cluster Drift Distance;D (mm)",150,0.0,2.6);
  dw0x->SetFillColor(kRed);
  TH1F* dw1x = new TH1F("dw1x","First Cluster Drift Distance;D (mm)",150,0.0,2.6);
  dw1x->SetFillColor(kCyan);
  TH1F* dw2x = new TH1F("dw2x","First Cluster Drift Distance;D (mm)",150,0.0,2.6);
  dw2x->SetFillColor(kBlue);
  sw->Project("dw0x","xddist","tmin>300&&nxing==0");
  dstack->Add(dw0x);
  sw->Project("dw1x","xddist","tmin>300&&nxing==1");
  dstack->Add(dw1x);
  sw->Project("dw2x","xddist","tmin>300&&nxing>=2");
  dstack->Add(dw2x);
//  TH1F* nclus = new TH1F("nclus","N clusters",50,-0.5,49.5);
  TH2F* nclusvs = new TH2F("nclusvs","N clusters vs steplength;steplength (mm);N clusters",50,0,10.0,50,-0.5,49.5);
  nclusvs->SetStats(0);
//  sw->Project("nclus","nhitlet","npAart==1&&nstep==1");
  sw->Project("nclusvs","nhitlet:slen","npart==1&&nstep==1");
  TH1F* tvmax = new TH1F("tvmax","Time of maximum voltage;ns",100,0,50);
  tvmax->SetStats(0);
  sw->Project("tvmax","tvmax-tmin","nxing>0");

  TH2F* vvsdist = new TH2F("vvsdist","Voltage vs Drift Distance",50,0,2.51,50,0,400);
  sw->Project("vvsdist","vmax:xddist","nxing>0");

  TCanvas* swcan = new TCanvas("swcan","swcan",1200,800);
  swcan->Divide(2,2);
  swcan->cd(1);
  vstack->Draw();
  TLegend * vleg = new TLegend(0.5,0.7,0.9,0.9);
  TLine* tc = new TLine(12.0,0.0,12.0,vstack->GetMaximum());
  tc->SetLineColor(kBlack);
  tc->SetLineWidth(3.0);
  tc->Draw();
  vmce->SetStats(1);
  TFitResultPtr gfit = vmce->Fit("gaus","srq0","same",20,180);
  char label[100];
  snprintf(label,100,"CE only in straw, mean =%4.1f",gfit->Parameter(1));
  vleg->AddEntry(vmce,label,"F");
  vleg->AddEntry(vmo,"#delta-ray and/or CE in straw","F");
  vleg->AddEntry(tc,"Threshold","L");
  vleg->Draw();
  swcan->cd(2);
//  npart->Draw();
//  vvsdist->Draw();
//  swcan->cd(3);
//  nclusvs->Draw("box");
  tvmax->Draw();
  swcan->cd(3);
  estack->Draw();
  TLegend * eleg = new TLegend(0.5,0.7,0.9,0.9);
  eleg->AddEntry(g4e0x,"0 threshold crossing","F");
  eleg->AddEntry(g4e1x,"1 Threshold crossing","F");
  eleg->AddEntry(g4e2x,"#geq 2 Threshold crossings","F");
  eleg->Draw();
  swcan->cd(4);
  dstack->Draw();
}
