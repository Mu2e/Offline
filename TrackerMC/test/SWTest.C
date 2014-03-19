#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"


void SWTest(TTree* sw){
  THStack* vstack = new THStack("vm","Maximum Voltage at Discriminator;Vmax (mV)");
  TH1F* vmo = new TH1F("vmo","Vmax (x20 gain);Vmax (mV)",100,0,1001);
  vmo->SetFillColor(kBlue);
  sw->Project("vmo","vmax","nxing>0&&(mce<100 || npart>1)");
  vstack->Add(vmo);
  TH1F* vmce = new TH1F("vmce","Vmax (x20 gain);Vmax (mV)",100,0,1001);
  vmce->SetFillColor(kRed);
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
//  TH1F* nclus = new TH1F("nclus","N clusters",50,-0.5,49.5);
  TH2F* nclusvs = new TH2F("nclusvs","N clusters vs steplength;steplength (mm);N clusters",50,0,10.0,50,-0.5,49.5);
  nclusvs->SetStats(0);
//  sw->Project("nclus","nhitlet","npAart==1&&nstep==1");
  sw->Project("nclusvs","nhitlet:slen","npart==1&&nstep==1");

  TCanvas* swcan = new TCanvas("swcan","swcan",1200,800);
  swcan->Divide(2,2);
  swcan->cd(1);
  vstack->Draw();
  TLegend * vleg = new TLegend(0.5,0.7,0.9,0.9);
  vleg->AddEntry(vmce,"CE only in straw","F");
  vleg->AddEntry(vmo,"#delta-ray and/or CE in straw","F");
  vmce->Fit("gaus","","same",20,180);
  vleg->Draw();
  swcan->cd(2);
  npart->Draw();
  swcan->cd(3);
  nclusvs->Draw("box");
  swcan->cd(4);
  estack->Draw();
  TLegend * eleg = new TLegend(0.5,0.7,0.9,0.9);
  eleg->AddEntry(g4e0x,"0 threshold crossing","F");
  eleg->AddEntry(g4e1x,"1 Threshold crossing","F");
  eleg->AddEntry(g4e2x,"#geq 2 Threshold crossings","F");
  eleg->Draw();
}
