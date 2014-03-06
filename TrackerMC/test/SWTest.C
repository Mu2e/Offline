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

//  TH1F* nclus = new TH1F("nclus","N clusters",50,-0.5,49.5);
  TH2F* nclusvs = new TH2F("nclusvs","N clusters vs steplength;steplength (mm);N clusters",50,0,10.0,25,-0.5,24.5);
  nclusvs->SetStats(0);
//  sw->Project("nclus","nhitlet","npAart==1&&nstep==1");
  sw->Project("nclusvs","nhitlet:slen","npart==1&&nstep==1");

  TCanvas* can = new TCanvas("vcan","vcan",1200,800);
  can->Divide(2,2);
  can->cd(1);
  vstack->Draw();
  TLegend * leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(vmce,"CE only in straw","F");
  leg->AddEntry(vmo,"#delta-ray and/or CE in straw","F");
  leg->Draw();
  can->cd(2);
  npart->Draw();
  can->cd(3);
  nclusvs->Draw("box");
  can->cd(4);
  
}
