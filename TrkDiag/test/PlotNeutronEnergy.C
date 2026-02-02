#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TCut.h"
#include "TGraph.h"
#include "TLegend.h"
#include <iostream>

void neutrons(TTree* sh) {
  TCut pneut("mcgpdg==2112");
  TCut pneutcap = pneut + TCut("mcpproc==35");
  TCut pneutscat = pneut + TCut("mcpproc==36");
  TCut pneutother = pneut + TCut("mcpproc<35||mcpproc>36");

  TH1F* necap = new TH1F("necap","neutron producing tracker hit Kinetic Energy;KE (MeV)",100,0,50);
  TH1F* nescat = new TH1F("nescat","neutron producing tracker hit Kinetic Energy;KE (MeV)",100,0,50);
  TH1F* neother = new TH1F("neother","neutron producing tracker hit Kinetic Energy;KE (MeV)",100,0,50);
  necap->SetLineColor(kRed);
  nescat->SetLineColor(kBlue);
  neother->SetLineColor(kGreen);
  necap->SetStats(0);

  TH2F* npcap = new TH2F("npcap","neutron-induced #gamma origin;z (mm);#rho (mm)",50,-7000,4000,50,0,1500);
  TH2F* npscat = new TH2F("npscat","neutron-induced #gamma origin;z (mm);#rho (mm)",50,-7000,4000,50,0,1500);
  TH2F* npother = new TH2F("npother","neutron-induced #gamma origin;z (mm);#rho (mm)",50,-7000,4000,50,0,1500);

  npcap->SetLineColor(kRed);
  npscat->SetLineColor(kBlue);
  npother->SetLineColor(kGreen);
  npcap->SetStats(0);

  sh->Project("necap","mcge-939.565",pneutcap);
  sh->Project("nescat","mcge-939.565",pneutscat);
  sh->Project("neother","mcge-939.565",pneutother);

  sh->Project("npcap","sqrt(mcpopos.x^2+mcpopos.y^2):mcpopos.z",pneutcap);
  sh->Project("npscat","sqrt(mcpopos.x^2+mcpopos.y^2):mcpopos.z",pneutscat);
  sh->Project("npother","sqrt(mcpopos.x^2+mcpopos.y^2):mcpopos.z",pneutother);

  TGraph* ngen = new TGraph("Offline/EventGenerator/data/neutronSpectrum.txt");
  ngen->SetMarkerStyle(20);
  ngen->SetMarkerColor(kCyan);
  ngen->SetMarkerSize(1);
  double integral = 0.0;
  for(size_t ipt = 0;ipt<ngen->GetN();++ipt){
    double x,y;
    if(ngen->GetPoint(ipt,x,y))
      integral += y;
  }
  double scale = necap->Integral()/integral;
  cout << "scaling factor = " << scale << endl;
  for(size_t ipt = 0;ipt<ngen->GetN();++ipt){
    double x,y;
    if(ngen->GetPoint(ipt,x,y))
      ngen->SetPoint(ipt,x,y*scale);
  }

  TCanvas* ncan = new TCanvas("ncan","neutrons",1200,800);
  ncan->Divide(2,1);
  ncan->cd(1);
  necap->Draw();
  nescat->Draw("same");
  neother->Draw("same");
  ngen->Draw("P");
  TLegend* leg = new TLegend(0.3,0.7,0.9,0.9);
  leg->AddEntry(necap,"Captured Neutrons","L");
  leg->AddEntry(nescat,"Scattered Neutrons","L");
  leg->AddEntry(neother,"Other Neutrons","L");
  leg->AddEntry(ngen,"Generated Neutrons","P");
  leg->Draw();
  ncan->cd(2);
  npcap->Draw("box");
  npscat->Draw("boxsame");
  npother->Draw("boxsame");
}
