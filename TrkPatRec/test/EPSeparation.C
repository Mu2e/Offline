#include "TTree.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCut.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>

struct EffStruct {
  double _ecut;
  double _decut;
  double _eeff, _eefferr;
  double _deeff, _deefferr;
  double _peff, _pefferr;
  double _dpeff, _dpefferr;
};

void EPSeparation(TTree* sh, EffStruct& eff,const char* fname=0) {
  double escale(1.0);
  TCut conv("mcgen==2");
  TCut prot("mcgen==28&&mcpdg==2212");
  TH1F* ee = new TH1F("ee","Energy Deposition;Straw energy (KeV)",100,0,20.0);
  TH1F* pe = new TH1F("pe","Energy Deposition;Straw energy (KeV)",100,0,20.0);
  ee->SetLineColor(kRed);
  pe->SetLineColor(kBlack);
  ee->SetStats(0);
  pe->SetStats(0);
//  ee->Sumw2();
  sh->Project("ee","1000.0*edep",conv);
  sh->Project("pe","1000.0*edep",prot);
  ee->Scale(escale);

  TH1F* edist = new TH1F("edist","MC true DOCA; Distance to wire (mm)",100,0,2.51);
  TH1F* pdist = new TH1F("pdist","MC true DOCA; Distance to wire (mm)",100,0,2.51);
  TH1F* epath = new TH1F("epath","MC true gas path; gas path (mm)",100,0,5.01);
  TH1F* ppath = new TH1F("ppath","MC true gas path; gas path (mm)",100,0,5.01);
  edist->SetLineColor(kRed);
  pdist->SetLineColor(kBlack);
  epath->SetLineColor(kRed);
  ppath->SetLineColor(kBlack);
  epath->SetStats(0);
  ppath->SetStats(0);
//  edist->Sumw2();
//  epath->Sumw2();
  sh->Project("edist","mcshd",conv);
  sh->Project("pdist","mcshd",prot);
  sh->Project("epath","2*sqrt(2.5^2-mcshd^2)",conv);
  sh->Project("ppath","2*sqrt(2.5^2-mcshd^2)",prot);
  edist->Scale(escale);
  epath->Scale(escale);
  edist->SetStats(0);
  pdist->SetStats(0);

  TProfile *ed = new TProfile("ed","Average energy vs pathlength",100,0,5.01);
  TProfile *pd = new TProfile("pd","Average energy vs pathlength",100,0,5.011);
  ed->SetLineColor(kRed);
  pd->SetLineColor(kBlack);
  sh->Project("ed","edep:2*sqrt(2.5^2-mcshd^2)",conv);
  sh->Project("pd","edep:2*sqrt(2.5^2-mcshd^2)",prot);
  ed->Scale(escale);
  ed->SetStats(0);
  pd->SetStats(0);

  TH1F* ededx = new TH1F("ededx","dE/dx;Straw energy/pathlength (KeV/mm)",100,0,5);
  TH1F* pdedx = new TH1F("pdedx","dE/dx;Straw energy/pathlength (KeV/mm)",100,0,5);
  ededx->SetLineColor(kRed);
  pdedx->SetLineColor(kBlack);
//  ededx->Sumw2();
  sh->Project("ededx","500.*edep/sqrt(2.5^2-mcshd^2)",conv);
  sh->Project("pdedx","500.*edep/sqrt(2.5^2-mcshd^2)",prot);
  ededx->Scale(escale);
  ededx->SetStats(0);
  pdedx->SetStats(0);

  int ibin = ee->FindFixBin(eff._ecut);
  double eint = ee->Integral(1,ibin);
  double etot = ee->Integral();
  eff._eeff = eint/etot; 
  eff._eefferr = sqrt(eint*(etot-eint)/pow(etot,3));
  double pint = pe->Integral(1,ibin);
  double ptot = pe->Integral();
  eff._peff = 1.0-pint/ptot;
  eff._pefferr = sqrt(pint*(ptot-pint)/pow(ptot,3));
  std::cout << "Electron efficiency = " << eff._eeff << " +- " << eff._eefferr
    << " Proton rejection = " << eff._peff << " +- " << eff._pefferr << std::endl;

  int jbin = ededx->FindFixBin(eff._decut);
  double deint = ededx->Integral(1,jbin);
  double detot = ededx->Integral();
  eff._deeff = deint/detot;
  eff._deefferr = sqrt(deint*(detot-deint)/pow(detot,3));
  double dpint = pdedx->Integral(1,jbin);
  double dptot = pdedx->Integral();
  eff._dpeff = 1.0-dpint/dptot;
  eff._dpefferr = sqrt(dpint*(dptot-dpint)/pow(dptot,3));
  std::cout << "Electron efficiency = " << eff._deeff << " +- " << eff._deefferr
    << " Proton rejection = " << eff._dpeff << " +- " << eff._dpefferr << std::endl;

  char elab[40];
  char plab[40];
  snprintf(elab,40,"e^{-} eff=%3.3f",eff._eeff);
  snprintf(plab,40,"P^{+} rej=%3.3f",eff._peff);
  TLegend* eleg = new TLegend(0.1,0.7,0.5,0.9);
  eleg->AddEntry(ee,elab,"l");
  eleg->AddEntry(pe,plab,"l");
  snprintf(elab,40,"e^{-} eff=%3.3f",eff._deeff);
  snprintf(plab,40,"P^{+} rej=%3.3f",eff._dpeff);
  TLegend* deleg = new TLegend(0.1,0.7,0.5,0.9);
  deleg->AddEntry(ededx,elab,"l");
  deleg->AddEntry(pdedx,plab,"l");

  TCanvas* ecan = new TCanvas("ecan","Energy",1200,800);
  ecan->Divide(3,2);
  ecan->cd(1);
  pe->Draw();
  ee->Draw("same");
  eleg->Draw();
  TLine *lecut = new TLine(eff._ecut,0.0,eff._ecut,pe->GetMaximum());
  lecut->SetLineStyle(3);
  lecut->Draw();
  ecan->cd(2);
  pdist->Draw();
  edist->Draw("same");
  ecan->cd(3);
  ppath->Draw();
  epath->Draw("same");
  ecan->cd(4);
  pd->Draw();
  ed->Draw("same");
  ecan->cd(5);
  pdedx->Draw();
  ededx->Draw("same");
  deleg->Draw();
  TLine *ldecut = new TLine(eff._decut,0.0,eff._decut,pdedx->GetMaximum());
  ldecut->SetLineStyle(3);
  ldecut->Draw();
  if(fname !=0)
    ecan->SaveAs(fname);
//  delete ecan;
}
