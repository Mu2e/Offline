#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
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
  TCut conv("mcproc==56&&mcpdg==11&&mcoe>100");
  TCut prot("mcpdg==2212");

  TH2F* eeve = new TH2F("eeve","Straw Hit Energy vs G4 Energy;G4 Energy (KeV);Reco Energy(KeV)",100,0,12.0,100,0,12.0);
  TH2F* peve = new TH2F("peve","Straw Hit Energy vs G4 Energy;G4 Energy (KeV);Reco Energy(KeV)",100,0,12.0,100,0,12.0);
  eeve->SetStats(0);
  peve->SetStats(0);
  eeve->SetMarkerStyle(6);
  peve->SetMarkerStyle(6);
  eeve->SetFillColor(kRed);
  peve->SetFillColor(kBlack);
  eeve->SetMarkerColor(kRed);
  peve->SetMarkerColor(kBlack);
  sh->Project("eeve","1000.0*edep:1000.0*mcedep",conv);
  sh->Project("peve","1000.0*edep:1000.0*mcedep",prot);

  TH1F* ee = new TH1F("ee","Energy Deposition;Straw energy (KeV)",100,0,10.0);
  TH1F* pe = new TH1F("pe","Energy Deposition;Straw energy (KeV)",100,0,10.0);
  ee->SetLineColor(kRed);
  pe->SetLineColor(kBlack);
  ee->SetStats(0);
  pe->SetStats(0);
//  ee->Sumw2();
  sh->Project("ee","1000.0*edep",conv);
  sh->Project("pe","1000.0*edep",prot);
  ee->Scale(escale);

  TH1F* ededx = new TH1F("ededx","dE/dx;Straw energy/pathlength (KeV/mm)",100,0,6);
  TH1F* pdedx = new TH1F("pdedx","dE/dx;Straw energy/pathlength (KeV/mm)",100,0,6);
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
  TLegend* eleg = new TLegend(0.6,0.7,0.9,0.9);
  eleg->AddEntry(ee,elab,"l");
  eleg->AddEntry(pe,plab,"l");
  snprintf(elab,40,"e^{-} eff=%3.3f",eff._deeff);
  snprintf(plab,40,"P^{+} rej=%3.3f",eff._dpeff);
  TLegend* deleg = new TLegend(0.6,0.7,0.9,0.9);
  deleg->AddEntry(ededx,elab,"l");
  deleg->AddEntry(pdedx,plab,"l");

  TCanvas* evecan = new TCanvas("evecan","E vs E", 600,600);
  peve->Draw();
  eeve->Draw("same");
  TLegend* eveleg = new TLegend(0.1,0.7,0.5,0.9);
  eveleg->AddEntry(eeve,"Ce e^{-}","p");
  eveleg->AddEntry(peve,"Protons","p");
  eveleg->Draw();

  TCanvas* ecan = new TCanvas("ecan","Energy",800,400);
  ecan->Divide(2,1);
  ecan->cd(1);
  pe->Draw();
  ee->Draw("same");
  eleg->Draw();
  TLine *lecut = new TLine(eff._ecut,0.0,eff._ecut,pe->GetMaximum());
  lecut->SetLineStyle(3);
  lecut->Draw();
  ecan->cd(2);
  pdedx->Draw();
  ededx->Draw("same");
  deleg->Draw();
  TLine *ldecut = new TLine(eff._decut,0.0,eff._decut,pdedx->GetMaximum());
  ldecut->SetLineStyle(3);
  ldecut->Draw();
  
  TCanvas* e2can = new TCanvas("e2can","e2can",700,700);
  pe->Draw();
  ee->Draw("same");
  eleg->Draw();
  lecut->SetLineStyle(3);
  lecut->Draw();

  
  if(fname !=0)
    e2can->SaveAs(fname);

}

void test(TTree* sh){
  EffStruct eff;
  eff._ecut = 3.0;
  eff._decut =  0.6;
  EPSeparation(sh,eff);
}

