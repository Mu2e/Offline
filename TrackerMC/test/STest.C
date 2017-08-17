#include "THStack.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

void STest(TTree* sdiag, TTree* cdiag, const char* page ="G4") {
  TString spage(page);
  if(spage == "G4") {
    THStack* sstack = new THStack("sstack","G4 Step Length in Straw Gas;step (mm)");
    THStack* estack = new THStack("estack","G4 Ionization Energy Left in Straw Gas;E_{Ion} (KeV)");

//    TH1F* slen = new TH1F("slen","G4 Step Length in Straw Gas;step (mm)",100,0,10.0);
    TH1F* slenh = new TH1F("slenh","G4 Step Length in Straw Gas;step (mm)",100,0,10.0);
    TH1F* slenl = new TH1F("slenl","G4 Step Length in Straw Gas;step (mm)",100,0,10.0);
//    TH1F* se = new TH1F("se","G4 Ionization Energy Left in Straw Gas;E_{Ion} (KeV)",100,0,5.0);
    TH1F* seh = new TH1F("seh","G4 Ionization Energy Left in Straw Gas;E_{Ion} (KeV)",100,0,5.0);
    TH1F* sel = new TH1F("sel","G4 Ionization Energy Left in Straw Gas;E_{Ion} (KeV)",100,0,5.0);
    TH2F* selen = new TH2F("selen","G4 Ionization Energy vs Step Length;step (mm);E_{Ion} (KeV)",50,0.5,10.0,50,0,5.0);
    TH1F* pe = new TH1F("pe","Particle Energy at Step;Energy (MeV)",100,0,106);

    selen->SetStats(0);
//    slen->SetStats(0);
    slenh->SetStats(0);
    slenl->SetStats(0);
    slenh->SetLineColor(kRed);
    slenh->SetFillColor(kRed);
    slenl->SetLineColor(kGreen);
    slenl->SetFillColor(kGreen);
//    se->SetStats(0);
    seh->SetStats(0);
    sel->SetStats(0);
    seh->SetLineColor(kRed);
    seh->SetFillColor(kRed);
    sel->SetLineColor(kGreen);
    sel->SetFillColor(kGreen);
    pe->SetStats(0);
//    sdiag->Project("slen","steplen");
    sdiag->Project("slenh","steplen","partP>100");
    sdiag->Project("slenl","steplen","partP<5");
//    sdiag->Project("se","stepE*1000");
    sdiag->Project("seh","stepE*1000","partP>100");
    sdiag->Project("sel","stepE*1000","partP<5");
    sdiag->Project("selen","stepE*1000:steplen");
    sdiag->Project("pe","partP");

    sstack->Add(slenl);
    sstack->Add(slenh);
 
    estack->Add(sel);
    estack->Add(seh);

    TCanvas* g4can = new TCanvas("g4can","g4can",1000,1000);
    g4can->Divide(2,2);
    g4can->cd(1);
    sstack->Draw();
//    slen->Draw();
//    slenh->Draw("same");
//    slenl->Draw("same");
    g4can->cd(2);
    selen->Draw("colorz");
    g4can->cd(3);
//    se->Draw();
//    seh->Draw("same");
//    sel->Draw("same");
    estack->Draw();
    TLegend* eleg = new TLegend(0.5,0.6,0.9,0.9);
    eleg->AddEntry(seh,"E_{e} > 100 MeV");
    eleg->AddEntry(sel,"E_{e} < 5 MeV");
//    eleg->AddEntry(se,"All E_{e}");
    eleg->Draw();
    g4can->cd(4);
    pe->Draw();

  } else if(spage == "cluster") {
    THStack* cstack = new THStack("cstack","Number of clusters");
    THStack* slstack = new THStack("slstack","Step Length Between Ionizations;step (mm)");
//    TH1F* nc = new TH1F("nc","Number of clusters",50,-0.5,49.5);
    TH1F* nch = new TH1F("nch","Number of clusters",50,-0.5,49.5);
    TH1F* ncl = new TH1F("ncl","Number of clusters",50,-0.5,49.5);
    TH2F* ncs = new TH2F("ncs","Number of clusters vs G4 step length;step (mm)",50,0.5,10.0,50,-0.5,49.5);
    TH1F* ne = new TH1F("ne","Number of electrons/cluster",15,-0.5,14.5);
    TH1F* ee = new TH1F("ee","Energy per electrons;ev",100,0.0,150.0);
    TH1F* gp = new TH1F("gp","Straw Gain",100,0,4e5);
    TH1F* slh = new TH1F("slh","Step Length Between Ionizations;step (mm)",100,0.0,2.0);
    TH1F* sll = new TH1F("sll","Step Length Between Ionizations;step (mm)",100,0.0,2.0);

//    nc->SetStats(0);
    nch->SetStats(0);
//    ncl->SetStats(0);
    slh->SetStats(0);
//    sll->SetStats(0);
    nch->SetFillColor(kRed);
    slh->SetFillColor(kRed);
    ncl->SetFillColor(kGreen);
    sll->SetFillColor(kGreen);
    ncs->SetStats(0);
    ne->SetStats(1);
    ee->SetStats(1);

//    sdiag->Project("nc","nsubstep");
    sdiag->Project("nch","nsubstep","partP>100");
    sdiag->Project("ncl","nsubstep","partP<5");
    sdiag->Project("ncs","nsubstep:steplen");
    cdiag->Project("ne","nion");
    sdiag->Project("ee","1.0e6*stepE/niontot");
    cdiag->Project("gp","gain");
    sdiag->Project("slh","steplen/nsubstep","partP>100");
    sdiag->Project("sll","steplen/nsubstep","partP<5");
    

    cstack->Add(ncl);
    cstack->Add(nch);
    slstack->Add(sll);
    slstack->Add(slh);

    TCanvas* ccan = new TCanvas("ccan","ccan",1200,800);
    ccan->Divide(3,2);
    ccan->cd(1);
    ne->Draw();
    ccan->cd(2);
//    nc->Draw();
//    nch->Draw("same");
//    ncl->Draw("same");
    cstack->Draw();
    TLegend* cleg = new TLegend(0.5,0.6,0.9,0.9);
    cleg->AddEntry(nch,"E_{e} > 100 MeV");
    cleg->AddEntry(ncl,"E_{e} < 5 MeV");
//    cleg->AddEntry(nc,"All E_{e}");
    cleg->Draw();
    ccan->cd(3);
    ncs->Draw("colorz");
    ccan->cd(4);
    ee->Draw();
    ccan->cd(5);
    gp->Draw();
    ccan->cd(6);
    slstack->Draw();

  } else if(spage == "gain") {

  }
}
