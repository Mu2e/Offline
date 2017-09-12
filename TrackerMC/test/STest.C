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

void STest(TTree* sdiag, const char* page ="G4") {
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

    sstack->Add(slenh);
    sstack->Add(slenl);
 
    estack->Add(seh);
    estack->Add(sel);

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
    TH2F* nch = new TH2F("nch","N Clusters vs G4 Step Length;Step Length (mm)",50,0.0,10.0,60,-0.5,59.5);
    TH2F* ncl = new TH2F("ncl","N Clusters vs G4 Step Length;Step Length (mm)",50,0.0,10.0,20,-0.5,19.5);
    TH1F* nech = new TH1F("nech","Number of electrons/cluster",25,-0.5,24.5);
    TH1F* necl = new TH1F("necl","Number of electrons/cluster",50,-0.5,49.5);
    TH1F* slh = new TH1F("slh","Average Distance Between Clusters;step (mm)",100,0.0,2.0);
    TH1F* sll = new TH1F("sll","Average Distance Between Clusters;step (mm)",100,0.0,2.0);
    TH2F* eeh = new TH2F("eeh","Sum electron energy vs G4 Step Energy;G4 Step Energy (KeV);e Energy (KeV)",100,0,5.0,100,0.0,5.0);
    TH2F* eel = new TH2F("eel","Sum electron energy vs G4 Step Energy;G4 Step Energy (KeV);e Energy (KeV)",100,0,5.0,100,0.0,5.0);

    nch->SetStats(0);
//    eeh->SetStats(0);
    nech->SetFillColor(kRed);
    slh->SetFillColor(kRed);
    ncl->SetStats(0);
    necl->SetFillColor(kGreen);
    sll->SetFillColor(kGreen);
    TCut mini("partPDG==11&&partP>100");
    TCut highi("partPDG==2212||partPDG==11&&partP<1.0");
    TCut bigstep("steplen>1.0");

    sdiag->Project("nch","nclust:steplen",mini);
    sdiag->Project("ncl","nclust:steplen",highi);
    sdiag->Project("nech","clusters._ne",mini);
    sdiag->Project("necl","clusters._ne",highi);
    sdiag->Project("slh","steplen/nclust",mini&&bigstep);
    sdiag->Project("sll","steplen/nclust",highi);
    sdiag->Project("eeh","eesum*1.0e3:stepE*1.0e3",mini);
    sdiag->Project("eel","eesum*1.0e3:stepE*1.0e3",highi);

    TCanvas* hcan = new TCanvas("hcan","hcan",600,600);
    hcan->Divide(2,2);
    hcan->cd(1);
    nch->Draw("colorz");
    hcan->cd(2);
    slh->Draw();
    hcan->cd(3);
    nech->Draw();
    hcan->cd(4);
    eeh->Draw("colorz");

    TCanvas* lcan = new TCanvas("lcan","lcan",600,600);
    lcan->Divide(2,2);
    lcan->cd(1);
    ncl->Draw("colorz");
    lcan->cd(2);
    sll->Draw();
    lcan->cd(3);
    necl->Draw();
    lcan->cd(4);
    eel->Draw("colorz");
  }
}
