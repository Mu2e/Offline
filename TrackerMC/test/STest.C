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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraph.h"

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
    TProfile* nchp = new TProfile("nchp","N Clusters vs G4 Step Length;Step Length (mm)",50,0.0,10.0,-0.5,99.5,"S");
    TH1F* ncl = new TH1F("ncl","N Clusters",12,-0.5,11.5);
    TH1F* nech = new TH1F("nech","Number of electrons/cluster",25,-0.5,24.5);
    TH1F* necl = new TH1F("necl","Number of electrons/cluster",300,-0.5,299.5);
    TH2F* eeh = new TH2F("eeh","Sum electron energy vs G4 Step Energy;G4 Step Energy (KeV);#Sigma N_{e}#timesE_{avg} (KeV)",100,0,5.0,100,0.0,5.0);
    TH2F* eel = new TH2F("eel","Sum electron energy vs G4 Step Energy;G4 Step Energy (KeV);#Sigma N_{e}#timesE_{avg} Energy (KeV)",100,0,100.0,100,0.0,100.0);

    nch->SetStats(0);
    nchp->SetStats(0);
    necl->SetStats(0);
    eel->SetStats(0);
//    eeh->SetStats(0);
    nech->SetFillColorAlpha(kRed,0.5);
    ncl->SetStats(0);
    eeh->SetStats(0);
    necl->SetFillColor(kGreen);
    TCut mini("partPDG==11&&partP>100");
    TCut highi("partPDG==2212||partPDG==11&&partP<0.3");

    sdiag->Project("nch","nclust:steplen",mini);
    sdiag->Project("nchp","nclust:steplen",mini);
    sdiag->Project("ncl","nclust",highi);
    sdiag->Project("nech","clusters._ne",mini);
    sdiag->Project("necl","clusters._ne",highi);
    sdiag->Project("eeh","eesum*1.0e3:stepE*1.0e3",mini);
    sdiag->Project("eel","eesum*1.0e3:stepE*1.0e3",highi);

    // create a histogram from the prob. distribution and scale it
    std::vector<double> nProb{0.656,0.15,0.064,0.035,0.0225,0.0155,0.0105,
      0.0081,0.0061, 0.0049, 0.0039, 0.0030, 0.0025, 0.0020, 0.0016, 0.0012, 0.00095, 0.00075}; // Blum, table 1.4
    TH1F* neB = new TH1F("neB","Number of electrons/cluster",25,-0.5,24.5);
    for(size_t iprob=0;iprob<nProb.size();++iprob){
      neB->Fill(iprob+1.0,nProb[iprob]);
    }
    neB->SetStats(0);
    neB->Scale(nech->GetEntries());
    neB->SetLineColor(kBlack);
//    neB->SetFillColor(kYellow);
    TFitResultPtr fp = nchp->Fit("pol1","+S");
    double slope = fp->Parameter(1);

    std::vector<float> nmean;
    std::vector<float> nspread;
    for(int ibin=2;ibin<nchp->GetNbinsX();++ibin){
      nmean.push_back(sqrt(nchp->GetBinContent(ibin)));
      nspread.push_back(nchp->GetBinError(ibin));
    }
    TGraph* ncs = new TGraph(nmean.size(),nmean.data(),nspread.data());
    ncs->SetTitle("N_{c} Spread vs sqrt(<N_{c}>);sqrt(<N_{c}>);#sigma_{Nc}");

    char legtit[80];
    TCanvas* hcan = new TCanvas("hcan","hcan",600,600);
    hcan->Divide(2,2);
    hcan->cd(1);
    nch->Draw("colorz");
    nchp->Draw("same");
    TLegend* fleg = new TLegend(0.2,0.7,0.6,0.9);
    snprintf(legtit,80,"Fit Slope = %3.2f mm^{-1}",slope);
    fleg->AddEntry(nchp,legtit,"");
    fleg->Draw();
    hcan->cd(2);
    ncs->Draw("ACP");
    hcan->cd(3);
    neB->Draw("HIST");
    nech->Draw("same");
    neB->Draw("HISTsame");
    TLegend* leg = new TLegend(0.4,0.7,0.9,0.9);
    snprintf(legtit,80,"Input <N_{e}> = %3.2f",neB->GetMean());
    leg->AddEntry(neB,legtit,"LF");
    snprintf(legtit,80,"Sim Result <N_{e}> = %3.2f",nech->GetMean());
    leg->AddEntry(nech,legtit,"F");
    leg->Draw();
    hcan->cd(4);
    eeh->Draw("colorz");
    TLine* diag = new TLine(0.0,0.0,5.0,5.0);
    diag->SetLineColor(kRed);
    diag->Draw();
    
    TCanvas* lcan = new TCanvas("lcan","lcan",600,600);
    lcan->Divide(2,2);
    lcan->cd(1);
    ncl->Draw();
    lcan->cd(2);
    necl->Draw();
    lcan->cd(3);
    gPad->SetLogz();
    eel->Draw("colorz");
    TLine* diag2 = new TLine(0.0,0.0,100.0,100.0);
    diag2->SetLineColor(kRed);
    diag2->Draw();
  }
}
