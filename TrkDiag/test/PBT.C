#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"

void PBTO(TTree* tree) {
  TProfile* meanprof = new TProfile("meanprof","Proton Bunch Mean Time vs EWM Offset;EWM offset (ns);PBTO (ns)",50,-12.5,12.5,-50,50);
  TH2F* mean2d = new TH2F("mean2d","Proton Bunch Mean Time vs EWM Offset;EWM offset (ns);PBTO (ns)",50,-12.5,12.5,30,-50,50);
  TH1F* res = new TH1F("res","Proton Bunch Mean Time Resolution;#Delta t (ns)",50,-50,50);
  TH1F* pull = new TH1F("pull","Proton Bunch Mean Time Pull;#Delta t/#sigma_t",50,-10,10);
  TH1F* terr = new TH1F("terr","Estimated Error on Proton Bunch Mean Time;$sigma_t (ns)",50,0,30);
  tree->Project("mean2d","tmean:-ewmoffset");
  tree->Project("meanprof","tmean:-ewmoffset");
  tree->Project("res","tmean+ewmoffset");
  tree->Project("pull","(tmean+ewmoffset)/terr");
  tree->Project("terr","terr");
  mean2d->SetStats(0);

  TCanvas* can = new TCanvas("can","can",1600,1200);
  can->Divide(2,2);
  can->cd(1);
  mean2d->Draw("colorz");
  meanprof->Fit("pol1","","sames");
  can->cd(2);
  res->Fit("gaus","","",-20,20);
  can->cd(3);
  pull->Fit("gaus");
  can->cd(4);
  terr->Draw();
  //  median2d->Draw();
  //  medianprof->Fit("pol1","","sames");

}

