
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"

void PBTO(TTree* tree) {
  TProfile* meanprof = new TProfile("meanprof","Early Straw Hit Time vs EWM Offset;EWM offset (ns);PBTO (ns)",30,-12.5,12.5,-50,200);
  TH2F* mean2d = new TH2F("mean2d","Early Straw Hit Time vs EWM Offset;EWM offset (ns);PBTO (ns)",30,-12.5,12.5,30,-50,200);
  TH1F* ewmres = new TH1F("ewmres","Early Straw Hit Time Resolution;#Delta t (ns)",30,-50,200);
  TH1F* mcres = new TH1F("mcres","Early Straw Hit Time Resolution;#Delta t (ns)",30,-50,50);
  TH1F* mcewmres = new TH1F("mcewmres","MC Hit Time Resolution;#Delta tt (ns)",30,-50,200);
  TCut earlyhit("correcttime-ewmoffset<200");
  tree->Project("mean2d","correcttime-ewmoffset:-ewmoffset",earlyhit);
  tree->Project("meanprof","correcttime-ewmoffset:-ewmoffset",earlyhit);
  tree->Project("ewmres","correcttime",earlyhit);
  tree->Project("mcres","correcttime-ewmoffset-mcsptime",earlyhit);
  tree->Project("mcewmres","mcsptime",earlyhit);
  mean2d->SetStats(0);

  TCanvas* can = new TCanvas("can","can",1600,1200);
  can->Divide(2,2);
  can->cd(1);
  mean2d->Draw("colorz");
  meanprof->Fit("pol1","","sames");
  can->cd(2);
  ewmres->Fit("gaus");
  can->cd(3);
  mcres->Fit("gaus");
  can->cd(4);
  mcewmres->Fit("gaus");
  //  median2d->Draw();
  //  medianprof->Fit("pol1","","sames");

}

