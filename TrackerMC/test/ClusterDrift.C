#include "TTree.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPaveStats.h"

void ClusterDrift(TTree* nt) {
  TH2F* d2t = new TH2F("d2t","Cluster drift time vs drift distance;distance (mm);time (ns)",50,0.0,2.5,210,-1.0,50.5);
  nt->Project("d2t","1000*t:10*r");
  d2t->Draw();
  d2t->Draw("box");
  d2t->FitSlicesY();

  TCanvas* cdcan = new TCanvas("cdcan","Cluster Drift",1200,800);
  cdcan->Divide(2,2);
  cdcan->cd(1);
  d2t->SetStats(0);
  d2t->Draw("colorZ");
  cdcan->cd(2);
  TH1D* d2t_1 =(TH1D*)gDirectory->Get("d2t_1"); 
  TH1D* d2t_2 =(TH1D*)gDirectory->Get("d2t_2"); 
  d2t_1->SetTitle("Mean Cluster Drift time vs drift distance;distance (mm);mean time (ns)");
  d2t_1->Fit("pol2");
  TPaveStats* d2t_1_stats = (TPaveStats*)d2t_1->FindObject("stats");
  gPad->Update();
  if(d2t_1_stats != 0){
    d2t_1_stats->SetX1NDC(0.0);
    d2t_1_stats->SetX2NDC(0.4);
  }
  cdcan->cd(3);
  d2t_2->SetTitle("#sigma Cluster Drift time vs drift distance;distance (mm);#sigma time (ns)");
  d2t_2->Fit("pol2");
  TPaveStats * d2t_2_stats = (TPaveStats*)d2t_2->FindObject("stats");
  gPad->Update();
  if(d2t_2_stats != 0){
    d2t_2_stats->SetX1NDC(0.0);
    d2t_2_stats->SetX2NDC(0.4);
  }  
  cdcan->SaveAs("clusterdrift.png");
}

