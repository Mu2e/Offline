#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TTree.h"
#include "THStack.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TRandom3.h"
#include <iostream>
#include <math.h>
#include <vector>

void StereoRecoTest(TTree* sdiag) {
  TCut mcrel("mcrel==0");
  TCut unrel("mcrel<0");
  TH1F* chi2r = new TH1F("chi2r","Transverse Matching Chisquared",200,-1.0,101.0);
  TH1F* chi2u = new TH1F("chi2u","Transverse Matching Chisquared",200,-1.0,101.0);
  chi2r->SetLineColor(kBlue);
  chi2u->SetLineColor(kRed);
  chi2r->SetStats(0);
  chi2u->SetStats(0);

  TH1F* timer = new TH1F("timer","Hit Time Difference;#Delta t (ns)",200,-1.0,41.0);
  TH1F* timeu = new TH1F("timeu","Hit Time Difference;#Delta t (ns)",200,-1.0,41.0);
  timer->SetLineColor(kBlue);
  timeu->SetLineColor(kRed);
  timer->SetStats(0);
  timeu->SetStats(0);

  TH1F* der = new TH1F("der","Fractional Energy Difference;#Delta E/#Sigma E",200,-0.01,1.01);
  TH1F* deu = new TH1F("deu","Fractional Energy Difference;#Delta E/#Sigma E",200,-0.01,1.01);
  der->SetLineColor(kBlue);
  deu->SetLineColor(kRed);
  der->SetStats(0);
  deu->SetStats(0);

  TH1F* mvar = new TH1F("mvar","MVA Output",200,-0.1,1.5);
  TH1F* mvau = new TH1F("mvau","MVA Output",200,-0.1,1.5);
  mvar->SetLineColor(kBlue);
  mvau->SetLineColor(kRed);
  mvar->SetStats(0);
  mvau->SetStats(0);

  sdiag->Project("chi2r","chi2",mcrel);
  sdiag->Project("chi2u","chi2",unrel);
  sdiag->Project("timer","dt",mcrel);
  sdiag->Project("timeu","dt",unrel);
  sdiag->Project("der","de",mcrel);
  sdiag->Project("deu","de",unrel);
  sdiag->Project("mvar","mvaout",mcrel);
  sdiag->Project("mvau","mvaout",unrel);

  TLegend* sdleg = new TLegend(0.6,0.7,0.9,0.9);
  sdleg->AddEntry(chi2r,"MC match","L");
  sdleg->AddEntry(chi2u,"MC mis-match","L");


  TCanvas* sdcan = new TCanvas("sdcan","Stereo Matching",1200,800);
  sdcan->Divide(2,2);
  sdcan->cd(1);
  gPad->SetLogy();
  chi2r->Draw();
  chi2u->Draw("same");
  sdleg->Draw();
  sdcan->cd(2);
  timer->Draw();
  timeu->Draw("same");
  sdcan->cd(3);
  der->Draw();
  deu->Draw("same");
  sdcan->cd(4);
  gPad->SetLogy();
  mvar->Draw();
  mvau->Draw("same");
}
