#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "Math.h"
#include "Math/ChebyshevPol.h"
// Chebyshev polynomial expansion up to order 5
// input goes from -1 to 1, range from -1 to 1
Double_t chebyshev5(Double_t *x, Double_t *par) {
  double arg = par[7]*(x[0] + par[6]);
  return ROOT::Math::ChebyshevN(5,arg,par);
}

void DrawCheby() {
  TF1* cheby5 = new TF1("cheby5",chebyshev5,0,50,8);
  cheby5->SetParameters(1.25,2.5,0.1,0.0,0.0,0.0,-22,0.022);
  cheby5->Draw();
}

void FitT2D(TTree* trks) {
  TF1* cheby5 = new TF1("cheby5",chebyshev5,0,50,8);
  cheby5->SetParameters(1.25,2.5,0.1,0.0,0.0,0.0,-22,0.022);

  TH2F* rt2d = new TH2F("rt2d","Reco drift distance vs #Delta t;Hit t - MC t_{0} (nsec);Reco Drift Distance (mm)",100,0,50.0,100,-0.05,2.55);
  TH2F* tt2d = new TH2F("tt2d","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",25,0,50.0,50,-0.05,2.55);
  TProfile* rt2dp = new TProfile("rt2dp","Reco drift distance vs #Delta t;Hit t - MC t_{0} (nsec);Reco Drift Distance (mm)",100,0,50.0,-0.05,2.55);
  TProfile* tt2dp = new TProfile("tt2dp","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",100,0,50.0,-0.05,2.55);
  TProfile* rt2dp2 = new TProfile("rt2dp2","Reco drift distance vs #Delta t;Hit t - MC t_{0} (nsec);Reco Drift Distance (mm)",100,0,50.0,-0.05,2.55,"s");
  TProfile* tt2dp2 = new TProfile("tt2dp2","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",100,0,50.0,-0.05,2.55,"s");
  trks->Project("tt2d","_mcdist:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");
  trks->Project("rt2d","_rdrift:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");
  trks->Project("tt2dp","_mcdist:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");
  trks->Project("rt2dp","_rdrift:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");
  trks->Project("tt2dp2","_mcdist:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");
  trks->Project("rt2dp2","_rdrift:_ht-_mct0","fitstatus>0&&_active&&fitcon>1e-2");

  TCanvas* t2dcan = new TCanvas("t2dcan","t2dcan",1200,800);
  t2dcan->Divide(2,2);
  t2dcan->cd(1);
  rt2d->Draw(); 
  t2dcan->cd(2);
  tt2d->Draw(); 
  t2dcan->cd(3);
  rt2dp->Fit("pol1","","",0,43);
  t2dcan->cd(4);
  tt2dp->Fit(cheby5,"","",0,43);

  TH1F* dresid = new TH1F("dresid","Drift residual RMS vs #Delta t;Hit t - MC t_{0};Drift RMS",100,0,50.0);
  for(unsigned ibin=0;ibin<100;++ibin){
    dresid->Fill(tt2dp2->GetBinCenter(ibin+1),tt2dp2->GetBinError(ibin+1));
  }
  tt2d->FitSlicesY();
  TH1D *tt2d_1 = (TH1D*)gDirectory->Get("tt2d_1");
  TH1D *tt2d_2 = (TH1D*)gDirectory->Get("tt2d_2");

  TCanvas* t2drcan = new TCanvas("t2drcan","t2drcan",800,800);
  t2drcan->Divide(2,1);
  t2drcan->cd(1);
  dresid->Draw();
  t2drcan->cd(2);
  tt2d_2->Draw();
}

