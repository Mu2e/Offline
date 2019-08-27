#include "TMath.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

void FitPBI(TFile* file){
  gStyle->SetOptFit(1111);
  TH1D* pbi = (TH1D*)file->Get("hln/rpbi");
  TF1* ln = new TF1("ln","[0]*TMath::LogNormal(x,[1],[2],[3])",0,4.0);
  ln->SetParameters(pbi->GetEntries(), 0.391, 0.0, 1.0);
  ln->SetParNames("Normalization","Sigma","Theta","Median");
  ln->FixParameter(2,0.0);
//  ln->FixParameter(3,1.0);
  pbi->Fit(ln);
  std::cout << "PBI Fit result lognormal_sigma = " << ln->GetParameter(1) << " lognormal_mu = " << TMath::Log(ln->GetParameter(3)) << std::endl;
  TTree ihep("ihep","ihep");
  // data is at bin lower edges, not all bins were measured
  ihep.ReadFile("EventGenerator/test/IHEP_PBI.dat","Irel/F:freq/F");
  TH1F* RInten = new TH1F("RInten","Relative Beam Intensity",46,0.0,3.0);
  RInten->Sumw2();
  ihep.Project("RInten","Irel+0.03243","freq");
  for(unsigned ibin=1;ibin <= RInten->GetNbinsX();++ibin){
    RInten->SetBinError(ibin,0.02);
  }
  RInten->SetMarkerStyle(20);
  RInten->SetMarkerColor(kGreen);
  RInten->Scale(pbi->Integral()*RInten->GetXaxis()->GetNbins()/(RInten->Integral()*pbi->GetXaxis()->GetNbins()));
  RInten->Draw("same");
  TLegend* leg = new TLegend(0.4,0.7,0.6,0.9);
  leg->AddEntry(RInten,"IHEP Data","P");
  leg->AddEntry(pbi,"MDC2018 Generator","L");
  leg->Draw();
}
