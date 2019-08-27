#include "TMath.h"
void FitPBI(TFile* file){
  gStyle->SetOptFit(1111);
  TH2D* pbi = (TH2D*)file->Get("hln/rpbi");
  TF1* ln = new TF1("ln","[0]*TMath::LogNormal(x,[1],[2],[3])",0,4.0);
  ln->SetParameters(pbi->GetEntries(), 0.391, 0.0, 1.0);
  ln->SetParNames("Normalization","Sigma","Theta","Median");
  ln->FixParameter(2,0.0);
//  ln->FixParameter(3,1.0);
  pbi->Fit(ln);
  std::cout << "PBI Fit result lognormal_sigma = " << ln->GetParameter(1) << " lognormal_mu = " << TMath::Log(ln->GetParameter(3)) << std::endl;

}
