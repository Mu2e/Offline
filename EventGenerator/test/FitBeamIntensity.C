#include "TMath.h"
{
  gStyle->SetOptFit(1111);
// bin size = 0.06486
  TTree ihep("ihep","ihep");
  // data is at bin lower edges, not all bins were measured
  ihep.ReadFile("EventGenerator/test/IHEP_PBI.dat","Irel/F:freq/F");
  RInten = new TH1F("RInten","Relative Beam Intensity",46,0.0,3.0);
  RInten->Sumw2();
  ihep.Project("RInten","Irel+0.03243","freq");
  for(unsigned ibin=1;ibin <= RInten->GetNbinsX();++ibin){
    RInten->SetBinError(ibin,0.02);
  }
  ln = new TF1("ln","[0]*TMath::LogNormal(x,[1],[2],[3])",0,4.0);
  ln->SetParameters(1.0, 0.391, 0.0, 1.0);
  ln->SetParNames("Normalization","Sigma","Theta","Median");
  ln->FixParameter(2,0.0);
//  ln->FixParameter(3,1.0);
  RInten->Fit(ln);
  std::cout << "PBI Fit result lognormal_sigma = " << ln->GetParameter(1) << " lognormal_mu = " << TMath::Log(ln->GetParameter(3)) << std::endl;

}
