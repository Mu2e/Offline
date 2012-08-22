void MomRes(TFile *d, int nGenEv=-1, TString fitOpt="L", TString fold="kalmanFit"){
  //TFile d(fName.Data());
  //TTree* t = ReadKalFits->Get("trkdiag");
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("Users/ignatov/KalmanFitTests/test/KalFit.C");
  KalFitRes(t,nGenEv,fitOpt);
}
