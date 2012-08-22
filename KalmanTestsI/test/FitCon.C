void FitCon(TFile *d){
  //TFile d(fName.Data());
  //TTree* t = ReadKalFits->Get("trkdiag");
  TTree* t = (TTree *) d->Get("kalmanFit/trfit");
  gROOT->LoadMacro("Users/ignatov/KalmanFitTests/test/KalFit.C");
  KalFitCon(t);
}
