void AccPlots(TFile *d, int nGenEv=-1, TString addCut="", size_t cutType=0, double liveTime=-1, int useExtapMom=0, double pCorr=0.0, TString fold="trackRecoFit", Long64_t NEntries=-1, Long64_t skipentries=0){
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
  KalFitAcc(t,nGenEv,addCut,cutType,liveTime,useExtapMom,pCorr,NEntries,skipentries);
}
