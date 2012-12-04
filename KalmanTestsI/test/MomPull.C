void MomPull(TFile *d, int fType=1, int nGenEv=-1, int nCuts=4, TString addCut="", bool useExtapMom=false, TString fitOpt="L", TString fold="trackRecoFit", Long64_t NEntries=-1, Long64_t skipentries=0){
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
  KalFitMomPull(t,nGenEv,nCuts,fitOpt,addCut,fType,useExtapMom,NEntries,skipentries);
}
