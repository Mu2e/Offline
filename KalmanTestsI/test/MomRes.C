void MomRes(TFile *d, int fType=2, int nGenEv=-1, int nCuts=4, TString addCut="", int useExtapMom=0, TString fitOpt="L", TString fold="trackRecoFit", Long64_t NEntries=-1, Long64_t skipentries=0){
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("KalmanTestsI/test/KalFit.C");
  KalFitRes(t,nGenEv,nCuts,fitOpt,addCut,fType,useExtapMom,NEntries,skipentries);
}
