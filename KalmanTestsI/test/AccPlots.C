void AccPlots(TFile *d, int nGenEv=-1, TString addCut="", size_t cutType=0, double liveTime=-1, int useExtapMom=0, double pCorr=0.0, TString fold="kalmanFit", Long64_t NEntries=-1, Long64_t skipentries=0){
  //TFile d(fName.Data());
  //TTree* t = ReadKalFits->Get("trkdiag");
  if (nGenEv<=0) nGenEv=((TH1F*)d->Get("g4run/totalMultiplicity"))->GetEntries();
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
  KalFitAcc(t,nGenEv,addCut,cutType,liveTime,useExtapMom,pCorr,NEntries,skipentries);
}
