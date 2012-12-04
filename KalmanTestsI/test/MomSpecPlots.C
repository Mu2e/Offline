void MomSpecPlots(TFile *d, TString addCut="", size_t cutType=0, double scale=1.0, double liveTime=-1, int useExtapMom=0, double pCorr=0.0, Long64_t NEntries=-1, Long64_t skipentries=0, TString fold="kalmanFit"){
  //TFile d(fName.Data());
  //TTree* t = ReadKalFits->Get("trkdiag");
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
  KalPlotMomSpec(t,addCut,cutType,scale,liveTime,useExtapMom,pCorr,NEntries,skipentries);
}
