void MomPull(TFile *d, int fType=1, int nGenEv=-1, int nCuts=4, TString addCut="", bool useExtapMom=false, TString fitOpt="L", TString fold="kalmanFit"){
  //TFile d(fName.Data());
  //TTree* t = ReadKalFits->Get("trkdiag");
  if (nGenEv<=0) nGenEv=((TH1F*)d->Get("g4run/totalMultiplicity"))->GetEntries();
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
  KalFitMomPull(t,nGenEv,nCuts,fitOpt,addCut,fType,useExtapMom);
}
