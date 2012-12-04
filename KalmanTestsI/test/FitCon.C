void FitCon(TFile *d, TString fold="trackRecoFit", Long64_t NEntries=-1, Long64_t skipentries=0){
        TString treeName = fold+"/trfit";
        TTree* t = (TTree *) d->Get(treeName.Data());
        gROOT->LoadMacro("$MU2E_BASE_RELEASE/KalmanTestsI/test/KalFit.C");
        KalFitCon(t,0,NEntries,skipentries);
}
