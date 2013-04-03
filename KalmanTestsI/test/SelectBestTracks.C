    struct {
      int fit;
      float chi2,chi2out;
      int nhits,nhitstot,nturn;
      float momin,momout;
      float fitmom, fitmomerr;
      float t0,t0fit,errt0;
      int nHitSwire, nhitFwire;
      float pathInSwire, pathInFwire, pathInGas, simPathInGas, simPathInGasPerLoop;
      void clear(){
	fit=0;chi2=-1;chi2out=-1;nhits=0;nhitstot=0;nturn=0;momin=0;momout=0;
	fitmom=0; fitmomerr=0;t0=0;t0fit=0;errt0=-1; nHitSwire=0; nhitFwire=0;
	pathInSwire=0.0; pathInFwire=0.0; pathInGas=0.0; simPathInGas=0.0; simPathInGasPerLoop=0.0;};
    } recoinfo;


void SelectBestTracks(TFile *d, int nGenEv=-1, TString addCut="", double liveTime=-1, TString fold="kalmanFit"){

  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());

  TFile outf("out_tmp.root","recreate");

  TTree *treefitNew= new  TTree("trfit","track fit params");
  treefitNew->Branch("fitinfo",&recoinfo,"fit/I:chi2/F:chi2out/F:nhits/I:nhitstot/I:nturn/I:momin/F:momout/F:fitmom/F:fitmomerr/F:t0/F:t0fit/F:errt0/F:nHitSwire/I:nhitFwire/I:pathInSwire/F:pathInFwire/F:pathInGas/F:simPathInGas/F:simPathInGasPerLoop/F");

}
