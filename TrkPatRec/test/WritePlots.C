void WritePlots(const char* filename) {
  TFile diagfile(filename);
  TTree* t = TrkPatRec->Get("trkdiag");
  gROOT->LoadMacro("TrkPatRec/test/FitTest.C");
  TFile resplots("resplots.root","RECREATE");
  MomRes(t);
  resplots.Write();
  resplots.Close();
}

