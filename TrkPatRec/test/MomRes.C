{
  char* filename = gSystem->Getenv("FILE");
  TFile d(filename);
  TTree* t = TrkPatRec.Get("trkdiag");
  gROOT->LoadMacro("TrkPatRec/test/FitTest.C");
  MomRes(t);
}
