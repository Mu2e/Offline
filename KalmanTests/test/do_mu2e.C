{  
  gROOT->LoadMacro("Scripts/FillChain.C");
  gROOT->LoadMacro("KalmanTests/test/mu2e.C+");
  TChain* con(0);
  TChain* dio(0);
  char condir[50] = {"/data/mixConversion.6418952"};
  char diodir[50] = {"/data/mixDIO.6430701"};
  FillChain(con,condir,"mixConversion.root",10,"ReadKalFits/trkdiag");
  FillChain(dio,diodir,"mixDIO.root",40,"ReadKalFits/trkdiag");
  mu2e(dio,con,10.0,1000000,100000,103.4,104.8,true);
}

