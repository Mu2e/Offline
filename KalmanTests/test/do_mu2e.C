{  
  gROOT->LoadMacro("Scripts/FillChain.C");
  gROOT->LoadMacro("KalmanTests/test/mu2e.C+");
  gROOT->LoadMacro("KalmanTests/test/KalFit.C+");
  TChain* con(0);
  TChain* dio(0);
  char condir[50] = {"/data/mixc_4.5441689"};
  char diodir[50] = {"/data/mixDIO4.5473031"};
  FillChain(con,condir,"mixConversion.root",10);
  FillChain(dio,diodir,"mixDIO.root",40);
  mu2e(dio,con,10.0,1000000,100000,103.4,104.8,true,"_4sig.png");
  KalCuts();
  KalFitAcc(con);
  acan->SaveAs("acc_4sig.png");
}

