{  
  gROOT->LoadMacro("Scripts/FillChain.C");
  gROOT->LoadMacro("KalmanTests/test/mu2e.C+");
  TChain* con;
  TChain* dio;
//  char conname[100] = {"/data/HD2/brownd/mu2e/mixConversion.8675879"};
//  char dioname[100] = {"/data/HD2/brownd/mu2e/mixDIO.8675884"};
  char conname[100] = {"/data/HD2/brownd/mu2e/mixConversion.8684286"};
  char conname2[100] = {"/data/HD2/brownd/mu2e/mixConversion.8746909"};
  char dioname[100] = {"/data/HD2/brownd/mu2e/mixDIO.8684287"};
  char dioname2[100] = {"/data/HD2/brownd/mu2e/mixDIO.8711081"};
  char dioname3[100] = {"/data/HD2/brownd/mu2e/mixDIO.8736887"};
  double momlow(103.4);
  double momhigh(104.8);
  cout << conname << endl;
  FillChain(con,conname,"mixConversion.root",10);
  FillChain(con,conname2,"mixConversion.root",30);
  FillChain(dio,dioname,"mixDIO.root",50);
  FillChain(dio,dioname2,"mixDIO.root",50);
  FillChain(dio,dioname3,"mixDIO.root",100);
  mu2e* m2e = new mu2e(dio,con,10.0,4000000,400000);
  m2e->drawmu2e(momlow,momhigh);
}

