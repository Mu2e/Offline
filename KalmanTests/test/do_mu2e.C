{  
  gROOT->LoadMacro("Scripts/FillChain.C");
  gROOT->LoadMacro("KalmanTests/test/mu2e.C+");
  TChain* con;
  TChain* dio;
//  char conname[100] = {"/data/HD2/brownd/mu2e/mixConversion.8675879"};
//  char dioname[100] = {"/data/HD2/brownd/mu2e/mixDIO.8675884"};
  char conname[100] = {"/data/mixConversion.8684286"};
  char conname2[100] = {"/data/mixConversion.8746909"};
  char dioname[100] = {"/data/mixDIO.8684287"};
  char dioname2[100] = {"/data/mixDIO.8711081"};
  char dioname3[100] = {"/data/mixDIO.8736887"};
  double momlow(103.4);
  double momhigh(104.8);
  cout << conname << endl;
  FillChain(con,conname,"mixConversion.root",10);
  FillChain(con,conname2,"mixConversion.root",30);
  FillChain(dio,dioname,"mixDIO.root",50);
  FillChain(dio,dioname2,"mixDIO.root",50);
  FillChain(dio,dioname3,"mixDIO.root",100);
  mu2e* m2e = new mu2e(dio,con,10.0,4000000,400000);
//  mu2e* m2e = new mu2e(dio,con,10.0,1000000,100000);
  m2e->fillmu2e(251,98.0,106.0);
  m2e->drawmu2e(momlow,momhigh);
  m2e->fitReco(2);
  m2e->smearDIO(100000,100);
  m2e->doExperiments(103.4,104.8,1e-16,2,1000,0);
}

