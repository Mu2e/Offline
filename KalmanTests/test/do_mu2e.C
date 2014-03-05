{
  gROOT->LoadMacro("Scripts/FillChain.C");
  gROOT->LoadMacro("KalmanTests/test/mu2e.C");
  TChain* con;
  TChain* dio;
  FillChain(con,"/data/MC.15285782","MC.root",10);
  FillChain(dio,"/data/MDIO.15308328/","MDIO.root",50);
  mu2e* m2e = new mu2e(dio,con,10.0,10000,10000);
  m2e->fillmu2e(251,98.0,106.0);
  double momlow(103.4);
  double momhigh(104.8);
  m2e->drawmu2e(momlow,momhigh,true);
//  m2e->fitReco(2);
//  m2e->smearDIO(100000,100);
//  m2e->doExperiments(103.4,104.8,1e-16,2,1000,0);
}

