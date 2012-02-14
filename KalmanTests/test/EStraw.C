double nmicro_total(1.06e13); // total # of microbunches for the experiment
double iongain(15); // number of ionizations per primary electron
double eele(0.0025); // MPV straw energy deposition of an electron (MeV)
double ggain(1e4); // gas amplification
double qe(1.6e-19); // electron charge (C)
unsigned nsect(6); // # sectors
unsigned nlayer(2); // # layers
unsigned ndevstat(2); // # devices/station
double mamp(1e6); // microamps/amp
double secpernsec(1e-9); // seconds/nanosecond
double mblength(1700); // time length of a microbunch (nanoseconds)
unsigned maxdevice(10); // maximum device to include
unsigned maxstraw(2); // maximum straw to include
unsigned ntbins(50);
unsigned nxbins(30);
double xrange(600.0); // cm

void estraw(TTree* estraw,double nmicro) {
  double factor = nmicro_total*iongain*ggain*qe*nxbins*10/(nmicro*nsect*nlayer*ndevstat*eele*2*xrange);
  TH2F* epanel[18];
  unsigned ican(0);
  unsigned ipad(0);
  unsigned npad(4);
  TCanvas* ecan(0);
  for(unsigned istation = 0;istation<18;istation++){
    char scut[100];
    char snam[20];
    char stit[200];
    snprintf(scut,100,"(device==%i||device==%i)",2*istation,2*istation+1);
    snprintf(snam,20,"epanel%i",istation);
    snprintf(stit,200,"Charge by panel station %i;straw;Position WRT wire center (mm);C/cm/straw",istation);
    epanel[istation] = new TH2F(snam,stit,51,-0.5,50.5,nxbins,-xrange,xrange);
    TCut tscut(scut);
    // count energy/straw/cm
    estraw->Project(snam,"x:straw","energy"*tscut);
    epanel[istation]->Scale(factor);
    epanel[istation]->SetStats(0);
    if(ipad - npad*floor(ipad/npad) == 0){
      if(ecan != 0){
	char efcan[20];
	snprintf(efcan,20,"ecan%i.png",ican);
	ecan->SaveAs(efcan);
      }
      ican++;
      char cnam[20];
      snprintf(cnam,20,"ecan%i",ican);
      ecan = new TCanvas(cnam,"ecan",1200,800);
      ecan->Clear();
      ecan->Divide(2,2);
      ipad = 1;
    } else
      ipad++;
    TPad* tp = (TPad*)ecan->cd(ipad);
    tp->SetRightMargin(0.15);
    epanel[istation]->Draw("CONTZ");
  }
}

void estraw_current(TTree* estraw,double nmicro) {
  double factor = iongain*ggain*qe*mamp*ntbins/(nmicro*nsect*nlayer*eele*maxdevice*maxstraw*mblength*secpernsec);
  double factor2 = iongain*ggain*qe*mamp*ntbins*nxbins*10/(nmicro*nsect*nlayer*eele*maxdevice*maxstraw*mblength*secpernsec*2*xrange);
//  cout << "factor = " << factor << endl;
  TH1F* cvt = new TH1F("cvt","Current on wire 0+1, station <=4;time(nsec);#muA/wire",ntbins,0,mblength);
  TH2F* cvt2 = new TH2F("cvt2","Current on wire 0+1, station <=4;time(nsec);Position WRT wire center (mm);#muA/wire/cm",ntbins,0,mblength,
  nxbins,-xrange,xrange);
  cvt2->Sumw2();
  cvt2->SetStats(0);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&device<%i)",maxstraw,maxdevice);
  TCut tcut(cut);
//  cout << "cut = " << cut << endl;
  estraw->Project("cvt","time","energy"*tcut);
  estraw->Project("cvt2","x:time","energy"*tcut);
  cvt->Scale(factor);
  cvt2->Scale(factor2);
  TCanvas* ccan = new TCanvas("ccan","background current",1200,800);
  ccan->Clear();
  ccan->Divide(1,2);
  ccan->cd(1);
  cvt->Draw();
  TPad* tp = (TPad*)ccan->cd(2);
  tp->SetRightMargin(0.15);
  cvt2->Draw("CONTZ");
}

void estraw_rate(TTree* estraw,double nmicro) {
  double factor = (1e-6)*ntbins/(nmicro*nsect*nlayer*maxdevice*maxstraw*mblength*secpernsec);
  double factor2 = (1e-6)*ntbins*nxbins*10/(nmicro*nsect*nlayer*maxdevice*maxstraw*mblength*secpernsec*2*xrange);
  cout << "factor = " << factor << endl;
  TH1F* hr = new TH1F("hr","Hit Rate on wire 0+1, station <=4;time(nsec);MHz/wire",ntbins,0,mblength);
  TH2F* hr2 = new TH2F("hr2","Hit Rate on wire 0+1, station <=4;time(nsec);Position WRT wire center (mm);MHz/wire/cm",ntbins,0,mblength,
  nxbins,-xrange,xrange);
  hr2->Sumw2();
  hr2->SetStats(0);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&device<%i)",maxstraw,maxdevice);
  TCut tcut(cut);
  estraw->Project("hr","time",tcut);
  hr->Scale(factor);
  estraw->Project("hr2","x:time",tcut);
  hr2->Scale(factor2);
  TCanvas* rcan = new TCanvas("rcan","background hit rate",1200,800);
  rcan->Clear();
  rcan->Divide(1,2);
  rcan->cd(1);
  hr->Draw();
  TPad* tp = (TPad*)rcan->cd(2);
  tp->SetRightMargin(0.15);
  hr2->Draw("CONTZ");
}

