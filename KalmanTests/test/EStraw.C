double nmicro_total(1.06e13); // total # of microbunches for the experiment
double nproton_microbunch(3.7e7); // # of protons on target/microbunch
double iongain(15); // number of ionizations per primary electron
double eele(0.0025); // MPV straw energy deposition of an electron.
double ggain(1e4); // gas amplification
double qe(1.6e-19); // electron charge (C)
unsigned nsect(6); // # sectors
unsigned nlayer(2); // # layers
unsigned ndevstat(2); // # devices/station
double mamp(1e6); // microamps/amp
double nsec(1e-9); // seconds/nanosecond
double mblength(1700); // time length of a microbunch (nanoseconds)
unsigned maxdevice(10); // maximum device to include
unsigned maxstraw(2); // maximum straw to include

void estraw(TTree* estraw,unsigned nmicro) {
  double factor = nmicro_total*iongain*ggain*qe/(nmicro*nsect*nlayer*ndevstat*eele);
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
    snprintf(stit,200,"Charge by panel station %i from muons;straw;Position from wire center (mm);C/cm/straw",istation);
    epanel[istation] = new TH2F(snam,stit,51,-0.5,50.5,120,-600,600);
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

void estraw_current(TTree* estraw,unsigned nmicro) {
  unsigned nbins(85); // number of time bins
  double factor = iongain*ggain*qe*mamp*nbins/(nmicro*nsect*nlayer*eele*maxdevice*maxstraw*mblength*nsec);
  cout << "factor = " << factor << endl;
  TH1F* cvt = new TH1F("cvt","Current on wire 0+1, station <=4;time(nsec);#muA/wire",nbins,0,mblength);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&device<%i)",maxstraw,maxdevice);
  TCut tcut(cut);
  estraw->Project("cvt","time","energy"*tcut);
  cvt->Scale(factor);
  TCanvas* ccan = new TCanvas("ccan","background current",1200,800);
  ccan->Clear();
  ccan->Divide(1,1);
  ccan->cd(1);
  cvt->Draw();
}

void estraw_rate(TTree* estraw,unsigned nmicro) {
  unsigned nbins(85); // number of time bins
  double factor = nbins/(nmicro*nsect*nlayer*maxdevice*maxstraw*mblength*nsec);
  cout << "factor = " << factor << endl;
  TH1F* hr = new TH1F("hr","Hit Rate on wire 0+1, station <=4;time(nsec);#hits/wire/sec",nbins,0,mblength);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&device<%i)",maxstraw,maxdevice);
  TCut tcut(cut);
  estraw->Project("hr","time",tcut);
  cvt->Scale(factor);
  TCanvas* rcan = new TCanvas("rcan","background hit rate",1200,800);
  rcan->Clear();
  rcan->Divide(1,1);
  rcan->cd(1);
  hr->Draw();
}



void estraw_flash(TTree* estraw,double nproton) {
  double nmicro_total(1.06e13); // total # of microbunches for the experiment
  double nproton_microbunch(3.7e7); // # of protons on target/microbunch
  double iongain(15); // number of ionizations per primary electron
  double eele(0.0025); // MPV straw energy deposition of an electron.
  double ggain(1e4); // gas amplification
  double qe(1.6e-19); // electron charge (C)
  double phisym(12); // phi symmetry factor (in 1 station)
  double nlayer(2); // number of layers in a panel
  double factor = nmicro_total*nproton_microbunch*iongain*ggain*qe/(nproton*phisym*nlayer*eele*12);
  cout << "factor = " << factor << endl;
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
    snprintf(stit,200,"Charge by panel station %i from beam flash;straw;Position from wire center (mm);C/cm/straw",istation);
    epanel[istation] = new TH2F(snam,stit,26,-1,51,20,-600,600);
    TCut tscut(scut);
// count energy/straw/cm
    estraw->Project(snam,"x:straw","energy"*tscut);
    epanel[istation]->Scale(factor);
    epanel[istation]->SetStats(0);
    if(ipad - npad*floor(ipad/npad) == 0){
      if(ecan != 0){
	char efcan[20];
	snprintf(efcan,20,"ecan_flash%i.png",ican);
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


void estraw_flash_current(TTree* estraw,double nproton) {

  double factor = nproton_microbunch*iongain*ggain*qe*1e6/(nproton*phisym*eele*10*2*2e-8);
  cout << "factor = " << factor << endl;
  TH1F* cvt = new TH1F("cvt","Current on wire 0+1, station <=4;time(nsec);#muA/wire",20,0,400);
  estraw->Project("cvt","time","energy*(straw<2&&device<10)");
  cvt->Scale(factor);
  TCanvas* fcan = new TCanvas("fcan","flash current",1200,800);
  fcan->Clear();
  fcan->Divide(1,1);
  fcan->cd(1);
  cvt->Draw();
}

