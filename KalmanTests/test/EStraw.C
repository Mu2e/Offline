void estraw(TTree* estraw,unsigned nmicro) {
  double nmicro_total(1.06e13); // total # of microbunches for the experiment
  double iongain(15); // number of ionizations per primary electron
  double eele(0.0025); // MPV straw energy deposition of an electron.
  double ggain(1e4); // gas amplification
  double qe(1.6e-19); // electron charge
  double phisym(12); // phi symmetry factor (in 1 station)
  double nlayer(2); // layer factor
  double factor = nmicro_total*iongain*ggain*qe/(nmicro*phisym*nlayer*eele);
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
  double iongain(15); // number of ionizations per primary electron
  double eele(0.0025); // MPV straw energy deposition of an electron.
  double ggain(1e4); // gas amplification
  double qe(1.6e-19); // electron charge (C)
  double phisym(12); // phi symmetry factor (in 1 station)
  double factor = iongain*ggain*qe*1e6/(nmicro*phisym*eele*10*2e-8);
  cout << "factor = " << factor << endl;
  TH1F* cvt = new TH1F("cvt","Current on wire 0+1, station <=4;time(nsec);#muA/wire",85,0,1700);
  estraw->Project("cvt","time","energy*(straw<2&&device<10)");
  cvt->Scale(factor);
  TCanvas* ccan = new TCanvas("ccan","background current",1200,800);
  ccan->Clear();
  ccan->Divide(1,1);
  ccan->cd(1);
  cvt->Draw();
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
  double nproton_microbunch(3.7e7); // # of protons on target/microbunch
  double iongain(15); // number of ionizations per primary electron
  double eele(0.0025); // MPV straw energy deposition of an electron.
  double ggain(1e4); // gas amplification
  double qe(1.6e-19); // electron charge (C)
  double phisym(12); // phi symmetry factor (in 1 station)
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

