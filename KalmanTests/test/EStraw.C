double nmicro_total(1.06e13); // total # of microbunches for the experiment
double iongain(15); // number of ionizations per primary electron
double eele(0.0025); // MPV straw energy deposition of an electron (MeV)
double ggain(1e4); // gas amplification
double qe(1.6e-19); // electron charge (C)
unsigned npanel(6); // # panels
unsigned nlayer(2); // # layers
unsigned ndevstat(2); // # planes/station
double mamp(1e6); // microamps/amp
double secpernsec(1e-9); // seconds/nanosecond
double mblength(1700); // time length of a microbunch (nanoseconds)
unsigned maxplane(10); // maximum plane to include
unsigned maxstraw(2); // maximum straw to include
unsigned ntbins(50);
unsigned nxbins(30);
double xrange(600.0); // cm

void estraw(TTree* estraw,double nmicro) {
  double factor = nmicro_total*iongain*ggain*qe*nxbins*10/(nmicro*npanel*nlayer*ndevstat*eele*2*xrange);
  TH2F* epanel[18];
  unsigned ican(0);
  unsigned ipad(0);
  unsigned npad(4);
  TCanvas* ecan(0);
  for(unsigned istation = 0;istation<18;istation++){
    char scut[100];
    char snam[20];
    char stit[200];
    snprintf(scut,100,"(plane==%i||plane==%i)",2*istation,2*istation+1);
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
  double factor = iongain*ggain*qe*mamp*ntbins/(nmicro*npanel*nlayer*eele*maxplane*maxstraw*mblength*secpernsec);
  double factor2 = iongain*ggain*qe*mamp*ntbins*nxbins*10/(nmicro*npanel*nlayer*eele*maxplane*maxstraw*mblength*secpernsec*2*xrange);
//  cout << "factor = " << factor << endl;
  TH1F* cvt = new TH1F("cvt","Current on wire 0+1, station <=4;time(nsec);#muA/wire",ntbins,0,mblength);
  TH2F* cvt2 = new TH2F("cvt2","Current on wire 0+1, station <=4;time(nsec);Position WRT wire center (mm);#muA/wire/cm",ntbins,0,mblength,
  nxbins,-xrange,xrange);
  cvt2->Sumw2();
  cvt2->SetStats(0);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&plane<%i)",maxstraw,maxplane);
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
  double factor = (1e-6)*ntbins/(nmicro*npanel*nlayer*maxplane*maxstraw*mblength*secpernsec);
  double factor2 = (1e-6)*ntbins*nxbins*10/(nmicro*npanel*nlayer*maxplane*maxstraw*mblength*secpernsec*2*xrange);
  cout << "factor = " << factor << endl;
  TH1F* hr = new TH1F("hr","Hit Rate on wire 0+1, station <=4;time(nsec);MHz/wire",ntbins,0,mblength);
  TH2F* hr2 = new TH2F("hr2","Hit Rate on wire 0+1, station <=4;time(nsec);Position WRT wire center (mm);MHz/wire/cm",ntbins,0,mblength,
  nxbins,-xrange,xrange);
  hr2->Sumw2();
  hr2->SetStats(0);
  char cut[80];
  snprintf(cut,80,"(straw<%i&&plane<%i)",maxstraw,maxplane);
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

void estraw_deadtime(TTree* estraw,double nmicro) {

  TH1F* eloge = new TH1F("eloge","Background straw hit deposited energy;log_{10}(MeV);hits/event",100,-4.5,0);
  TH1F* ploge = new TH1F("ploge","Background straw hit deposited energy;log_{10}(MeV);hits/event",100,-4.5,0);
  TH1F* gloge = new TH1F("gloge","Background straw hit deposited energy;log_{10}(MeV);hits/event",100,-4.5,0);
  TH1F* nloge = new TH1F("nloge","Background straw hit deposited energy;log_{10}(MeV);hits/event",100,-4.5,0);
  eloge->SetLineColor(kRed);
  ploge->SetLineColor(kBlue);
  gloge->SetLineColor(kCyan);
  nloge->SetLineColor(kOrange);
  eloge->SetStats(0);
  ploge->SetStats(0);
  gloge->SetStats(0);
  nloge->SetStats(0);

  estraw->Project("eloge","log10(energy)","abs(mcpdg)==11");
  estraw->Project("ploge","log10(energy)","mcpdg==2212");
  estraw->Project("gloge","log10(energy)","mcpdg==22");
  estraw->Project("nloge","log10(energy)","mcpdg==2112");

  double factor=1.0/nmicro;
  eloge->Scale(factor);
  ploge->Scale(factor);
  gloge->Scale(factor);
  nloge->Scale(factor);

  double scale=ggain*iongain*qe/eele;

  TH1F* elogq = new TH1F("elogq","Background straw hit deposited charge;nC;hits/event",500,-0.0001,0.005);
  TH1F* plogq = new TH1F("plogq","Background straw hit deposited charge;nC;hits/event",500,-0.0001,0.005);
  TH1F* glogq = new TH1F("glogq","Background straw hit deposited charge;nC;hits/event",500,-0.0001,0.005);
  TH1F* nlogq = new TH1F("nlogq","Background straw hit deposited charge;nC;hits/event",500,-0.0001,0.005);
  elogq->SetLineColor(kRed);
  plogq->SetLineColor(kBlue);
  glogq->SetLineColor(kCyan);
  nlogq->SetLineColor(kOrange);
  elogq->SetStats(0);
  plogq->SetStats(0);
  glogq->SetStats(0);
  nlogq->SetStats(0);


// conversion to nano-Coulombs
  estraw->Project("elogq","energy*9.6e-3","abs(mcpdg)==11");
  estraw->Project("plogq","energy*9.6e-3","mcpdg==2212");
  estraw->Project("glogq","energy*9.6e-3","mcpdg==22");
  estraw->Project("nlogq","energy*9.6e-3","mcpdg==2112");

  double factor=1.0/nmicro;
  elogq->Scale(factor);
  plogq->Scale(factor);
  glogq->Scale(factor);
  nlogq->Scale(factor);

// JPs data
  double fc[7] = {1,2,5,7,10,15,20};
  double dead[7] = {60,210,370,410,450,490,520};
  TGraph* jp = new TGraph(7,fc,dead);

  TCanvas* dcan = new TCanvas("dcan","dcan",1200,800);
  dcan->Divide(2,2);
  dcan->cd(1);
  gPad->SetLogy();
  eloge->Draw();
  ploge->Draw("same");
  gloge->Draw("same");
  nloge->Draw("same");
  dcan->cd(2);
  gPad->SetLogy();
  elogq->Draw();
  plogq->Draw("same");
  glogq->Draw("same");
  nlogq->Draw("same");
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(eloge,"Electrons","l");
  leg->AddEntry(ploge,"Protons","l");
  leg->AddEntry(gloge,"gammas","l");
  leg->AddEntry(nloge,"Neutrons","l");
  leg->Draw(); 

}

