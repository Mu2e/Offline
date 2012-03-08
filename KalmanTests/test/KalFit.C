
TCut tpitch("mcenttd>0.5&&mcenttd<1.1");
TCut tt0("mct0>700.0");
TCut tmom("mcentmom>100");
TCut nmch("nchits>=20");

TCut reco("fitstatus>0");
TCut livegate("t0>710");
TCut rpitch("td>0.57735&&td<1.0");
TCut goodfit("fitcon>1e-4&&nactive>=20&&t0err<1.5&&fitmomerr<0.2");
TCut cosmic("abs(d0)<105 && d0+2/om>460 && d0+2/om<660");
TCut rmom("fitmom>103.5&&fitmom<104.7");

Double_t splitgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Double_t core;
  Double_t tail;
  Float_t xval = x[0];
  if(xval > par[1]) {
    core = exp(-0.5*pow((xval-par[1])/par[2],2))/par[2];
    tail = par[4]*exp(-0.5*pow((xval-par[1])/par[5],2))/par[5];
  } else {
    core = exp(-0.5*pow((xval-par[1])/par[3],2))/par[3];
    tail = (1/par[2]-1/par[3]+par[4]/par[5])*exp(-0.5*pow((xval-par[1])/par[6],2));
  }
  retval = par[0]*0.398942*(core+tail);
// add a tail Gaussian
  return retval;
}

void KalFitHit (TTree* hits ) {
    TCanvas* dcan = new TCanvas("driftcan","driftcan",1200,800);
    TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
    TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-10,10);
    TH2F* drad = new TH2F("drad","Drift radius;true drift radius (mm);reco drift radius (mm)",
      100,-0.3,2.8,100,-0.3,2.8);
    TH1F* rpull = new TH1F("rpull","residual pull",100,-10,10);
    hits->Project("dres","rdrift-mcrdrift","active");
    hits->Project("dpull","(rdrift-mcrdrift)/rdrifterr","active");
    hits->Project("rpull","resid/residerr","active");
    dcan->Clear();
    dcan->Divide(2,2);
    dcan->cd(1);
    hits->Draw("rdrift:mcrdrift>>drad","active");
    dcan->cd(2);
    dres->Fit("gaus");
    dcan->cd(3);
    dpull->Fit("gaus");
    dcan->cd(4);
    rpull->Fit("gaus");
    
    TCanvas* tcan = new TCanvas("ht0can","hit_t0can",1200,800);
    TH1F* t0res = new TH1F("t0res","hit t0 resolution;nsec",100,-10,10);
    TH1F* t0pull = new TH1F("t0pull","hit t0 pull",100,-10,10);
    TH2F* dt0 = new TH2F("dt0","Hit t0;true t0 (nsec);reco t0 (nsec)",
        100,500,4000,100,500,4000);
    hits->Project("t0res","hitt0-mchitt0","active");
    hits->Project("t0pull","(hitt0-mchitt0)/hitt0err","active");
    tcan->Clear();
    tcan->Divide(2,2);
    tcan->cd(1);
    hits->Draw("hitt0:mchitt0>>dt0","active");
    tcan->cd(2);
    t0res->Fit("gaus");
    tcan->cd(3);
    t0pull->Fit("gaus");
    
    TCanvas* tdcan = new TCanvas("tdcan","tdcan",1200,800);
    TH1F* tdres = new TH1F("tdres","#Deltat V resolution;mm",100,-200,200);
    TH1F* tdpull = new TH1F("tdpull","#Deltat V pull",100,-10,10);
    TH2F* dtd = new TH2F("dtd","Hit V position;true V (mm);#Deltat V (mm)",
      100,-600,600,100,-600,600);
    TH2F* pocatd = new TH2F("pocatd","Hit POCA V;true V (mm);POCA V (mm)",
      100,-600,600,100,-600,600);
        
    hits->Project("tdres","dmid-mcdmid","active");
    hits->Project("tdpull","(dmid-mcdmid)/dmiderr","active");
    tdcan->Clear();
    tdcan->Divide(2,2);
    tdcan->cd(1);
    hits->Draw("dmid:mcdmid>>dtd","active");
    tdcan->cd(2);
    hits->Draw("hflt:mcdmid>>pocatd","active");
    tdcan->cd(3);
    tdres->Fit("gaus");
    tdcan->cd(4);
    tdpull->Fit("gaus");
    
}

void KalFitTrk (TTree* trks ) {
  
  TCanvas* tcan = new TCanvas("tt0can","trk_t0can",1200,800);
  TH1F* t00res = new TH1F("t00res","Initial t0 resolution;nsec",100,-20,20);
  TH1F* t0res = new TH1F("t0res","Final t0 resolution;nsec",100,-10,10);
  TH1F* t0pull = new TH1F("t0pull","Track t0 pull",100,-10,10);
  TH2F* dt0 = new TH2F("dt0","Track t0;true t0 (nsec);Initial t0 (nsec)",
    100,500,4000,100,500,4000);
  trks->Project("t00res","t00-mct0","fitstatus>0");
  trks->Project("t0res","t0-mct0","fitstatus>0");
  trks->Project("t0pull","(t0-mct0)/t0err","fitstatus>0");
  tcan->Clear();
  tcan->Divide(2,2);
  tcan->cd(1);
  trks->Draw("t00:mct0>>dt0","fitstatus>0");
  tcan->cd(2);
  t00res->Fit("gaus");
  tcan->cd(3);
  t0res->Fit("gaus");
  tcan->cd(4);
  t0pull->Fit("gaus");
  
  
  TCanvas* pcan = new TCanvas("pullcan","pullcan",1200,800);
  TH1F* d0pull = new TH1F("d0pull","d0 pull",100,-10,10);
  TH1F* p0pull = new TH1F("p0pull","#phi0 pull",100,-10,10);
  TH1F* ompull = new TH1F("ompull","#omega pull",100,-10,10);
  TH1F* z0pull = new TH1F("z0pull","z0 pull",100,-10,10);
  TH1F* tdpull = new TH1F("tdpull","tan(#lambda) pull",100,-10,10);
  trks->Project("d0pull","(d0-mcd0)/d0err","fitstatus>0");
  trks->Project("p0pull","(p0-mcp0)/p0err","fitstatus>0");
  trks->Project("ompull","(om-mcom)/omerr","fitstatus>0");
  trks->Project("z0pull","(z0-mcz0)/z0err","fitstatus>0");
  trks->Project("tdpull","(td-mctd)/tderr","fitstatus>0");
  pcan->Clear();
  pcan->Divide(3,2);
  pcan->cd(1);
  d0pull->Fit("gaus");
  pcan->cd(2);
  p0pull->Fit("gaus");
  pcan->cd(3);
  ompull->Fit("gaus");
  pcan->cd(4);
  z0pull->Fit("gaus");
  pcan->cd(5);
  tdpull->Fit("gaus");

  TCanvas* fcan = new TCanvas("fitcan","fitcan",1200,800);
  TH2F* mom = new TH2F("mom","momentum at first hit;true mom (MeV);reco mom (MeV)",
    100,80,110,100,80,110);
  TH1F* mres = new TH1F("mres","momentum resolution at first hit;MeV",100,-2,2);
  TH1F* mpull = new TH1F("mpull","momentum pull at first hit",100,-10,10);
  TH1F* chisq = new TH1F("chisq","Chisq/NDof",100,0,10);
  trks->Project("mres","fitmom-mcmom","fitstatus>0");
  trks->Project("mpull","(fitmom-mcmom)/fitmomerr","fitstatus>0");
  trks->Project("chisq","chisq/ndof","fitstatus>0");
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  trks->Draw("fitmom:mcmom>>mom","fitstatus>0");
  fcan->cd(2);
  mres->Draw();
  fcan->cd(3);
  mpull->Fit("gaus");
  fcan->cd(4);
  chisq->Draw();
}

void KalFitAccPlots(TTree* trks) {

  TH1F* nmc = new TH1F("nmc","N Straw Hits from CE;N straws",81,-0.5,80.5);
  TH1F* mcmom = new TH1F("mcmom","CE true momentum at tracker;CE momentum (MeV)",57,49,106);
  
  TH1F* fitcon = new TH1F("fitcon","log_{10} fit consistency",101,-8,0);
  TH1F* momerr = new TH1F("momerr","Fit momentum error;momentum error (MeV)",100,0,0.5);
  TH1F* t0err = new TH1F("t0err","Fit t_{0} error; t_{0} error (ns)",100,0,2.0);
  TH1F* na = new TH1F("na","Fit N active hits;N hits",71,-0.5,70.5);

  TH1F* t0 = new TH1F("t0","Track Fit t_{0};t_{0} (ns)",100,300,1800);
  TH1F* td = new TH1F("td","Track Fit tan(#lambda);tan(#lambda)",100,0.5,1.5);
  TH1F* d0 = new TH1F("d0","Track fit d_{0};d_{0} (mm)",100,-150,150);
  TH1F* rmax = new TH1F("rmax","Track fit rmax;d_{0}+2/#omega (mm)",100,300,900);

  TH1F* fitmom = new TH1F("fitmom","Track fit momentum;fit momentum (MeV)",100,98,107);

  trks->Project("nmc","nchits");
  trks->Project("mcmom","mcentmom",nmch);
  
  trks->Project("fitcon","log10(fitcon)",reco+nmch+tmom);
  trks->Project("momerr","fitmomerr",reco+nmch+tmom);
  trks->Project("t0err","t0err",reco+nmch+tmom);
  trks->Project("na","nactive",reco+nmch+tmom);

  trks->Project("t0","t0",reco+nmch+tmom+goodfit);
  trks->Project("td","td",reco+nmch+tmom+goodfit+livegate);

  trks->Project("d0","d0",reco+nmch+tmom+goodfit+livegate);
  trks->Project("rmax","d0+2.0/om",reco+nmch+tmom+goodfit+livegate);

  trks->Project("fitmom","fitmom",reco+nmch+tmom+goodfit+livegate+cosmic);

  TCanvas* pcan = new TCanvas("pcan","Pre-acceptance",1200,800);
  pcan->Clear();
  pcan->Divide(1,2);
  pcan->cd(1);
  gPad->SetLogy();
  nmc->Draw();
  TLine* nmccut = new TLine(20,0.0,20,nmc->GetMaximum());
  nmccut->SetLineColor(kBlack);
  nmccut->SetLineStyle(2);
  nmccut->SetLineWidth(2);
  nmccut->Draw();

  pcan->cd(2);
  gPad->SetLogy();
  mcmom->Draw();
  TLine* tmomcut = new TLine(100,0.0,100,mcmom->GetMaximum());
  tmomcut->SetLineColor(kBlack);
  tmomcut->SetLineStyle(2);
  tmomcut->SetLineWidth(2);
  tmomcut->Draw();

  TCanvas* fcan = new TCanvas("fcan","Fit Acceptance",1200,800);
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  fitcon->Draw();
  TLine* fitconcut = new TLine(-4,0.0,-4,fitcon->GetMaximum());
  fitconcut->SetLineColor(kBlack);
  fitconcut->SetLineStyle(2);
  fitconcut->SetLineWidth(2);
  fitconcut->Draw();

  fcan->cd(2);
  na->Draw();
  TLine* nacut = new TLine(20,0.0,20,na->GetMaximum());
  nacut->SetLineColor(kBlack);
  nacut->SetLineStyle(2);
  nacut->SetLineWidth(2);
  nacut->Draw();

  fcan->cd(3);
  t0err->Draw();
  TLine* t0errcut = new TLine(1.5,0.0,1.5,t0err->GetMaximum());
  t0errcut->SetLineColor(kBlack);
  t0errcut->SetLineStyle(2);
  t0errcut->SetLineWidth(2);
  t0errcut->Draw();
  
  fcan->cd(4);
  momerr->Draw();
  TLine* momerrcut = new TLine(0.2,0.0,0.2,momerr->GetMaximum());
  momerrcut->SetLineColor(kBlack);
  momerrcut->SetLineStyle(2);
  momerrcut->SetLineWidth(2);
  momerrcut->Draw();

  TCanvas* tcan = new TCanvas("tcan","Track parameter acceptance",1200,800);
  tcan->Divide(2,2);
  tcan->cd(1);
  t0->Draw();
  TLine* t0cut = new TLine(710,0.0,710,t0->GetMaximum());
  t0cut->SetLineColor(kBlack);
  t0cut->SetLineStyle(2);
  t0cut->SetLineWidth(2);
  t0cut->Draw();

  tcan->cd(2);
  td->Draw();
  TLine* tdcut_l = new TLine(0.57735,0.0,0.57735,td->GetMaximum());
  tdcut_l->SetLineColor(kBlack);
  tdcut_l->SetLineStyle(2);
  tdcut_l->SetLineWidth(2);
  tdcut_l->Draw();
  TLine* tdcut_h = new TLine(1.0,0.0,1.0,td->GetMaximum());
  tdcut_h->SetLineColor(kBlack);
  tdcut_h->SetLineStyle(2);
  tdcut_h->SetLineWidth(2);
  tdcut_h->Draw();
  
  tcan->cd(3);
  d0->Draw();
  TLine* d0cut = new TLine(105,0.0,105,d0->GetMaximum());
  d0cut->SetLineColor(kBlack);
  d0cut->SetLineStyle(2);
  d0cut->SetLineWidth(2);
  d0cut->Draw();

  tcan->cd(4);
  rmax->Draw();
  TLine* rmaxcut = new TLine(660,0.0,660,rmax->GetMaximum());
  rmaxcut->SetLineColor(kBlack);
  rmaxcut->SetLineStyle(2);
  rmaxcut->SetLineWidth(2);
  rmaxcut->Draw();
 
  TCanvas* mcan = new TCanvas("mcan","momentum",800,600);
  mcan->Divide(1,1);
  mcan->cd(1);
  fitmom->Draw();
  TLine* fitmomcut_l = new TLine(103.5,0.0,103.5,fitmom->GetMaximum());
  fitmomcut_l->SetLineColor(kBlack);
  fitmomcut_l->SetLineStyle(2);
  fitmomcut_l->SetLineWidth(2);
  fitmomcut_l->Draw();
  TLine* fitmomcut_h = new TLine(104.7,0.0,104.7,fitmom->GetMaximum());
  fitmomcut_h->SetLineColor(kBlack);
  fitmomcut_h->SetLineStyle(2);
  fitmomcut_h->SetLineWidth(2);
  fitmomcut_h->Draw();
 

} 

void KalFitAcc(TTree* trks) {
  unsigned nbins(9);
  double bmax = nbins-0.5;
  TH1F* acc = new TH1F("acc","CE Acceptance;;cummulative acceptance",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","CE Acceptance;;relative acceptance",nbins,-0.5,bmax);
//  acc->Sumw2();
//  racc->Sumw2();
  unsigned i(1);
  acc->GetXaxis()->SetBinLabel(1,"All CE");
  acc->GetXaxis()->SetBinLabel(2,">=20 CE SH");
  acc->GetXaxis()->SetBinLabel(3,"CE p>100 MeV/c");
//  acc->GetXaxis()->SetBinLabel(4,"CE pitch");
  acc->GetXaxis()->SetBinLabel(4,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(5,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(6,"Livegate");
  acc->GetXaxis()->SetBinLabel(7,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(8,"Cosmic Rejection");
  acc->GetXaxis()->SetBinLabel(9,"Momentum window");

  racc->GetXaxis()->SetBinLabel(1,"All CE");
  racc->GetXaxis()->SetBinLabel(2,">=20 CE SH");
  racc->GetXaxis()->SetBinLabel(3,"CE p>100 MeV/c");
//  racc->GetXaxis()->SetBinLabel(4,"CE pitch");
  racc->GetXaxis()->SetBinLabel(4,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(5,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(6,"Livegate");
  racc->GetXaxis()->SetBinLabel(7,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(8,"Cosmic Rejection");
  racc->GetXaxis()->SetBinLabel(9,"Momentum window");
  
  
  trks->Project("acc","0.0");
  trks->Project("+acc","1.0",nmch);
  trks->Project("+acc","2.0",nmch+tmom);
//  trks->Project("+acc","3.0",nmch+tmom+tpitch);
  trks->Project("+acc","3.0",nmch+tmom+reco);
  trks->Project("+acc","4.0",nmch+tmom+reco+goodfit);
  trks->Project("+acc","5.0",nmch+tmom+reco+goodfit+livegate);
  trks->Project("+acc","6.0",nmch+tmom+reco+goodfit+livegate+rpitch);
  trks->Project("+acc","7.0",nmch+tmom+reco+rpitch+livegate+cosmic+goodfit);
  trks->Project("+acc","8.0",nmch+tmom+reco+rpitch+livegate+cosmic+goodfit+rmom);

  double all = acc->GetBinContent(1);
  double prev = all;
  for(unsigned ibin=1;ibin<=nbins;++ibin){
    racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    prev = acc->GetBinContent(ibin);
  }
  racc->SetMaximum(1.1);
  acc->Scale(1.0/all);
  acc->SetMaximum(1.1);
  acc->GetXaxis()->SetLabelSize(0.06);
  racc->GetXaxis()->SetLabelSize(0.06);
  acc->SetMarkerSize(2.0);
  racc->SetMarkerSize(2.0);
  acc->GetYaxis()->SetTitleSize(0.05);
  racc->GetYaxis()->SetTitleSize(0.05);

  gStyle->SetPaintTextFormat("5.4f");
  TCanvas* acan = new TCanvas("acan","Acceptance",1200,800);
  acan->Clear();
  acan->Divide(1,2);
  acan->cd(1);
  TPad* tp = (TPad*)acan->cd(1);
  tp->SetBottomMargin(0.15);
  acc->Draw("histtext0");
  acan->cd(2);
  tp = (TPad*)acan->cd(2);
  tp->SetBottomMargin(0.15);
  racc->Draw("histtext0");
}

void KalFitRes(TTree* trks) {
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
 
  TH1F* momres[4];
// basic cuts
  TCut reco("fitstatus>0");
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(710);
  char ctext[80];
  snprintf(ctext,80,"td>%f&&td<%f",tdlow,tdhigh);
  TCut pitch(ctext);
  snprintf(ctext,80,"t0>%f",t0min);
  TCut livegate(ctext);
// mc selection cuts for efficiency
  TCut mcsel("mcentmom>100&&mcenttd<1.0&&mcenttd>0.5774&&nchits>=20&&mct0>670");
// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = "nactive>=20";
  ncuts[1] = "nactive>=20";
  ncuts[2] = "nactive>=25";
  ncuts[3] = "nactive>=30";
  t0cuts[0] = "";
  t0cuts[1] = "t0err<1.5";
  t0cuts[2] = "t0err<1.0";
  t0cuts[3] = "t0err<0.9";
  momcuts[0] = "";
  momcuts[1] = "fitmomerr<0.2";
  momcuts[2] = "fitmomerr<0.18";
  momcuts[3] = "fitmomerr<0.15";
  fitcuts[0] = "";
  fitcuts[1] = "fitcon>1e-4";
  fitcuts[2] = "fitcon>1e-3";
  fitcuts[3] = "fitcon>1e-2";

  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","mcentmom",mcsel);
 
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
  rcan->Clear();
  rcan->Divide(2,2);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat(1111111);
  for(unsigned ires=0;ires<4;ires++){
    rcan->cd(ires+1);
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",151,-2.5,2.5);
//  momres[ires]->SetStats(0);
    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
    TCut final = (reco+pitch+livegate+quality);
    trks->Project(mname,"fitmom-mcentmom",final);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momres[ires]->GetRMS(),momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
    sgau->SetParLimits(5,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.8);
    momres[ires]->Fit("sgau","L");


    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
    TPaveText* rtext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    rtext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    rtext->AddText(line);
    sprintf(line,"%s",ncuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"Eff=%3.2f",keff);
    rtext->AddText(line);
    rtext->Draw();
 
  }
  rcan->cd(0);
}


