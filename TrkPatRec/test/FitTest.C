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

void KalTest (TTree* trk) {
  TCanvas* kcan = new TCanvas("kcan","Kalman Fit",1200,800);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TH1F* t0 = new TH1F("t0","t0 resolution;nsec",100,-5,5);
  trk->Project("t0","t0-mcmidt0","kalfail==0&&nactive>=20");
  TH1F* nh = new TH1F("nh","N hits",66,-0.5,65.5);
  TH1F* na = new TH1F("na","N hits",66,-0.5,65.5);
  TH1F* nd = new TH1F("nd","N hits",66,-0.5,65.5);
  nh->SetLineColor(kBlack);
  na->SetLineColor(kRed);
  nd->SetLineColor(kGreen);
  TH1F* fstat = new TH1F("fstat","fit status",22,-1.5,20.5);
  TH1F* momr = new TH1F("momr","momentum resolution at start of tracker;MeV",100,-2.5,2.5);
  TH2F* mom = new TH2F("mom","momentum at start of tracker;true momentum (MeV);fit momentum (MeV)",
    100,90,107,100,90,107);
  mom->SetStats(0);
  trk->Project("fstat","fitstatus","mcentmom>100");
  trk->Project("nh","nhits","kalfail==0");
  trk->Project("na","nactive","kalfail==0");
  trk->Project("nd","nhits-nactive","kalfail==0");
//  trk->Project("momr","fitmom-mcmom","kalfail==0");
//  trk->Project("momr","fitmom-mcmom","kalfail==0&&nactive>=20&&fitmom>100&&t0err<1&&chisq/ndof<5&&fitmomerr<0.2");
  trk->Project("momr","fitmom-mcentmom","mcentmom>100&&kalfail==0");
  //"&&t0err<0.8&&fitmomerr<0.1&&chisq/ndof<2");
  kcan->Clear();
  kcan->Divide(2,2);
  kcan->cd(1);
  nd->Draw();
  na->Draw("same");
  nh->Draw("same");
  TLegend* leg = new TLegend(0.3,0.7,0.8,0.9);
  leg->AddEntry(nh,"All hits","l");
  leg->AddEntry(na,"Active hits","l");
  leg->AddEntry(nd,"Disabled hits","l");
  leg->Draw();
  kcan->cd(2);
  t0->SetStats(1);
  t0->Fit("gaus");
  kcan->cd(3);
  trk->Draw("fitmom:mcentmom>>mom","kalfail==0&&nactive>=20");
  kcan->cd(4);
  momr->SetStats(1);
  double integral = momr->GetEntries()*momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,momr->GetRMS(),momr->GetRMS(),0.01,2*momr->GetRMS(),2*momr->GetRMS());
  sgau->SetParLimits(5,1.0*momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.8);
  momr->Fit("sgau","L");
  
  double keff = 0.5*momr->GetEntries()/fstat->GetEntries();
  TPaveText* text = new TPaveText(0.1,0.7,0.4,0.8,"NDC");  
  char line[40];
  sprintf(line,"P>100MeV/c Eff. = %3.2f",keff);
  text->AddText(line);
  text->Draw();

  }  
  void HelixTest (TTree* trk) {
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111111);

  TCanvas* ccan = new TCanvas("ccan","Circle Fit",1200,800);
  TH2F* cx = new TH2F("cx","Circle X center;MC true (mm);Fit (mm)",50,-400,400,50,-400,400);
  TH2F* cy = new TH2F("cy","Circle Y center;MC true (mm);Fit (mm)",50,-400,400,50,-400,400);
  TH2F* cr = new TH2F("cr","Circle radius;MC true (mm);Fit (mm)",50,180,310,50,180,310);
  TH1F* cxr = new TH1F("cxr","X center resolution;mm",100,-200,200);
  TH1F* cyr = new TH1F("cyr","Y center resolution;mm",100,-200,200);
  TH1F* crr = new TH1F("crr","Radius resolution;mm",100,-200,200);
  cx->SetStats(0);
  cy->SetStats(0);
  cr->SetStats(0);
  TCut helix("helixfail==0");

  trk->Project("cxr","hcx-mccx",helix);
  trk->Project("cyr","hcy-mccy",helix);
  trk->Project("crr","hr-mcr",helix);
  ccan->Clear();
  ccan->Divide(3,2);
  ccan->cd(1);
  trk->Draw("hcx:mccx>>cx",helix);
  ccan->cd(2);
  trk->Draw("hcy:mccy>>cy",helix);
  ccan->cd(3);
  trk->Draw("hr:mcr>>cr",helix);
  ccan->cd(4);
  cxr->Fit("gaus");
  ccan->cd(5);
  cyr->Fit("gaus");
  ccan->cd(6);
  crr->Fit("gaus");
  
  TCanvas* hcan = new TCanvas("hcan","Helix Fit",1200,800);
  
  TH2F* dfdz = new TH2F("dfdz","Helix pitch (d#phi/dZ);MC true (radians/mm);helix fit (radians/mm)",50,0.0035,0.0065,50,0.0035,0.0065);
  TH2F* fz0 = new TH2F("fz0","Heliz #phi intercept;MC true (radians); helix fit (radians)",50,-1,12,50,-1,12.);
  TH1F* dfdzr = new TH1F("dfdzr","d#phi/dZ resolution;radians/mm",100,-0.0005,0.0005);
  TH1F* fz0r = new TH1F("fz0r","#phi intercept resolution;radians",100,-0.4,0.4);
  trk->Project("dfdzr","hdfdz-mcdfdz",helix);
  trk->Project("fz0r","hfz0-mcfz0",helix);
  dfdz->SetStats(0);
  fz0->SetStats(0);
  
  hcan->Clear();
  hcan->Divide(2,2);
  hcan->cd(1);
  trk->Draw("hdfdz:mcdfdz>>dfdz",helix);
  hcan->cd(2);
  trk->Draw("hfz0:mcfz0>>fz0",helix);
  hcan->cd(3);
  dfdzr->Fit("gaus");
  hcan->cd(4);
  fz0r->Fit("gaus");
  
  TCanvas* scan = new TCanvas("scan","Helix Seed",1200,800);
  
  TH2F* d0 = new TH2F("d0","d0;MC true (mm);helix fit (mm)",50,-200,200,50,-200,200);
  TH2F* phi0 = new TH2F("phi0","#phi0;MC true (radians);helix fit (radians)",50,-3.5,3.5,50,-3.5,3.5);
  TH2F* om = new TH2F("om","#omega;MC true (1/mm);helix fit (1/mm)",50,0.002,0.005,50,0.002,0.005);
  TH2F* z0 = new TH2F("z0","z0;MC true (mm);helix fit (mm)",50,-1000,1000,50,-1000,1000);
  TH2F* td = new TH2F("td","tan(#lambda);MC true;helix fit",50,0.5,1.2,50,0.5,1.2);
  d0->SetStats(0);
  phi0->SetStats(0);
  om->SetStats(0);
  z0->SetStats(0);
  td->SetStats(0);

  
  scan->Clear();
  scan->Divide(3,2);
  scan->cd(1);
  trk->Draw("hd0:mcmidd0>>d0");
  scan->cd(2);
  trk->Draw("hp0:mcmidp0>>phi0");
  scan->cd(3);
  trk->Draw("hom:mcmidom>>om");
  scan->cd(4);
  trk->Draw("hz0:mcmidz0>>z0");
  scan->cd(5);
  trk->Draw("htd:mcmidtd>>td");
}

void MomRes(TTree* trk) {
  TCanvas* mcan = new TCanvas("mcan","Momentum",1200,800);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
// should have pitch angle and generated hit cuts here, FIXME!!!
  TCut mcsel("mcentmom>100&&mcenttd<1.0&&mcenttd>0.5774&&nchit>=20");
  TCut tsel = mcsel +TCut("kalfail==0&&nhit>=25");
// selection cuts
  TCut cuts[4];

  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trk->Project("effnorm","mcentmom",mcsel);
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TH1F* momres[4];
// cuts for different tightness of selection
  TCut cuts[4];
  cuts[0] = tsel;
  cuts[1] = tsel+"t0err<2.8&&fitmomerr<0.22&&chisq/ndof<5.0";
  cuts[2] = tsel+"t0err<1.4&&fitmomerr<0.15&&chisq/ndof<3.0";
  cuts[3] = tsel+"t0err<1.0&&fitmomerr<0.1&&chisq/ndof<2.0";

  const char* labels[4] = {"A","B","C","D"};

  mcan->Clear();
  mcan->Divide(2,2);
  for(unsigned ires=0;ires<4;ires++){
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",151,-2.5,2.5);
    trk->Project(mname,"fitmom-mcentmom",cuts[ires]);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momres[ires]->GetRMS(),momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
    sgau->SetParLimits(5,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momres[ires]->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.8);
    mcan->cd(ires+1);
    momres[ires]->Fit("sgau","L");

    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
    TPaveText* text = new TPaveText(0.1,0.7,0.4,0.8,"NDC");  
    char line[40];
    sprintf(line,"%s   Efficiency = %3.2f",labels[ires],keff);
    text->AddText(line);
    text->Draw();


  }  
  
}
  
void SeedTest (TTree* trk) {
  
  TCanvas* srcan = new TCanvas("srcan","Helix Seed Resolution",1200,800);
  
  TH1F* d0r = new TH1F("d0r","seed d0 resolution;mm",100,-300,300);
  TH1F* phi0r = new TH1F("phi0r","seed #phi0 resolution;radians",100,-0.25,0.25);
  TH1F* omr = new TH1F("omr","seed #omega resolution;1/mm",100,-0.002,0.002);
  TH1F* z0r = new TH1F("z0r","seed z0 resolution;mm",100,-150,150);
  TH1F* tdr = new TH1F("tdr","seed tan(#lambda) resolution",100,-0.3,0.3);


  trk->Project("d0r","hd0-mcmidd0");
  trk->Project("phi0r","hp0-mcmidp0");
  trk->Project("omr","hom-mcmidom");
  trk->Project("z0r","hz0-mcmidz0");
  trk->Project("tdr","htd-mcmidtd");
  
  srcan->Clear();
  srcan->Divide(3,2);
  srcan->cd(1);
  d0r->Fit("gaus");
  srcan->cd(2);
  phi0r->Fit("gaus");
  srcan->cd(3);
  omr->Fit("gaus");
  srcan->cd(4);
  z0r->Fit("gaus");
  srcan->cd(5);
  tdr->Fit("gaus");
  
  
}
