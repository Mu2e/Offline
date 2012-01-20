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

KalFitAcc(TTree* trks) {
  TCanvas* acan = new TCanvas("acan","Acceptance",1200,800);
  TH1F* pcost = new TH1F("pcost","e^{-} cos(#theta) at production",105,-1.02,1.02);
  TH1F* pcosta = new TH1F("pcosta","e^{-} cos(#theta) at production",105,-1.02,1.02);
  pcosta->SetLineColor(kBlue);

  TCut hittrk("mcentmom>0.0");
  TCut pitch("mcenttd/sqrt(1+mcenttd^2)>0.5&&mcenttd/sqrt(1+mcenttd^2)<0.707107");
  TCut mom("mcentmom>102");
  TCut reco("fitstatus>0");
  TCut goodfit("fitcon>1e-2&&nactive>=25");

  trks->Project("pcost","mccost");
  trks->Project("pcosta","mccost",hittrk);
//  pcosta->Divide(pcost);

  TH1F* dcost = new TH1F("dcost","e^{-} cos(#theta) at tracker",105,0,1.02);
  TH1F* dmom = new TH1F("dmom","e^{-} momentum at tracker",100,90,108);
  TH1F* rmom = new TH1F("rmom","recsonstructed e^{-} momentum",100,90,108);
  TH1F* fitcon = new TH1F("fitcon","fit consistency",205,-0.01,1.01);
  TH1F* acc = new TH1F("acc","Acceptance",6,-0.5,5.5);
  acc->GetXaxis()->SetBinLabel(1,"All");
  acc->GetXaxis()->SetBinLabel(2,"Reaches Tracker");
  acc->GetXaxis()->SetBinLabel(3,"Pitch");
  acc->GetXaxis()->SetBinLabel(4,"Momentum");
  acc->GetXaxis()->SetBinLabel(5,"Reconstructed");
  acc->GetXaxis()->SetBinLabel(6,"Fit Quality");

  trks->Project("dcost","mcenttd/sqrt(1+mcenttd^2)",hittrk);
  trks->Project("dmom","mcentmom",hittrk+pitch);
  trks->Project("fitcon","fitcon",hittrk+pitch+mom+reco);
  trks->Project("rmom","fitmom",hittrk+pitch+mom+reco+goodfit);
//  trks->Project("+fitcon","-0.05",hittrk+pitch+mom+"fitstatus<0");

  trks->Project("acc","0.0");
  trks->Project("+acc","1.0",hittrk);
  trks->Project("+acc","2.0",hittrk+pitch);
  trks->Project("+acc","3.0",hittrk+pitch+mom);
  trks->Project("+acc","4.0",hittrk+pitch+mom+reco);
  trks->Project("+acc","5.0",hittrk+pitch+mom+reco+goodfit);

  acan->Clear();
  acan->Divide(3,2);
  acan->cd(1);
  pcost->Draw();
  pcosta->Draw("same");
  acan->cd(2);
  dcost->Draw();
  TLine* costcut_h = new TLine(0.707107,0.0,0.707107,dcost->GetMaximum());
  costcut_h->SetLineColor(kBlack);
  costcut_h->SetLineStyle(2);
  costcut_h->SetLineWidth(2);
  TLine* costcut_l = new TLine(0.5,0.0,0.5,dcost->GetMaximum());
  costcut_l->SetLineColor(kBlack);
  costcut_l->SetLineStyle(2);
  costcut_l->SetLineWidth(2);
  costcut_h->Draw();
  costcut_l->Draw();
  acan->cd(3);
  dmom->Draw();
  TLine* momcut_l = new TLine(102,0.0,102,dmom->GetMaximum());
  momcut_l->SetLineColor(kBlack);
  momcut_l->SetLineStyle(2);
  momcut_l->SetLineWidth(2);
  momcut_l->Draw();
  
  acan->cd(4);
  fitcon->Draw();
  TLine* concut_l = new TLine(1e-2,0.0,1e-2,fitcon->GetMaximum());
  concut_l->SetLineColor(kBlack);
  concut_l->SetLineStyle(2);
  concut_l->SetLineWidth(2);
  concut_l->Draw();
 
  acan->cd(5);
  rmom->Draw();

  acan->cd(6);
  acc->Draw();

  double all = acc->GetBinContent(1);
  double reach = acc->GetBinContent(2)/all;
  double gpitch = acc->GetBinContent(3)/all;
  double gmom = acc->GetBinContent(4)/all;
  double greco = acc->GetBinContent(5)/all;
  double good = acc->GetBinContent(6)/all;


  cout << "Acceptance: Reach tracker " << reach
  << " Pitch " << gpitch
  << " Momentum " << gmom
  << " Reconstruction " << greco
  << " Fit Quality " << good << std::endl;



}
