void KalFitHit (TTree* hits ) {
  TCut goodhit("active");
    TCanvas* dcan = new TCanvas("driftcan","driftcan",1200,800);
    TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
    TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-8,8);
    TH2F* drad = new TH2F("drad","Drift radius;true drift radius (mm);reco drift radius (mm)",
      100,-0.2,2.7,100,-0.2,2.7);
    drad->SetStats(0);
    TH1F* rpull = new TH1F("rpull","residual pull",100,-8,8);
    TCut gooddrift("0.5*(rdrift+mcrdrift)>0.2&&0.5*(rdrift+mcrdrift)<2.3");
    hits->Project("dres","rdrift-mcrdrift",goodhit+gooddrift);
    hits->Project("dpull","(rdrift-mcrdrift)/rdrifterr",goodhit+gooddrift);
    hits->Project("rpull","resid/residerr",goodhit+gooddrift);
    dcan->Clear();
    dcan->Divide(2,2);
    dcan->cd(1);
    hits->Draw("rdrift:mcrdrift>>drad",goodhit);
    dcan->cd(2);
    dres->Fit("gaus");
    dcan->cd(3);
    dpull->Fit("gaus");
    dcan->cd(4);
    rpull->Fit("gaus");
    
    TCanvas* tcan = new TCanvas("hitt0can","hitt0can",1200,800);
    TH1F* t0res = new TH1F("t0res","hit t0 resolution;nsec",100,-5,5);
    TH1F* t0pull = new TH1F("t0pull","hit t0 pull",100,-8,8);
    TH2F* dt0 = new TH2F("dt0","Hit t0;true t0 (nsec);reco t0 (nsec)",
        100,500,4000,100,500,4000);
    hits->Project("t0res","hitt0-mchitt0",goodhit);
    hits->Project("t0pull","(hitt0-mchitt0)/hitt0err",goodhit);
    tcan->Clear();
    tcan->Divide(2,2);
    tcan->cd(1);
    hits->Draw("hitt0:mchitt0>>dt0",goodhit);
    tcan->cd(2);
    t0res->Fit("gaus");
    tcan->cd(3);
    t0pull->Fit("gaus");
    
    
    TCanvas* tdcan = new TCanvas("tdcan","tdcan",1200,800);
    TH1F* tdres = new TH1F("tdres","#Deltat V resolution;mm",100,-200,200);
    TH1F* tdpull = new TH1F("tdpull","#Deltat V pull",100,-8,8);
    TH2F* dtd = new TH2F("dtd","Hit V position;true V (mm);#Deltat V (mm)",
      100,-600,600,100,-600,600);
    TH2F* pocatd = new TH2F("pocatd","Hit POCA V;true V (mm);POCA V (mm)",
      100,-600,600,100,-600,600);
        
    hits->Project("tdres","dmid-mcdmid",goodhit);
    hits->Project("tdpull","(dmid-mcdmid)/dmiderr",goodhit);
    tdcan->Clear();
    tdcan->Divide(2,2);
    tdcan->cd(1);
    hits->Draw("dmid:mcdmid>>dtd",goodhit);
    tdcan->cd(2);
    hits->Draw("hflt:mcdmid>>pocatd",goodhit);
    tdcan->cd(3);
    tdres->Draw();
    tdcan->cd(4);
    tdpull->Fit("gaus");
    
}

void KalFitTrk (TTree* trks ) {
  
  TCanvas* ttcan = new TCanvas("trkt0can","trkt0can",1200,800);
  TH1F* t00res = new TH1F("t00res","Initial t0 resolution;nsec",100,-20,20);
  TH1F* t0res = new TH1F("t0res","Final t0 resolution;nsec",100,-5,5);
  TH1F* t0pull = new TH1F("t0pull","Track t0 pull",100,-8,8);
  TH2F* dt0 = new TH2F("dt0","Track t0;true t0 (nsec);Initial t0 (nsec)",
    100,500,4000,100,500,4000);
  trks->Project("t00res","t00-mct0","fitstatus>0");
  trks->Project("t0res","t0-mct0","fitstatus>0");
  trks->Project("t0pull","(t0-mct0)/t0err","fitstatus>0");
  ttcan->Clear();
  ttcan->Divide(2,2);
  ttcan->cd(1);
  trks->Draw("t00:mct0>>dt0","fitstatus>0");
  ttcan->cd(2);
  t00res->Draw();
  ttcan->cd(3);
  t0res->Fit("gaus");
  ttcan->cd(4);
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
