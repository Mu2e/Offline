void CDRHits(TTree* hits) {
  gStyle->SetOptStat(0);
  TCanvas* ecan = new TCanvas("ecan","energy",1200,800);
  TCanvas* rcan = new TCanvas("rcan","radius",1200,800);
    
  TH1F* econv = new TH1F("econv","Straw Hit Energy;MeV",200,-0.01,0.15);
  TH1F* edio = new TH1F("edio","Straw Hit Energy;MeV",200,-0.01,0.15);
  TH1F* edelta = new TH1F("edelta","Straw Hit Energy;MeV",200,-0.01,0.15);
  TH1F* ep = new TH1F("ep","Straw Hit Energy;MeV",200,-0.01,0.15);
  econv->SetLineColor(kRed);
  edio->SetLineColor(kGreen);
  edelta->SetLineColor(kCyan);
  ep->SetLineColor(kBlue);
  econv->SetStats(0);
  edio->SetStats(0);
  edelta->SetStats(0);
  ep->SetStats(0);

  TH1F* rconv = new TH1F("rconv","StrawHit Radius;mm",100,360,780);
  TH1F* rdio = new TH1F("rdio","StrawHit Radius;mm",100,360,780);
  TH1F* rdelta = new TH1F("rdelta","StrawHit Radius;mm",100,360,780);
  TH1F* rp = new TH1F("rp","StrawHit Radius;mm",100,360,780);
  rconv->SetLineColor(kRed);
  rdio->SetLineColor(kGreen);
  rdelta->SetLineColor(kCyan);
  rp->SetLineColor(kBlue);

  TCut ecut("edep<0.0045");
  TCut rmin("sqrt(shpos.x^2+shpos.y^2)>410");

  hits->Project("econv","edep","mcpdg==11&&mcgen==2");
  hits->Project("edio","edep","mcpdg==11&&mcgen==6");
  hits->Project("edelta","edep","mcpdg==11&&mcgen<0");
  hits->Project("ep","edep","mcpdg==2212");

  hits->Project("rconv","sqrt(shpos.y^2+shpos.x^2)","mcgen==2");
  hits->Project("rdio","sqrt(shpos.y^2+shpos.x^2)","mcgen==6");
  hits->Project("rdelta","sqrt(shpos.y^2+shpos.x^2)","mcgen<0");
  hits->Project("rp","sqrt(shpos.y^2+shpos.x^2)","mcpdg==2212");
    
  TLegend* leg2 = new TLegend(0.55,0.65,0.9,0.9);
  leg2->AddEntry(rconv,"Conv. Electrons","l");
  leg2->AddEntry(rdio,"DIO Electrons","l");
  leg2->AddEntry(rdelta,"Delta Electrons","l");
  leg2->AddEntry(rp,"Protons","l");

  ecan->Clear();
  ecan->Divide(1,1);
  ecan->cd(1);
  gPad->SetLogy();
  edelta->Draw();
  ep->Draw("same");
  econv->Draw("same");
  edio->Draw("same");
  leg2->Draw();


//  TLine* ecut_t = new TLine(0.004,0.0,0.004,econv->GetMaximum());
//  ecut_t->SetLineColor(kBlack);
//  ecut_t->SetLineStyle(2);
//  ecut_t->SetLineWidth(2);
//  TLine* ecut_l = new TLine(0.0055,0.0,0.0055,econv->GetMaximum());
//  ecut_l->SetLineColor(kBlack);
//  ecut_l->SetLineStyle(3);
//  ecut_l->SetLineWidth(2);
//  ecut_t->Draw();
//  ecut_l->Draw();

//  TLegend* leg3 = new TLegend(0.55,0.7,0.9,0.9);
//  leg3->AddEntry(ecut_t,"Tight cut","l");
//  leg3->AddEntry(ecut_l,"Loose cut","l");
//  leg3->Draw();

  rcan->Clear();
  rcan->Divide(1,1);
  rcan->cd(1);

  rdelta->Draw();
  rp->Draw("same");
  rconv->Draw("same");
  rdio->Draw("same");
  
//  TLine* rmin_t = new TLine(410,0.0,410,rp->GetMaximum());
//  rmin_t->SetLineColor(kBlack);
//  rmin_t->SetLineStyle(2);
//  rmin_t->SetLineWidth(2);
//  TLine* rmin_l = new TLine(390,0.0,390,rp->GetMaximum());
//  rmin_l->SetLineColor(kBlack);
//  rmin_l->SetLineStyle(3);
//  rmin_l->SetLineWidth(2);
//  
//  TLine* rmax_t = new TLine(630,0.0,630,rp->GetMaximum());
//  rmax_t->SetLineColor(kBlack);
//  rmax_t->SetLineStyle(2);
//  rmax_t->SetLineWidth(2);
//  TLine* rmax_l = new TLine(650,0.0,650,rp->GetMaximum());
//  rmax_l->SetLineColor(kBlack);
//  rmax_l->SetLineStyle(3);
//  rmax_l->SetLineWidth(2);
//  rmin_t->Draw();
//  rmin_l->Draw();
//  rmax_t->Draw();
//  rmax_l->Draw();
  
  leg2->Draw();
  
  TCanvas* ccan = new TCanvas("ccan","cleaned hits",1200,800);
  TCut clean("loose>0&&abs(time-tpeaks[0]s[0])<40.0");
  
  TH1F* tconv = new TH1F("tconv","Hit Time WRT Peak;nsec",101,-80,80);
  TH1F* tdio = new TH1F("tdio","Hit Time WRT Peak;nsec",101,-80,80);
  TH1F* tdelta = new TH1F("tdelta","Hit Time WRT Peak;nsec",101,-80,80);
  TH1F* tp = new TH1F("tp","Hit Time WRT Peak;nsec",101,-80,80);
  tconv->SetLineColor(kRed);
  tdio->SetLineColor(kGreen);
  tdelta->SetLineColor(kCyan);
  tp->SetLineColor(kBlue);
  tconv->SetStats(0);

  hits->Project("tconv","time-tpeaks[0]","mcpdg==11&&mcgen==2&&loose");
  hits->Project("tdio","time-tpeaks[0]","mcpdg==11&&mcgen==6&&loose");
  hits->Project("tdelta","time-tpeaks[0]","mcpdg==11&&mcgen<0&&loose");
  hits->Project("tp","time-tpeaks[0]","mcpdg==2212&&loose");
    
  TLine* tmin = new TLine(-40,0.0,-40,tconv->GetMaximum());
  tmin->SetLineColor(kBlack);
  tmin->SetLineStyle(2);
  tmin->SetLineWidth(2);
  TLine* tmax = new TLine(40,0.0,40,tconv->GetMaximum());
  tmax->SetLineColor(kBlack);
  tmax->SetLineStyle(2);
  tmax->SetLineWidth(2);
  
  ccan->Clear();
  ccan->Divide(1,1);
  ccan->cd(1);
  tconv->Draw();
  tdio->Draw("same");
  tdelta->Draw("same");
  tp->Draw("same");
  leg2->Draw();
  tmin->Draw();
  tmax->Draw();
}

void
CDROutliers(TTree* rhits) {
  TH1F* rconv = new TH1F("rconv","DOCA to robust helix",100,-0.1,400);
  TH1F* rp = new TH1F("rp","DOCA to robust helix",100,-0.1,400);
  TH1F* rdio = new TH1F("rdio","DOCA to robust helix",100,-0.1,400);
  TH1F* rdelta = new TH1F("rdelta","DOCA to robust helix",100,-0.1,400);
  rconv->SetLineColor(kRed);
  rdio->SetLineColor(kGreen);
  rdelta->SetLineColor(kCyan);
  rp->SetLineColor(kBlue);
  rconv->SetStats(0);

  rhits->Project("rconv","resid","mcpdg==11&&mcgen==2");
  rhits->Project("rp","resid","mcpdg==2212");
  rhits->Project("rdio","resid","mcpdg==11&&mcgen==6");
  rhits->Project("rdelta","resid","mcpdg==11&&mcgen<0");

  TCanvas* rcan = new TCanvas("rcan","robust helix DOCA",1200,800);
  rcan->Divide(1,1);
  rcan->cd(1);
  gPad->SetLogy();
  rconv->Draw();
  rdio->Draw("same");
  rdelta->Draw("same");
  rp->Draw("same");
  TLine* tmin = new TLine(160,0.0,160,rconv->GetMaximum());
  tmin->SetLineColor(kBlack);
  tmin->SetLineStyle(2);
  tmin->SetLineWidth(2);
  tmin->Draw();

  TLegend* leg2 = new TLegend(0.55,0.65,0.9,0.9);
  leg2->AddEntry(rconv,"Conv. Electrons","l");
  leg2->AddEntry(rdio,"DIO Electrons","l");
  leg2->AddEntry(rdelta,"Delta Electrons","l");
  leg2->AddEntry(rp,"Protons","l");

  leg2->Draw();

  Int_t ilow = rconv->FindBin(0.0);
  Int_t ihi = rconv->FindBin(160.0);
  Double_t convint = rconv->Integral(ilow,ihi);
  Double_t pint = rp->Integral(ilow,ihi);
  Double_t dioint = rdio->Integral(ilow,ihi);
  Double_t deltaint = rdelta->Integral(ilow,ihi);
  double total = convint + dioint + deltaint + pint;
  double eff = convint/rconv->GetEntries();
  double pfrac = pint/total;
  double diofrac = dioint/total;
  double deltafrac = deltaint/total;
  double purity = convint/total;
  cout << "eff = " << eff << " purity = " << purity << " pfrac = " << pfrac << " diofrac = " << diofrac << " deltafrac = " << deltafrac << endl;

}

void
CDRT0 (TTree* trkd) {

  TH1F* t00res = new TH1F("t00res","Initial T_{0} resolution;nsec",100,-10,10);
  TH1F* t00pull = new TH1F("t00pull","Initial T_{0} pull",100,-20,20);
  TH1F* t0res = new TH1F("t0res","Final T_{0} resolution;nsec",100,-5,5);
  TH1F* t0pull = new TH1F("t0pull","Final T_{0} pull",100,-10,10);
  trkd->Project("t00res","t00-mcmidt0","kalfail==0");
  trkd->Project("t00pull","(t00-mcmidt0)/t00err","kalfail==0");
  trkd->Project("t0res","t0-mcmidt0","kalfail==0");
  trkd->Project("t0pull","(t0-mcmidt0)/t0err","kalfail==0");


  TCanvas* t0can = new TCanvas("t0can","T_{0} canvas",1200,800);
  t0can->Clear();
  t0can->Divide(2,2);
  t0can->cd(1);
  t00res->Fit("gaus");
  t0can->cd(2);
  t00pull->Draw();
  t0can->cd(3);
  t0res->Fit("gaus");
  t0can->cd(4);
  t0pull->Fit("gaus");

}
