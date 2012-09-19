void TestFlash(TTree* ntvd) {
  TH2F* xy = new TH2F("xy","Y vs X position at TS5 exit",100,-4200,-3600,100,-300,300);
  TH1F* pdg = new TH1F("pdg","PDG code",200,-250,2250);
  TH1F* emt = new TH1F("emt","Time at TS5",100,0,400);
  TH1F* ept = new TH1F("ept","Time at TS5",100,0,400);
  TH1F* pt = new TH1F("pt","Time at TS5",100,0,400);
  TH1F* mut = new TH1F("mut","Time at TS5",100,0,400);
  TH1F* pimt = new TH1F("pimt","Time at TS5",100,0,400);
  TH1F* pht = new TH1F("pht","Time at TS5",100,0,400);
  TH1F* nt = new TH1F("nt","Time at TS5",100,0,400);
  emt->SetLineColor(kRed);
  ept->SetLineColor(kBlack);
  pt->SetLineColor(kOrange);
  mut->SetLineColor(kBlue);
  pimt->SetLineColor(kGreen);
  pht->SetLineColor(kMagenta);
  nt->SetLineColor(kCyan);
  emt->SetStats(0);

  TH1F* emke = new TH1F("emke","Kinetic Energy at TS5",100,0,100);
  TH1F* epke = new TH1F("epke","Kinetic Energy at TS5",100,0,100);
  TH1F* pke = new TH1F("pke","Kinetic Energy at TS5",100,0,100);
  TH1F* muke = new TH1F("muke","Kinetic Energy at TS5",100,0,100);
  TH1F* pimke = new TH1F("pimke","Kinetic Energy at TS5",100,0,100);
  TH1F* phke = new TH1F("phke","Kinetic Energy at TS5",100,0,100);
  TH1F* nke = new TH1F("nke","Kinetic Energy at TS5",100,0,100);

  emke->SetLineColor(kRed);
  epke->SetLineColor(kBlack);
  pke->SetLineColor(kOrange);
  muke->SetLineColor(kBlue);
  pimke->SetLineColor(kGreen);
  phke->SetLineColor(kMagenta);
  nke->SetLineColor(kCyan);
  emke->SetStats(0);

  ntvd->Project("xy","y:x");
  ntvd->Project("pdg","pdg");

  ntvd->Project("emt","time","pdg==11");
  ntvd->Project("ept","time","pdg==-11");
  ntvd->Project("pt","time","pdg==2212");
  ntvd->Project("mut","time","pdg==13");
  ntvd->Project("pimt","time","pdg==-211");
  ntvd->Project("pht","time","pdg==22");
  ntvd->Project("nt","time","pdg==2112");
  emt->SetMinimum(1);

  ntvd->Project("emke","ke","pdg==11");
  ntvd->Project("epke","ke","pdg==-11");
  ntvd->Project("pke","ke","pdg==2212");
  ntvd->Project("muke","ke","pdg==13");
  ntvd->Project("pimke","ke","pdg==-211");
  ntvd->Project("phke","ke","pdg==22");
  ntvd->Project("nke","ke","pdg==2112");
  emke->SetMinimum(1);

  TCanvas* flashcan = new TCanvas("flashcan","flashcan",1200,800);
  flashcan->Divide(2,2);
  flashcan->cd(1);
  xy->Draw("box");
  flashcan->cd(2);
  gPad->SetLogy();
  pdg->Draw();
  flashcan->cd(3);
  gPad->SetLogy();
  emt->Draw();
  ept->Draw("same");
  mut->Draw("same");
  pt->Draw("same");
  pimt->Draw("same");
  pht->Draw("same");
  nt->Draw("same");
  TLegend* leg = new TLegend(0.6,0.4,0.8,0.9);
  leg->AddEntry(emke,"e^{-}","l");
  leg->AddEntry(epke,"e^{+}","l");
  leg->AddEntry(muke,"#mu^{-}","l");
  leg->AddEntry(pke,"P^{+}","l");
  leg->AddEntry(pimke,"#pi^{-}","l");
  leg->AddEntry(phke,"#gamma","l");
  leg->AddEntry(nke,"neutron","l");
  leg->Draw();

  flashcan->cd(4);
  gPad->SetLogy();
  emke->Draw();
  epke->Draw("same");
  muke->Draw("same");
  pke->Draw("same");
  pimke->Draw("same");
  phke->Draw("same");
  nke->Draw("same");


}


