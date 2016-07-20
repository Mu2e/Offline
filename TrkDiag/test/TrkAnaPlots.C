void Draw(TTree* c) {

  TH2F* evspep = new TH2F("evspep","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmp = new TH2F("evspmp","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspem = new TH2F("evspem","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmm = new TH2F("evspmm","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  evspep->SetMaximum(50);
  evspmp->SetMaximum(50);
  evspep->SetStats(0);
  evspmp->SetStats(0);
  evspep->SetLineColor(kOrange);
  evspmp->SetLineColor(kCyan);
  evspem->SetMaximum(50);
  evspmm->SetMaximum(50);
  evspem->SetStats(0);
  evspmm->SetStats(0);
  evspem->SetLineColor(kRed);
  evspmm->SetLineColor(kBlue);

  TCanvas* can = new TCanvas("can","can",800,800);

  c->Project("evspem","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==11");
  c->Project("evspmm","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==13");
  c->Project("evspep","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==-11");
  c->Project("evspmp","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==-13");

  evspem->Draw();
  evspem->Draw("box");
  evspmm->Draw("boxsame");
  TLegend* leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry(evspem,"True e^{-}","L");
  leg->AddEntry(evspmm,"True #mu^{-}","L");
  leg->AddEntry(evspep,"True e^{+}","L");
  leg->AddEntry(evspmp,"True #mu^{+}","L");
  leg->Draw();
}
