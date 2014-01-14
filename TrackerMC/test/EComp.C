{
  TH2F* ecomp = new TH2F("ecomp","Electron Hit Energy Comparison;MC Edep (MeV);StrawHit Edep (MeV)",100,0,0.006,100,0,0.006);
  se->Project("ecomp","edep:mcedep");
  ecomp->SetStats(0);
  TCanvas* ec = new TCanvas("ec","Energy comparison",1000,600);
  ec->Divide(2,2);
  ec->cd(1);
  ecomp->Draw("box");
  TH1F* eedep = new TH1F("eedep","Straw Hit Energy;MeV",100,0,0.04);
  TH1F* pedep = new TH1F("pedep","Straw Hit Energy;MeV",100,0,0.04);
  se->Project("eedep","edep");
  sp->Project("pedep","edep");
  eedep->SetLineColor(kRed);
  pedep->SetLineColor(kBlue);
  eedep->SetStats(0);
  pedep->SetStats(0);
  ec->cd(2);
  eedep->Draw();
  pedep->Draw("same");
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(eedep,"Electrons","l");
  leg->AddEntry(pedep,"Protons","l");
  leg->Draw();

  TH1F* emcedep = new TH1F("emcedep","MC True EDep;MeV",100,0,0.04);
  TH1F* pmcedep = new TH1F("pmcedep","MC True EDep;MeV",100,0,0.04);
  se->Project("emcedep","mcedep");
  sp->Project("pmcedep","mcedep");
  emcedep->SetLineColor(kRed);
  pmcedep->SetLineColor(kBlue);
  emcedep->SetStats(0);
  pmcedep->SetStats(0);
  ec->cd(3);
  emcedep->Draw();
  pmcedep->Draw("same");
  leg->Draw();

}
