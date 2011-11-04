void PlotTimeSpectra(TDirectory* tdir,unsigned nmax=20, unsigned nps=3){
  gStyle->SetOptStat(0);
  bool moreplots(true);
  unsigned ican(0);
  TCanvas* cans[50];
  unsigned nplots(0);
  TH1F* dummy = new TH1F("dummy","dummy",10,0,1);
  dummy->SetMarkerStyle(23);
  dummy->SetMarkerColor(kRed);
//  dummy->SetMarkerSize(1);
  while(moreplots && nplots < nmax) {
    bool first(true);
    char cname[20];
    snprintf(cname,20,"timecan%i",ican);
    cans[ican] = new TCanvas(cname,"Straw Hit Times",800,800);
    cans[ican]->Clear();
    cans[ican]->Divide(nps,nps);
    char rname[100];
    char tname[100];
    char lname[100];
    char cname[100];
    for(unsigned iplot=0;iplot<nps*nps;++iplot){
      nplots++;
      unsigned jplot = ican*nps*nps+1+iplot;
      snprintf(rname,100,"rawtspectrum%i",jplot);
      snprintf(tname,100,"tighttspectrum%i",jplot);
      snprintf(lname,100,"loosetspectrum%i",jplot);
      snprintf(cname,100,"convtspectrum%i",jplot);
      TH1F* rh = tdir->Get(rname);
      TH1F* th = tdir->Get(tname);
      TH1F* lh = tdir->Get(lname);
      TH1F* ch = tdir->Get(cname);
      if(rh != 0 && th != 0 &&lh != 0 && ch != 0){
        rh->SetStats(0);
	th->SetStats(0);
        lh->SetStats(0);
        ch->SetStats(0);
        cans[ican]->cd(iplot+1);
        rh->SetTitle("Time Spectrum");
        rh->GetXaxis()->SetTitle("nsec");
        rh->GetYaxis()->SetTitle("# hits");
        rh->Draw();
//        th->SetLineWidth(2);
	lh->SetFillColor(kYellow);
        lh->Draw("same");
	th->SetFillColor(kGreen);
        th->Draw("same");
	ch->SetLineStyle(1);
        ch->SetLineWidth(2);
        ch->Draw("same");
        if(first){
          first = false;
          TLegend* leg = new TLegend(0.4,0.7,0.9,0.9);
          leg->AddEntry(rh,"All hits","l");
	  leg->AddEntry(th,"Tight Selected hits","F");
          leg->AddEntry(lh,"Loose Selected hits","F");
          leg->AddEntry(ch,"Conversion hits","L");
          leg->AddEntry(dummy,"Reconstructed Time Peak","P");
          leg->Draw();
        }
      } else {
        moreplots = false;
        break;
      }
    }
    char fname[50];
    snprintf(fname,20,"shtime%i.png",ican);
    cans[ican]->SaveAs(fname);
    ican++;
  }
}
