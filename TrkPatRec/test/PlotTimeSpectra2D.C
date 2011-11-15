void PlotTimeSpectra2D(TDirectory* tdir,double sigma=2.0,double minn=5,unsigned nmax=20, unsigned nps=2){
  gStyle->SetOptStat(0);
  bool moreplots(true);
  unsigned ican(0);
  TCanvas* cans[50];
  TSpectrum2 tp2(100);
  unsigned nplots(0);
  TH1F* dummy = new TH1F("dummy","dummy",10,0,1);
  dummy->SetMarkerStyle(23);
  dummy->SetMarkerColor(kRed);
  dummy->SetMarkerSize(1);
  while(moreplots && nplots < nmax) {
    bool first(true);
    char cname[20];
    snprintf(cname,20,"timecan%i",ican);
    cans[ican] = new TCanvas(cname,"Straw Hit Times",800,800);
    cans[ican]->Clear();
    cans[ican]->Divide(nps,nps);
    char rname[100];
    char sname[100];
    char cname[100];
    for(unsigned iplot=0;iplot<nps*nps;++iplot){
      nplots++;
      unsigned jplot = ican*nps*nps+1+iplot;
      snprintf(rname,100,"rawptspectrum%i",jplot);
      snprintf(sname,100,"looseptspectrum%i",jplot);
      snprintf(cname,100,"convptspectrum%i",jplot);
      TH2F* rh = tdir->Get(rname);
      TH2F* th = tdir->Get(sname);
      TH2F* ch = tdir->Get(cname);
      if(rh != 0 && th != 0 && ch != 0){
        rh->SetStats(0);
        th->SetStats(0);
        ch->SetStats(0);
        cans[ican]->cd(iplot+1);
//        th->SetLineWidth(2);
	th->SetMarkerColor(kGreen);
        th->SetMarkerStyle(2);
// make an absolute threshold
	double maxn = th.GetMaximum();
	double thresh = minn/maxn;
	tp2.Search(th,sigma,"nobackground nomarkov nodraw",thresh);
 	th->SetTitle("Time Spectrum;nsec;#phi");
        th->Draw();
//       rh->GetXaxis()->SetTitle("nsec");
  //      rh->GetYaxis()->SetTitle("#phi");
        rh->Draw("same");
	ch->SetMarkerColor(kRed);
        ch->SetMarkerStyle(5);
        ch->Draw("same");
        if(first){
          first = false;
          TLegend* leg = new TLegend(0.0,0.7,0.3,0.9);
          leg->AddEntry(rh,"All hits","P");
          leg->AddEntry(th,"Selected hits","P");
          leg->AddEntry(ch,"Conversion hits","P");
          leg->AddEntry(dummy,"TSpectrum Peak","P");
          leg->Draw();
        }
      } else {
        moreplots = false;
        break;
      }
    }
    char fname[50];
    snprintf(fname,20,"shtime2d%i.png",ican);
    cans[ican]->Flush();
    cans[ican]->Update();
    cans[ican]->SaveAs(fname);
    ican++;
  }
}
