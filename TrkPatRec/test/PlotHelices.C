void PlotHelices(TDirectory* tdir,unsigned nmax=20, unsigned nps=3){
  gStyle->SetOptStat(0);
  unsigned ican(0);
  unsigned iplot(0);
  TCanvas* cans[50];
  unsigned nplots(0);
  unsigned ntest(0);
  while( nplots < nmax) {
    ++ntest;
    bool first(true);
    char cname[20];
    snprintf(cname,20,"helixcan%i",ican);
    cans[ican] = new TCanvas(cname,"Straw Hit Positions",800,600);
    cans[ican]->Clear();
    cans[ican]->Divide(nps,2);
    char gxyname[100];
    char bxyname[100];
    char gfzname[100];
    char bfzname[100];
    for(unsigned ievt=0;ievt<1000;++ievt){
      unsigned jplot = 10*ievt;
      snprintf(gxyname,100,"gshxy%i",jplot);
      snprintf(bxyname,100,"bshxy%i",jplot);
      snprintf(gfzname,100,"gshphiz%i",jplot);
      snprintf(bfzname,100,"bshphiz%i",jplot);
      TH1F* gshxy = tdir->Get(gxyname);
      TH1F* bshxy = tdir->Get(bxyname);
      TH1F* gshfz = tdir->Get(gfzname);
      TH1F* bshfz= tdir->Get(bfzname);
      if(gshxy != 0 && bshxy != 0) {
	++nplots;
        gshxy->SetStats(0);
        gshxy->GetXaxis()->SetTitle("mm");
        gshxy->GetYaxis()->SetTitle("mm");        
        gshfz->SetStats(0);
        gshfz->GetXaxis()->SetTitle("mm");
        gshfz->GetYaxis()->SetTitle("radians");
        if(iplot >= nps){
	  iplot = 0;
	  ++ican;
	  snprintf(cname,20,"helixcan%i",ican);
	  cans[ican] = new TCanvas(cname,"Straw Hit Positions",800,600);
	  cans[ican]->Clear();
	  cans[ican]->Divide(nps,2);
	} 
        cans[ican]->cd(iplot+1);
        gshxy->Draw();
        bshxy->Draw("same");
	if(gshfz != 0 && bshfz != 0){
	  cans[ican]->cd(iplot+nps+1);
	  gshfz->Draw();        
	  bshfz->Draw("same");
	}
	++iplot;
      }
    }
    char fname[50];
    snprintf(fname,20,"shpos%i.png",ican);
    cans[ican]->SaveAs(fname);
    ican++;
  }
}
