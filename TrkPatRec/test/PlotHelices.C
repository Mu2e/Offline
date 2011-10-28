void PlotHelices(TDirectory* tdir,unsigned nmax=20, unsigned nps=3){
  gStyle->SetOptStat(0);
  bool moreplots(true);
  unsigned ican(0);
  TCanvas* cans[50];
  unsigned nplots(0);
  while(moreplots && nplots < nmax) {
    bool first(true);
    char cname[20];
    snprintf(cname,20,"timecan%i",ican);
    cans[ican] = new TCanvas(cname,"Straw Hit Positions",800,600);
    cans[ican]->Clear();
    cans[ican]->Divide(nps,2);
    char xyname[100];
    char fzname[100];
    for(unsigned iplot=0;iplot<nps;++iplot){
      nplots++;
      unsigned jplot = ican*nps+1+iplot;
      snprintf(xyname,100,"shxy%i",jplot);
      snprintf(fzname,100,"shphiz%i",jplot);
      TH1F* shxy = tdir->Get(xyname);
      TH1F* shfz = tdir->Get(fzname);
      if(shxy != 0 && shfz != 0){
        shxy->SetStats(0);
        shxy->GetXaxis()->SetTitle("mm");
        shxy->GetYaxis()->SetTitle("mm");        
        shfz->SetStats(0);
        shfz->GetXaxis()->SetTitle("mm");
        shfz->GetYaxis()->SetTitle("radians");
        
        cans[ican]->cd(iplot+1);
        shxy->Draw();
        cans[ican]->cd(iplot+nps+1);
        shfz->Draw();        
      } else {
        moreplots = false;
        break;
      }
    }
    char fname[50];
    snprintf(fname,20,"shpos%i.png",ican);
    cans[ican]->SaveAs(fname);
    ican++;
  }
}
