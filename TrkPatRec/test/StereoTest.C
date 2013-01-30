void StereoTest(TTree* shdiag,size_t ifirst, size_t ilast) { 
  double rmax(750);
  const size_t nevt = ilast-ifirst+1;
  size_t nbins(200);

  TH2F* shpos[nevt]; 
  TH2F* stpos[nevt]; 
  TH2F* nstpos[nevt];
  TH2F* mcpos[nevt];
  TCanvas* shcans[nevt];
  for(size_t ind =0;ind<nevt;++ind){
    size_t ievt = ifirst+ind;
    char name[80];
    char title[80];

    snprintf(name,80,"shpos_%i",ievt);
    shpos[ind] = new TH2F(name,"Straw Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
    snprintf(name,80,"stpos_%i",ievt);
    stpos[ind] = new TH2F(name,"Stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
    snprintf(name,80,"nstpos_%i",ievt);
    nstpos[ind] = new TH2F(name,"Non-stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
    snprintf(name,80,"mcpos_%i",ievt);
    mcpos[ind] = new TH2F(name,"MC Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
    shpos[ind]->SetStats(0);
    stpos[ind]->SetStats(0);
    nstpos[ind]->SetStats(0);
    mcpos[ind]->SetStats(0);
  
    snprintf(name,80,"can_%i",ievt);
    snprintf(title,80,"Tracker Hit Positions for Event %i",ievt);
    shcans[ind] = new TCanvas(name,title,800,800);
    shcans[ind]->Clear();
    shcans[ind]->Divide(2,2);

    char cutn[80];
    snprintf(cutn,80,"eventid==%i",ievt);
    TCut evtcut(cutn);

    char val[100];
    shcans[ind]->cd(1);
    snprintf(val,100,"shpos.y:shpos.x>>shpos_%i",ievt);
    shdiag->Draw(val,evtcut);
    shcans[ind]->cd(2);
    snprintf(val,100,"stpos.y:stpos.x>>stpos_%i",ievt);
    shdiag->Draw(val,evtcut+"stereoh>0");
    shcans[ind]->cd(3);
    snprintf(val,100,"shpos.y:shpos.x>>nstpos_%i",ievt);
    shdiag->Draw(val,evtcut+"stereoh<=0");
    shcans[ind]->cd(4);
    snprintf(val,100,"mcshpos.y:mcshpos.x>>mcpos_%i",ievt);
    shdiag->Draw(val,evtcut);
  }
}


