void StereoTest(TTree* shdiag,const char* page="events") {
  TString spage(page);
  TCut stereoit("dz<100");
  TCut convhit("mcgen==2");
  TCut dhit("mcgen<0&&mcproc==12");
  unsigned ilast(4);
  unsigned ifirst(3);
  if(spage=="events"){
    double rmax(750);
    const size_t nevt = ilast-ifirst+1;
    size_t nbins(200);

    TH2F* shpos[nevt]; 
    TH2F* shpos[nevt]; 
    TH2F* nshpos[nevt];
    TH2F* mcpos[nevt];
    TCanvas* shcans[nevt];
    for(size_t ind =0;ind<nevt;++ind){
      size_t ievt = ifirst+ind;
      char name[80];
      char title[80];

      snprintf(name,80,"shpos_%i",ievt);
      shpos[ind] = new TH2F(name,"Straw Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shpos_%i",ievt);
      shpos[ind] = new TH2F(name,"Stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"nshpos_%i",ievt);
      nshpos[ind] = new TH2F(name,"Non-stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"mcpos_%i",ievt);
      mcpos[ind] = new TH2F(name,"MC Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      shpos[ind]->SetStats(0);
      shpos[ind]->SetStats(0);
      nshpos[ind]->SetStats(0);
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
      snprintf(val,100,"shpos.y:shpos.x>>shpos_%i",ievt);
      shdiag->Draw(val,evtcut+"stereo>0");
      shcans[ind]->cd(3);
      snprintf(val,100,"shpos.y:shpos.x>>nshpos_%i",ievt);
      shdiag->Draw(val,evtcut+"stereo<=0");
      shcans[ind]->cd(4);
      snprintf(val,100,"mcshpos.y:mcshpos.x>>mcpos_%i",ievt);
      shdiag->Draw(val,evtcut);
    }
  } else if(page=="dz") {
    TH1F* zsep = new TH1F("zsep","Stereo hit z separation;#Delta z (mm)",100,0,80);
    TProfile2D* deltaz = new TProfile2D("deltaz","Stereo #Delta z vs position;x (mm);y (mm)",100,-750,750,100,-750,750);
    deltaz->SetStats(0);
    s->Project("zsep","dz",stereoit);
    s->Project("deltaz","dz:shpos.y:shpos.x",stereoit);
    TCanvas* dz = new TCanvas("dz","dz",1200,600);
    dz->Divide(2,1);
    dz->cd(1);
    zsep->Draw();
    dz->cd(2);
    deltaz->Draw("colorz");


  } else if(page=="dzres"){
    TH1F* sdpc1 = new TH1F("sdpc1","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc2 = new TH1F("sdpc2","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc3 = new TH1F("sdpc3","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpd1 = new TH1F("sdpd1","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd2 = new TH1F("sdpd2","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd3 = new TH1F("sdpd3","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH2F* ceres = new TH2F("ceres","#phi resolution vs #Delta z, CE",20,10,60,100,-.3,.3);
    TH2F* deres = new TH2F("deres","#phi resolution vs #Delta z, #delta-ray",20,10,60,100,-.1,.1);
    
    sdpc1->SetLineColor(kRed);
    sdpd1->SetLineColor(kRed);
    sdpc2->SetLineColor(kBlue);
    sdpd2->SetLineColor(kBlue);
    sdpc3->SetLineColor(kGreen);
    sdpd3->SetLineColor(kGreen); 
    sdpc1->SetStats(0);
    sdpc2->SetStats(0);
    sdpc3->SetStats(0);
    sdpd1->SetStats(0);
    sdpd2->SetStats(0);
    sdpd3->SetStats(0); 
    ceres->SetStats(0); 
    deres->SetStats(0); 
    s->Project("sdpc1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+convhit+"dz<30");
    s->Project("sdpc2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+convhit+"dz>30&&dz<50");
    s->Project("sdpc3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+convhit+"dz>50");
    s->Project("sdpd1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+dhit+"dz<30");
    s->Project("sdpd2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+dhit+"dz>30&&dz<50");
    s->Project("sdpd3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereoit+dhit+"dz>50");

    s->Project("ceres","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x):dz",stereoit+convhit);
    s->Project("deres","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x):dz",stereoit+dhit);

    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(sdpc1,"dz < 30","l");
    leg->AddEntry(sdpc2,"30 < dz < 50","l");
    leg->AddEntry(sdpc3," dz > 50","l");
    sdpc1->Scale(1.0/sdpc1->GetEntries());
    sdpc2->Scale(1.0/sdpc2->GetEntries());
    sdpc3->Scale(1.0/sdpc3->GetEntries());
    sdpd1->Scale(1.0/sdpd1->GetEntries());
    sdpd2->Scale(1.0/sdpd2->GetEntries());
    sdpd3->Scale(1.0/sdpd3->GetEntries());

    TCanvas* srcan = new TCanvas("srcan","Stereo resolution",1000,1000);
    srcan->Divide(2,2);
    srcan->cd(1);
    sdpd1->Draw();
    sdpd2->Draw("same");
    sdpd3->Draw("same");
    srcan->cd(2);
    sdpc1->Draw();
    sdpc2->Draw("same");
    sdpc3->Draw("same");
    leg->Draw();
    srcan->cd(3);
    gPad->SetLogz();
    deres->Draw("colorz");
    srcan->cd(4);
    gPad->SetLogz();
    ceres->Draw("colorz");
  }

}
