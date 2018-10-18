void DeltaHitRes(TTree* ddiag) {
  TCut deltah("_mcgen<0&&abs(_mcpdg)==11");
  TCut deltap("nwide>50&&nconv==0&&nprimary==nwide");
  TCut stereo("_stereo");
  TCut nstereo("!_stereo");
  TH1F* sdp = new TH1F("sdp","Stereo #Delta#phi;#Delta#phi",100,-0.2,0.2);
  TH1F* ndp = new TH1F("ndp","Non-stereo #Delta#phi;#Delta#phi",100,-0.4,0.4);
  TH1F* sdr = new TH1F("sdr","Stereo #Delta#rho;#Delta#rho",100,-40,40);
  TH1F* ndr = new TH1F("ndr","Non-stereo #Delta#rho;#Delta#rho",100,-40,40);
  TH1F* sdt = new TH1F("sdt","Stereo #Delta t;#Delta t",100,-100,100);
  TH1F* ndt = new TH1F("ndt","Non-stereo #Delta t;#Delta t",100,-100,100);
  ddiag->Project("sdp","atan2(_pos.dy,_pos.dx)-pmed",deltap+deltah+stereo);
  ddiag->Project("ndp","atan2(_pos.dy,_pos.dx)-pmed",deltap+deltah+nstereo);
  ddiag->Project("sdr","sqrt(_pos.dy^2+_pos.dx^2)-rmed",deltap+deltah+stereo);
  ddiag->Project("ndr","sqrt(_pos.dy^2+_pos.dx^2)-rmed",deltap+deltah+nstereo);
  ddiag->Project("sdt","_time-tmed",deltap+deltah+stereo);
  ddiag->Project("ndt","_time-tmed",deltap+deltah+nstereo);

  TCanvas* dhcan = new TCanvas("dhcan","Delta hits",1200,800);
  dhcan->Divide(3,2);
  dhcan->cd(1);
  sdp->Fit("gaus");
  dhcan->cd(2);
  ndp->Fit("gaus");
  dhcan->cd(3);
  sdr->Fit("gaus");
  dhcan->cd(4);
  ndr->Fit("gaus");
  dhcan->cd(5);
  sdt->Fit("gaus");
  dhcan->cd(6);
  ndt->Fit("gaus");
}

