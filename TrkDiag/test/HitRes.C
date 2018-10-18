#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TCut.h"


void HitRes(TTree* shdiag) {
  TCut stereo("stereo&&mcproc<20");
  TCut nstereo("(!stereo)&&mcproc<20");
  TH1F* sdp = new TH1F("sdp","Stereo #Delta#phi;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
  TH1F* ndp = new TH1F("ndp","Non-stereo #Delta#phi;#phi_{sh}-#phi_{MC}",100,-0.4,0.4);
  TH1F* sdrp = new TH1F("sdrp","Stereo #Delta#rho#phi;#rho#phi_{sh}-#rho#phi_{MC}",100,-80.0,80.0);
  TH1F* ndrp = new TH1F("ndrp","Non-stereo #Delta#rho#phi;#rho#phi_{sh}-#rho#phi_{MC}",100,-150.0,150.0);
  TH1F* sdr = new TH1F("sdr","Stereo #Delta#rho;#rho_{sh}-#rho_{MC}",100,-40,40);
  TH1F* ndr = new TH1F("ndr","Non-stereo #Delta#rho;#rho_{sh}-#rho_{MC}",100,-100,100);
//  TH1F* sdt = new TH1F("sdt","Stereo #Delta t;t_{sh}-t_{MC}",100,-100,100);
//  TH1F* ndt = new TH1F("ndt","Non-stereo #Delta t;t_{sh}-t_{MC}",100,-100,100);
  TH1F* sdx = new TH1F("sdx","Stereo #Delta x;x_{sh}-x_{MC}",100,-50,50);
  TH1F* ndx = new TH1F("ndx","Non-stereo #Delta x;x_{sh}-x_{MC}",100,-150,150);
  shdiag->Project("sdp","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",stereo);
  shdiag->Project("ndp","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",nstereo);
  shdiag->Project("sdrp","atan2(shpos.y,shpos.x)*sqrt(shpos.y^2+shpos.x^2)-atan2(mcshpos.y,mcshpos.x)*sqrt(mcshpos.y^2+mcshpos.x^2)",stereo);
  shdiag->Project("ndrp","atan2(shpos.y,shpos.x)*sqrt(shpos.y^2+shpos.x^2)-atan2(mcshpos.y,mcshpos.x)*sqrt(mcshpos.y^2+mcshpos.x^2)",nstereo);
  shdiag->Project("sdr","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",stereo);
  shdiag->Project("ndr","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",nstereo);
//  shdiag->Project("sdt","time-mctime",stereo);
//  shdiag->Project("ndt","time-mctime",nstereo);
  shdiag->Project("sdx","shpos.x-mcshpos.x",stereo);
  shdiag->Project("ndx","shpos.x-mcshpos.x",nstereo);

  TCanvas* dhcan = new TCanvas("hcan","Hit Resolution",1200,800);
  dhcan->Divide(4,2);
  dhcan->cd(1);
  sdp->Fit("gaus");
  dhcan->cd(5);
  ndp->Fit("gaus");
  dhcan->cd(2);
  sdrp->Fit("gaus");
  dhcan->cd(6);
  ndrp->Fit("gaus");
  dhcan->cd(3);
  sdr->Fit("gaus");
  dhcan->cd(7);
  ndr->Fit("gaus");
  dhcan->cd(4);
  sdx->Fit("gaus");
  dhcan->cd(8);
  ndx->Fit("gaus");
}

