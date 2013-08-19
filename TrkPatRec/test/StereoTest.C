#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TRandom3.h"
#include <iostream>
#include <math.h>
#include <vector>

void StereoTest(TTree* shdiag,const char* page="events") {
  TString spage(page);
  TCut stereohit("stereo!=0");
  TCut convhit("mcgen==2");
  TCut dhit("mcgen<0&&mcproc==12");
  unsigned ilast(4);
  unsigned ifirst(3);
  if(spage=="events"){
    double rmax(750);
    const size_t nevt = ilast-ifirst+1;
    size_t nbins(200);

    TH2F* shpos[nevt]; 
    TH2F* nshpos[nevt];
    TH2F* mcpos[nevt];
    TCanvas* shcans[nevt];
    for(unsigned ind =0;ind<nevt;++ind){
      unsigned ievt = ifirst+ind;
      char name[80];
      char title[80];

      snprintf(name,80,"shpos_%u",ievt);
      shpos[ind] = new TH2F(name,"Straw Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shpos_%u",ievt);
      shpos[ind] = new TH2F(name,"Stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"nshpos_%u",ievt);
      nshpos[ind] = new TH2F(name,"Non-stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"mcpos_%u",ievt);
      mcpos[ind] = new TH2F(name,"MC Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      shpos[ind]->SetStats(0);
      shpos[ind]->SetStats(0);
      nshpos[ind]->SetStats(0);
      mcpos[ind]->SetStats(0);

      snprintf(name,80,"can_%u",ievt);
      snprintf(title,80,"Tracker Hit Positions for Event %u",ievt);
      shcans[ind] = new TCanvas(name,title,800,800);
      shcans[ind]->Clear();
      shcans[ind]->Divide(2,2);

      char cutn[80];
      snprintf(cutn,80,"eventid==%u",ievt);
      TCut evtcut(cutn);

      char val[100];
      shcans[ind]->cd(1);
      snprintf(val,100,"shpos.y:shpos.x>>shpos_%u",ievt);
      shdiag->Draw(val,evtcut);
      shcans[ind]->cd(2);
      snprintf(val,100,"shpos.y:shpos.x>>shpos_%u",ievt);
      shdiag->Draw(val,evtcut+"stereo>0");
      shcans[ind]->cd(3);
      snprintf(val,100,"shpos.y:shpos.x>>nshpos_%u",ievt);
      shdiag->Draw(val,evtcut+"stereo<=0");
      shcans[ind]->cd(4);
      snprintf(val,100,"mcshpos.y:mcshpos.x>>mcpos_%u",ievt);
      shdiag->Draw(val,evtcut);
    }
  } else if(spage=="perr") {
    TH1F* perr = new TH1F("perr","Hit #phi error;#phi error (mm)",100,0,45);
    perr->SetStats(0);
    TProfile2D* pperr = new TProfile2D("pperr","Hit #phi error vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    pperr->SetStats(0);
    shdiag->Project("perr","pres");
    shdiag->Project("pperr","pres:shpos.y:shpos.x");
    pperr->SetMaximum(45);
    TCanvas* cperr = new TCanvas("cperr","cperr",1200,600);
    cperr->Divide(2,1);
    cperr->cd(1);
    perr->Draw();
    cperr->cd(2);
    pperr->Draw("colorz");

  } else if(spage=="rerr"){
    TH1F* rerr = new TH1F("rerr","Hit #rho error ;#rho error (mm)",100,0,45);
    rerr->SetStats(0);
    TProfile2D* prerr = new TProfile2D("prerr","Hit #rho error vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    prerr->SetStats(0);
    shdiag->Project("rerr","rres");
    shdiag->Project("prerr","rres:shpos.y:shpos.x");
    prerr->SetMaximum(45);
    TCanvas* crerr = new TCanvas("crerr","crerr",1200,600);
    crerr->Divide(2,1);
    crerr->cd(1);
    rerr->Draw();
    crerr->cd(2);
    prerr->Draw("colorz");

  } else if(spage=="dz") {
    TH1F* dist = new TH1F("dist","Hit #Delta z;#Delta z (mm)",100,-10,45);
    dist->SetStats(0);
    TProfile2D* pdist = new TProfile2D("pdist","Hit #Delta z vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    pdist->SetStats(0);
    shdiag->Project("dist","dist");
    shdiag->Project("pdist","dist:shpos.y:shpos.x");
    pdist->SetMaximum(45);
    pdist->SetMinimum(-5);
    TCanvas* cdist = new TCanvas("cdist","cdist",1200,600);
    cdist->Divide(2,1);
    cdist->cd(1);
    dist->Draw();
    cdist->cd(2);
    pdist->Draw("colorz");

  } else if(spage=="sfrac"){

    TH2F* apos = new TH2F("apos","True hit position;x (mm);y (mm)",100,-750,750,100,-750,750);
    TH2F* spos = new TH2F("spos","Stereo fraction vs true hit position;x (mm);y (mm)",100,-750,750,100,-750,750);
    TH2F* nspos = new TH2F("nspos","Non-stereo fraction vs true hit position;x (mm);y (mm)",100,-750,750,100,-750,750);
  
    shdiag->Project("apos","mcshpos.y:mcshpos.x");
    shdiag->Project("spos","mcshpos.y:mcshpos.x",stereohit);
    shdiag->Project("nspos","mcshpos.y:mcshpos.x",!stereohit);
    spos->SetStats(0);
    nspos->SetStats(0);
    spos->Divide(apos);
    nspos->Divide(apos);
    spos->SetMaximum(1.0);
    spos->SetMinimum(0.0);
    nspos->SetMaximum(1.0);
    nspos->SetMinimum(0.0);
    TCanvas* scan = new TCanvas("scan","scan",1200,600);
    scan->Divide(2,1);
    scan->cd(1);
    spos->Draw("colorz");
    scan->cd(2);
    nspos->Draw("colorz");

  } else if(spage=="res"){
    TH1F* sdpc1 = new TH1F("sdpc1","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc2 = new TH1F("sdpc2","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc3 = new TH1F("sdpc3","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc4 = new TH1F("sdpc4","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpd1 = new TH1F("sdpd1","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd2 = new TH1F("sdpd2","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd3 = new TH1F("sdpd3","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd4 = new TH1F("sdpd4","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH2F* ceres = new TH2F("ceres","#phi resolution vs #Delta z, CE",20,-5,45,100,-.5,.5);
    TH2F* deres = new TH2F("deres","#phi resolution vs #Delta z, #delta-ray",20,-5,45,100,-.5,.5);
    
    sdpc1->SetLineColor(kRed);
    sdpd1->SetLineColor(kRed);
    sdpc2->SetLineColor(kBlue);
    sdpd2->SetLineColor(kBlue);
    sdpc3->SetLineColor(kGreen);
    sdpd3->SetLineColor(kGreen); 
    sdpc3->SetLineColor(kCyan);
    sdpd3->SetLineColor(kCyan); 
    sdpc1->SetStats(0);
    sdpc2->SetStats(0);
    sdpc3->SetStats(0);
    sdpd1->SetStats(0);
    sdpd2->SetStats(0);
    sdpd3->SetStats(0); 
    ceres->SetStats(0); 
    deres->SetStats(0); 

    TCut nst("dist<0");
    TCut s1("dist>0&&dist<18");
    TCut s2("dist>18&&dist<30");
    TCut s3("dist>30");

    shdiag->Project("sdpc1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+nst);
    shdiag->Project("sdpc2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+s1);
    shdiag->Project("sdpc3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+s2);
    shdiag->Project("sdpc4","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+s3);
    shdiag->Project("sdpd1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+nst);
    shdiag->Project("sdpd2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+s1);
    shdiag->Project("sdpd3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+s2);
    shdiag->Project("sdpd4","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+s3);

    shdiag->Project("ceres","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x):dist",convhit);
    shdiag->Project("deres","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x):dist",dhit);

    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(sdpc1,"non-stereo","l");
    leg->AddEntry(sdpc2,"#Delta z<18","l");
    leg->AddEntry(sdpc3,"18 < #Delta z < 30","l");
    leg->AddEntry(sdpc4,"#Delta z > 30","l");
    sdpc1->Scale(1.0/sdpc1->GetEntries());
    sdpc2->Scale(1.0/sdpc2->GetEntries());
    sdpc3->Scale(1.0/sdpc3->GetEntries());
    sdpc4->Scale(1.0/sdpc4->GetEntries());
    sdpd1->Scale(1.0/sdpd1->GetEntries());
    sdpd2->Scale(1.0/sdpd2->GetEntries());
    sdpd3->Scale(1.0/sdpd3->GetEntries());
    sdpd4->Scale(1.0/sdpd4->GetEntries());

    TCanvas* srcan = new TCanvas("srcan","Stereo resolution",1000,1000);
    srcan->Divide(2,2);
    srcan->cd(1);
    sdpd2->Draw();
    sdpd1->Draw("same");
    sdpd3->Draw("same");
    sdpd4->Draw("same");
    srcan->cd(2);
    sdpc2->Draw();
    sdpc1->Draw("same");
    sdpc3->Draw("same");
    sdpc4->Draw("same");
    leg->Draw();
    srcan->cd(3);
    gPad->SetLogz();
    deres->Draw("colorz");
    srcan->cd(4);
    gPad->SetLogz();
    ceres->Draw("colorz");
  }
}

