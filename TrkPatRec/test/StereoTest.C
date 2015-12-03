#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TTree.h"
#include "THStack.h"
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

void StereoTest(TTree* shdiag,const char* page="events",const char* cutstring="") {
  TString spage(page);
  TCut stereohit("stereo!=0");
  TCut convhit("mcgen==2&&pmom>100.0&&tsel");
  TCut dhit("mcproc<20");
  TCut addcut(cutstring);
  TCut evenstation("(plane/2)%2==0");
  TCut oddstation("(plane/2)%2==1");
  if(strcmp(cutstring,"")!= 0){
    convhit += addcut;
    dhit += addcut;
  }

  unsigned ilast(10);
  unsigned ifirst(1);
  if(spage=="bkgevents"){
    double rmax(750);
    const size_t nevt = ilast-ifirst+1;
    size_t nbins(200);

    TH2F* shpos[nevt]; 
    TH2F* shpose[nevt]; 
    TH2F* shposr[nevt];
    TH2F* shposd[nevt];
    TH2F* shposi[nevt];
    TH2F* shposc[nevt];
    TCanvas* shcans[nevt];
    for(unsigned ind =0;ind<nevt;++ind){
      unsigned ievt = ifirst+ind;
      char name[80];
      char title[80];

      snprintf(name,80,"shpos_%u",ievt);
      shpos[ind] = new TH2F(name,"Straw Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shpose_%u",ievt);
      shpose[ind] = new TH2F(name,"Straw Hit Transverse Position, esel;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shposr_%u",ievt);
      shposr[ind] = new TH2F(name,"Straw Hit Transverse Position, rsel;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shposd_%u",ievt);
      shposd[ind] = new TH2F(name,"Straw Hit Transverse Position, !delta;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shposi_%u",ievt);
      shposi[ind] = new TH2F(name,"Straw Hit Transverse Position, isolated;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"shposc_%u",ievt);
      shposc[ind] = new TH2F(name,"Straw Hit Transverse Position, clean;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      shpos[ind]->SetStats(0);
      shpose[ind]->SetStats(0);
      shposr[ind]->SetStats(0);
      shposd[ind]->SetStats(0);
      shposi[ind]->SetStats(0);
      shposc[ind]->SetStats(0);

      snprintf(name,80,"bkgevtcan_%u",ievt);
      snprintf(title,80,"Tracker Hit Positions for Event %u",ievt);
      shcans[ind] = new TCanvas(name,title,1000,660);
      shcans[ind]->Clear();
      shcans[ind]->Divide(3,2);

      char cutn[80];
      snprintf(cutn,80,"eventid==%u",ievt);
      TCut evtcut(cutn);

      char val[100];
      shcans[ind]->cd(1);
      snprintf(val,100,"shpos.y:shpos.x>>shpos_%u",ievt);
      shdiag->Draw(val,evtcut);
      shcans[ind]->cd(2);
      snprintf(val,100,"shpos.y:shpos.x>>shpose_%u",ievt);
      shdiag->Draw(val,evtcut+"esel");
      shcans[ind]->cd(3);
      snprintf(val,100,"shpos.y:shpos.x>>shposr_%u",ievt);
      shdiag->Draw(val,evtcut+"rsel");
      shcans[ind]->cd(4);
      snprintf(val,100,"shpos.y:shpos.x>>shposd_%u",ievt);
      shdiag->Draw(val,evtcut+"(!delta)");
      shcans[ind]->cd(5);
      snprintf(val,100,"shpos.y:shpos.x>>shposi_%u",ievt);
      shdiag->Draw(val,evtcut+"isolated");
      shcans[ind]->cd(6);
      snprintf(val,100,"shpos.y:shpos.x>>shposc_%u",ievt);
      shdiag->Draw(val,evtcut+"(!delta)&&(!isolated)&&esel&&rsel");
    }
  } else if(spage=="events"){
    double rmax(750);
    const size_t nevt = ilast-ifirst+1;
    size_t nbins(200);

    TH2F* shpos[nevt]; 
    TH2F* stpos[nevt]; 
    TH2F* nshpos[nevt];
    TH2F* mcpos[nevt];
    TCanvas* shcans[nevt];
    for(unsigned ind =0;ind<nevt;++ind){
      unsigned ievt = ifirst+ind;
      char name[80];
      char title[80];

      snprintf(name,80,"shpos_%u",ievt);
      shpos[ind] = new TH2F(name,"Straw Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"stpos_%u",ievt);
      stpos[ind] = new TH2F(name,"Stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"nshpos_%u",ievt);
      nshpos[ind] = new TH2F(name,"Non-stereo Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      snprintf(name,80,"mcpos_%u",ievt);
      mcpos[ind] = new TH2F(name,"MC Hit Transverse Position;x (mm);y (mm)",nbins,-rmax,rmax,nbins,-rmax,rmax);
      shpos[ind]->SetStats(0);
      stpos[ind]->SetStats(0);
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
      snprintf(val,100,"shpos.y:shpos.x>>stpos_%u",ievt);
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
    TH1F* hdz = new TH1F("hdz","Hit #Delta z;#Delta z (mm)",100,0,70);
//    dist->SetStats(0);
    TProfile2D* evposdz = new TProfile2D("evposdz","Even Station #Delta z vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    TProfile2D* odposdz = new TProfile2D("odposdz","Odd Station #Delta z vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    TProfile2D* alposdz = new TProfile2D("alposdz","All Stations #Delta z vs reco position;x (mm);y (mm)",100,-800,800,100,-800,800);
    evposdz->SetStats(0);
    odposdz->SetStats(0);
    alposdz->SetStats(0);
    shdiag->Project("hdz","dist",stereohit+addcut);
    shdiag->Project("evposdz","dist:shpos.y:shpos.x",stereohit+evenstation+addcut);
    shdiag->Project("odposdz","dist:shpos.y:shpos.x",stereohit+oddstation+addcut);
    shdiag->Project("alposdz","dist:shpos.y:shpos.x",stereohit+addcut);
    evposdz->SetMaximum(65);
    evposdz->SetMinimum(0);
    odposdz->SetMaximum(65);
    odposdz->SetMinimum(0);
    alposdz->SetMaximum(65);
    alposdz->SetMinimum(0);
    TCanvas* cdz = new TCanvas("cdz","cdz",1000,1000);
    cdz->Divide(2,2);
    cdz->cd(1);
    alposdz->Draw("colorz");
    cdz->cd(2);
    evposdz->Draw("colorz");
    cdz->cd(3);
    odposdz->Draw("colorz");
    cdz->cd(4);
    hdz->Draw();

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
    TH1F* sdpc0 = new TH1F("sdpc0","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc1 = new TH1F("sdpc1","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc2 = new TH1F("sdpc2","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc3 = new TH1F("sdpc3","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    TH1F* sdpc4 = new TH1F("sdpc4","#phi resolution, CE;#phi_{sh}-#phi_{MC}",100,-0.3,0.3);
    THStack* sdpc = new THStack("sdpc","#phi resolution, CE;#phi_{sh}-#phi_{MC}");
    sdpc->Add(sdpc0);
    sdpc->Add(sdpc4);
    sdpc->Add(sdpc3);
    sdpc->Add(sdpc2);
//    sdpc->Add(sdpc1);

    TH1F* sdpd0 = new TH1F("sdpd0","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd1 = new TH1F("sdpd1","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd2 = new TH1F("sdpd2","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd3 = new TH1F("sdpd3","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    TH1F* sdpd4 = new TH1F("sdpd4","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}",100,-0.1,0.1);
    THStack* sdpd = new THStack("sdpd","#phi resolution, #delta-ray;#phi_{sh}-#phi_{MC}");
    sdpd->Add(sdpd0);
    sdpd->Add(sdpd4);
    sdpd->Add(sdpd3);
    sdpd->Add(sdpd2);
//    sdpd->Add(sdpd1);

    TH1F* sdrc0 = new TH1F("sdrc0","#rho resolution, CE;#rho_{sh}-#rho_{MC}",100,-50,50);
    TH1F* sdrc1 = new TH1F("sdrc1","#rho resolution, CE;#rho_{sh}-#rho_{MC}",100,-50,50);
    TH1F* sdrc2 = new TH1F("sdrc2","#rho resolution, CE;#rho_{sh}-#rho_{MC}",100,-50,50);
    TH1F* sdrc3 = new TH1F("sdrc3","#rho resolution, CE;#rho_{sh}-#rho_{MC}",100,-50,50);
    TH1F* sdrc4 = new TH1F("sdrc4","#rho resolution, CE;#rho_{sh}-#rho_{MC}",100,-50,50);
    THStack* sdrc = new THStack("sdrc","#rho resolution, CE;#rho_{sh}-#rho_{MC}");
    sdrc->Add(sdrc0);
    sdrc->Add(sdrc4);
    sdrc->Add(sdrc3);
    sdrc->Add(sdrc2);
//    sdrc->Add(sdrc1);

    TH1F* sdrd0 = new TH1F("sdrd0","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}",100,-20,20);
    TH1F* sdrd1 = new TH1F("sdrd1","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}",100,-20,20);
    TH1F* sdrd2 = new TH1F("sdrd2","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}",100,-20,20);
    TH1F* sdrd3 = new TH1F("sdrd3","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}",100,-20,20);
    TH1F* sdrd4 = new TH1F("sdrd4","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}",100,-20,20);
    THStack* sdrd = new THStack("sdrd","#rho resolution, #delta-ray;#rho_{sh}-#rho_{MC}");
    sdrd->Add(sdrd0);
    sdrd->Add(sdrd4);
    sdrd->Add(sdrd3);
    sdrd->Add(sdrd2);
//    sdrd->Add(sdrd1);

    sdpc0->SetFillColor(kYellow);
    sdpd0->SetFillColor(kYellow);
    sdpc1->SetLineColor(kBlack);
    sdpd1->SetLineColor(kBlack);
    sdpc2->SetFillColor(kCyan);
    sdpd2->SetFillColor(kCyan);
    sdpc3->SetFillColor(kGreen);
    sdpd3->SetFillColor(kGreen); 
    sdpc4->SetFillColor(kRed);
    sdpd4->SetFillColor(kRed); 
    sdpc0->SetStats(0);

    sdrc0->SetFillColor(kYellow);
    sdrd0->SetFillColor(kYellow);
    sdrc1->SetLineColor(kBlack);
    sdrd1->SetLineColor(kBlack);
    sdrc2->SetFillColor(kCyan);
    sdrd2->SetFillColor(kCyan);
    sdrc3->SetFillColor(kGreen);
    sdrd3->SetFillColor(kGreen); 
    sdrc4->SetFillColor(kRed);
    sdrd4->SetFillColor(kRed); 
//    sdrc0->SetStats(0);
//    sdpc1->SetStats(0);
//    sdpc2->SetStats(0);
//    sdpc3->SetStats(0);
//    sdpc4->SetStats(0);
//    sdpd0->SetStats(0);
//    sdpd1->SetStats(0);
//    sdpd2->SetStats(0);
//    sdpd3->SetStats(0); 
//    sdpd4->SetStats(0); 

//    sdrc1->SetStats(0);
//    sdrc2->SetStats(0);
//    sdrc3->SetStats(0);
//    sdrc4->SetStats(0);
//    sdrd0->SetStats(0);
//    sdrd1->SetStats(0);
//    sdrd2->SetStats(0);
//    sdrd3->SetStats(0); 
//    sdrd4->SetStats(0); 

    TCut st0("dist<0");
    TCut st1("");
    TCut st2("dist>0&&dist<30");
    TCut st3("dist>30&&dist<50");
    TCut st4("dist>50");
    shdiag->Project("sdpc0","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+st0);
    shdiag->Project("sdpc1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+st1);
    shdiag->Project("sdpc2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+st2);
    shdiag->Project("sdpc3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+st3);
    shdiag->Project("sdpc4","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",convhit+st4);
    shdiag->Project("sdpd0","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+st0);
    shdiag->Project("sdpd1","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+st1);
    shdiag->Project("sdpd2","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+st2);
    shdiag->Project("sdpd3","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+st3);
    shdiag->Project("sdpd4","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)",dhit+st4);

    shdiag->Project("sdrc0","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",convhit+st0);
    shdiag->Project("sdrc1","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",convhit+st1);
    shdiag->Project("sdrc2","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",convhit+st2);
    shdiag->Project("sdrc3","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",convhit+st3);
    shdiag->Project("sdrc4","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",convhit+st4);
    shdiag->Project("sdrd0","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",dhit+st0);
    shdiag->Project("sdrd1","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",dhit+st1);
    shdiag->Project("sdrd2","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",dhit+st2);
    shdiag->Project("sdrd3","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",dhit+st3);
    shdiag->Project("sdrd4","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)",dhit+st4);

    TLegend* leg = new TLegend(0.1,0.6,0.4,0.9);
    leg->AddEntry(sdpc2,"#Delta z<30","f");
    leg->AddEntry(sdpc3,"30 < #Delta z < 50","f");
    leg->AddEntry(sdpc4,"#Delta z > 50","f");
    leg->AddEntry(sdpc0,"non-stereo","f");
    leg->AddEntry(sdpc1,"all","l");
//    sdpc0->Scale(1.0/sdpc0->GetEntries());
//    sdpc1->Scale(1.0/sdpc1->GetEntries());
//    sdpc2->Scale(1.0/sdpc2->GetEntries());
//    sdpc3->Scale(1.0/sdpc3->GetEntries());
//    sdpc4->Scale(1.0/sdpc4->GetEntries());
//    sdpd0->Scale(1.0/sdpd0->GetEntries());
//    sdpd1->Scale(1.0/sdpd1->GetEntries());
//    sdpd2->Scale(1.0/sdpd2->GetEntries());
//    sdpd3->Scale(1.0/sdpd3->GetEntries());
//    sdpd4->Scale(1.0/sdpd4->GetEntries());
//
//    sdrc0->Scale(1.0/sdrc0->GetEntries());
//    sdrc1->Scale(1.0/sdrc1->GetEntries());
//    sdrc2->Scale(1.0/sdrc2->GetEntries());
//    sdrc3->Scale(1.0/sdrc3->GetEntries());
//    sdrc4->Scale(1.0/sdrc4->GetEntries());
//    sdrd0->Scale(1.0/sdrd0->GetEntries());
//    sdrd1->Scale(1.0/sdrd1->GetEntries());
//    sdrd2->Scale(1.0/sdrd2->GetEntries());
//    sdrd3->Scale(1.0/sdrd3->GetEntries());
//    sdrd4->Scale(1.0/sdrd4->GetEntries());
//    double pmax(0.125);
//    sdpc0->SetMaximum(pmax);
//    sdpc1->SetMaximum(pmax);
//    sdpc2->SetMaximum(pmax);
//    sdpc3->SetMaximum(pmax);
//    sdpc4->SetMaximum(pmax);
//    sdpd0->SetMaximum(pmax);
//    sdpd1->SetMaximum(pmax);
//    sdpd2->SetMaximum(pmax);
//    sdpd3->SetMaximum(pmax);
//    sdpd4->SetMaximum(pmax);
//
//    double rmax(0.1);
//    sdrc0->SetMaximum(rmax);
//    sdrc1->SetMaximum(rmax);
//    sdrc2->SetMaximum(rmax);
//    sdrc3->SetMaximum(rmax);
//    sdrc4->SetMaximum(rmax);
//    sdrd0->SetMaximum(rmax);
//    sdrd1->SetMaximum(rmax);
//    sdrd2->SetMaximum(rmax);
//    sdrd3->SetMaximum(rmax);
//    sdrd4->SetMaximum(rmax);
 
    TCanvas* srcan = new TCanvas("srcan","Stereo resolution",1000,1000);
    srcan->Divide(2,2);
    srcan->cd(1);
    sdpd->Draw();
    sdpd1->Draw("sames");
//    sdpd0->Draw("same");
//    sdpd2->Draw("same");
//    sdpd3->Draw("same");
//    sdpd4->Draw("same");
    leg->Draw();
    srcan->cd(2);
    sdpc->Draw();
    sdpc1->Draw("sames");
//    sdpc0->Draw("same");
//    sdpc2->Draw("same");
//    sdpc3->Draw("same");
//    sdpc4->Draw("same");
    leg->Draw();
    srcan->cd(3);
    sdrd->Draw();
    sdrd1->Draw("sames");
//    sdrd0->Draw("same");
//    sdrd2->Draw("same");
//    sdrd3->Draw("same");
//    sdrd4->Draw("same");
    leg->Draw();
    srcan->cd(4);
    sdrc->Draw();
    sdrc1->Draw("sames");
//    sdrc0->Draw("same");
//    sdrc2->Draw("same");
//    sdrc3->Draw("same");
//    sdrc4->Draw("same");
    leg->Draw();
  } else if (spage == "delta") {
    TH2F* sdce = new TH2F("sdce","#delta-ray hit flaging, true CE;Reco hit #delta flag;Reco hit stereo flag",2,-0.5,1.5,2,-0.5,1.5);
    TH2F* sdde = new TH2F("sdde","#delta-ray hit flaging, true #delta ray;Reco hit #delta flag;Reco hit stereo flag",2,-0.5,1.5,2,-0.5,1.5);
    TH2F* sce = new TH2F("sce","hit flaging, true CE;Reco hit selection flag;Reco hit stereo flag",2,-0.5,1.5,2,-0.5,1.5);
    TH2F* sbk = new TH2F("sbk","hit flaging, non-CE;Reco hit selection flag;Reco hit stereo flag",2,-0.5,1.5,2,-0.5,1.5);
    sdce->GetXaxis()->SetBinLabel(1,"Not Delta");
    sdce->GetXaxis()->SetBinLabel(2,"Delta");
    sdce->GetYaxis()->SetBinLabel(1,"Not Stereo");
    sdce->GetYaxis()->SetBinLabel(2,"Stereo");
    sdde->GetXaxis()->SetBinLabel(1,"Not Delta");
    sdde->GetXaxis()->SetBinLabel(2,"Delta");
    sdde->GetYaxis()->SetBinLabel(1,"Not Stereo");
    sdde->GetYaxis()->SetBinLabel(2,"Stereo");
    sdce->GetXaxis()->SetLabelSize(0.06);
    sdce->GetYaxis()->SetLabelSize(0.06);
    sdde->GetXaxis()->SetLabelSize(0.06);
    sdde->GetYaxis()->SetLabelSize(0.06);
    sdce->GetXaxis()->SetTitleSize(0.05);
    sdce->GetYaxis()->SetTitleSize(0.05);
    sdde->GetXaxis()->SetTitleSize(0.05);
    sdde->GetYaxis()->SetTitleSize(0.05);
    sdce->GetYaxis()->SetTitleOffset(1.5);
    sdde->GetYaxis()->SetTitleOffset(1.5);
    sdce->SetFillColor(kRed);
    sdde->SetFillColor(kCyan);
    sdce->SetStats(0);
    sdde->SetStats(0);
    sdce->SetMarkerSize(3.0);
    sdde->SetMarkerSize(3.0);
    sce->GetXaxis()->SetBinLabel(1,"Selected");
    sce->GetXaxis()->SetBinLabel(2,"Not Selected");
    sce->GetYaxis()->SetBinLabel(1,"Not Stereo");
    sce->GetYaxis()->SetBinLabel(2,"Stereo");
    sbk->GetXaxis()->SetBinLabel(1,"Selected");
    sbk->GetXaxis()->SetBinLabel(2,"Not Selected");
    sbk->GetYaxis()->SetBinLabel(1,"Not Stereo");
    sbk->GetYaxis()->SetBinLabel(2,"Stereo");
    sce->GetXaxis()->SetLabelSize(0.06);
    sce->GetYaxis()->SetLabelSize(0.06);
    sbk->GetXaxis()->SetLabelSize(0.06);
    sbk->GetYaxis()->SetLabelSize(0.06);
    sce->GetXaxis()->SetTitleSize(0.05);
    sce->GetYaxis()->SetTitleSize(0.05);
    sbk->GetXaxis()->SetTitleSize(0.05);
    sbk->GetYaxis()->SetTitleSize(0.05);
    sce->GetYaxis()->SetTitleOffset(1.5);
    sbk->GetYaxis()->SetTitleOffset(1.5);
    sce->SetFillColor(kRed);
    sbk->SetFillColor(kCyan);
    sce->SetStats(0);
    sbk->SetStats(0);
    sce->SetMarkerSize(3.0);
    sbk->SetMarkerSize(3.0);
    TCut esel("esel>0");
    TCut tsel("tsel>0");
    TCut rsel("rsel>0");
    TCut isolated("isolated>0");
    TCut delta("delta>0");

    shdiag->Project("sdce","stereo:delta",convhit+esel+rsel+tsel+(!isolated));
    shdiag->Project("sdde","stereo:delta",dhit+esel+rsel+tsel+(!isolated));

    shdiag->Project("sce","stereo:!(esel&&rsel&&tsel&&(!isolated)&&(!delta))",convhit);
    shdiag->Project("sbk","stereo:!(esel&&rsel&&tsel&&(!isolated)&&(!delta))",!convhit);

    TCanvas* dcan = new TCanvas("dcan","dcan",800,600);
    dcan->Divide(2,2);
    dcan->cd(1);
    gPad->SetLeftMargin(0.15);
    sdce->Draw("box1:text");
    dcan->cd(3);
    gPad->SetLeftMargin(0.15);
    sdde->Draw("box1:text");
    dcan->cd(2);
    gPad->SetLeftMargin(0.15);
    sce->Draw("box1:text");
    dcan->cd(4);
    gPad->SetLeftMargin(0.15);
    sbk->Draw("box1:text");
  }
}

