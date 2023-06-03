#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBox.h"

void BDiag(TTree* bdiag, const char* page="rho") {
  TString spage(page);
  TCut cluster("nch>=3");
  TCut pri("prel==0&&mmom>100");
  TCut ebkg("prel<0&&abs(mpdg)==11&&mmom<10.0");
  // hit cuts
  TCut main("_mrel>=0");
  TCut norel("_mrel<0");
  TCut prihit("_prel==0");

  if(spage == "rho"){

    TH1F* crhocon = new TH1F("crhocon","Cluster Transverse Radius;#rho (mm)",100,350,700);
    TH1F* crhobkg = new TH1F("crhobkg","Cluster Transverse Radius;#rho (mm)",100,350,700);
    crhocon->SetLineColor(kRed);
    crhobkg->SetLineColor(kBlue);
    crhocon->SetStats(0);
    crhobkg->SetStats(0);
    bdiag->Project("crhocon","cpos.Rho()",cluster+pri);
    bdiag->Project("crhobkg","cpos.Rho()",cluster+ebkg);

    crhocon->Scale(10);

    TLegend* rleg = new TLegend(0.6,0.7,0.9,0.9);
    rleg->AddEntry(crhobkg,"Background clusters","L");
    rleg->AddEntry(crhocon,"Conversion clusters (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",400,400);
    crhobkg->Draw();
    crhocon->Draw("same");
    rleg->Draw();

  } else if(spage == "rms"){
    TH1F* rhoscon = new TH1F("rhoscon","Cluster #rho RMS;#rho RMS (mm)",100,0,100);
    TH1F* rhosbkg = new TH1F("rhosbkg","Cluster #rho RMS;#rho RMS(mm)",100,0,100);
    rhoscon->SetLineColor(kRed);
    rhosbkg->SetLineColor(kBlue);
    rhoscon->SetStats(0);
    rhosbkg->SetStats(0);
    bdiag->Project("rhoscon","sqrt(rmscposx^2+rmscposy^2)",cluster+pri);
    bdiag->Project("rhosbkg","sqrt(rmscposx^2+rmscposy^2)",cluster+ebkg);

    TH1F* timescon = new TH1F("timescon","Cluster time RMS;time RMS (ns)",100,0,50);
    TH1F* timesbkg = new TH1F("timesbkg","Cluster time RMS;time RMS(ns)",100,0,50);
    timescon->SetLineColor(kRed);
    timesbkg->SetLineColor(kBlue);
    timescon->SetStats(0);
    timesbkg->SetStats(0);
    bdiag->Project("timescon","rmsctime",cluster+pri);
    bdiag->Project("timesbkg","rmsctime",cluster+ebkg);

    rhoscon->Scale(10);
    timescon->Scale(10);

    TLegend* rmsleg = new TLegend(0.6,0.7,0.9,0.9);
    rmsleg->AddEntry(rhosbkg,"Background clusters","L");
    rmsleg->AddEntry(rhoscon,"Conversion clusters (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",800,400);
    rhocan->Divide(2,1);
    rhocan->cd(1);
    rhosbkg->Draw();
    rhoscon->Draw("same");
    rmsleg->Draw();
    rhocan->cd(2);
    timesbkg->Draw();
    timescon->Draw("same");



  } else if(spage == "z"){

    TH1F* zmincon = new TH1F("zmincon","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zminbkg = new TH1F("zminbkg","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zmaxcon = new TH1F("zmaxcon","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zmaxbkg = new TH1F("zmaxbkg","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zgapcon = new TH1F("zgapcon","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zgapbkg = new TH1F("zgapbkg","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zfraccon = new TH1F("zfraccon","peak zfrac;zfrac (mm)",100,0,1.0);
    TH1F* zfracbkg = new TH1F("zfracbkg","peak zfrac;zfrac (mm)",100,0,1.0);

    zmincon->SetLineColor(kRed);
    zminbkg->SetLineColor(kBlue);
    zmaxcon->SetLineColor(kRed);
    zmaxbkg->SetLineColor(kBlue);
    zgapcon->SetLineColor(kRed);
    zgapbkg->SetLineColor(kBlue);
    zfraccon->SetLineColor(kRed);
    zfracbkg->SetLineColor(kBlue);

    zmincon->SetStats(0);
    zminbkg->SetStats(0);
    zmaxcon->SetStats(0);
    zmaxbkg->SetStats(0);
    zgapcon->SetStats(0);
    zgapbkg->SetStats(0);
    zfraccon->SetStats(0);
    zfracbkg->SetStats(0);

    bdiag->Project("zmincon","zmin",cluster+pri);
    bdiag->Project("zminbkg","zmin",cluster+ebkg);
    bdiag->Project("zmaxcon","zmax",cluster+pri);
    bdiag->Project("zmaxbkg","zmax",cluster+ebkg);
    bdiag->Project("zgapcon","zgap",cluster+pri);
    bdiag->Project("zgapbkg","zgap",cluster+ebkg);
    bdiag->Project("zfraccon","zgap/(zmax-zmin)",cluster+pri);
    bdiag->Project("zfracbkg","zgap/(zmax-zmin)",cluster+ebkg);

    zmincon->Scale(10);
    zmaxcon->Scale(10);
    zgapcon->Scale(10);
    zfraccon->Scale(10);

    TLegend* zleg = new TLegend(0.2,0.7,0.8,0.9);
    zleg->AddEntry(zminbkg,"Background clusters","L");
    zleg->AddEntry(zmincon,"Conversion clusters (X10)","L");

    TCanvas* zcan = new TCanvas("zcan","zcan",1000,600);
    zcan->Divide(2,2);
    zcan->cd(1);
    zminbkg->Draw();
    zmincon->Draw("same");
    zleg->Draw();
    zcan->cd(2);
    zmaxbkg->Draw();
    zmaxcon->Draw("same");
    zcan->cd(3);
    zgapbkg->Draw();
    zgapcon->Draw("same");
    zcan->cd(4);
    zfracbkg->SetMaximum(max(zfracbkg->GetMaximum(),zfraccon->GetMaximum()));
    zfracbkg->Draw();
    zfraccon->Draw("same");

  } else if(spage == "planes"){

    TH1F* npmisscon = new TH1F("npmisscon","Missing Planes;N Missing",36,-0.5,36.5);
    TH1F* npmissbkg = new TH1F("npmissbkg","Missing Planes;N Missing",36,-0.5,36.5);
    TH1F* npcon = new TH1F("npcon","N Planes;N planes",37,-0.5,36.5);
    TH1F* npbkg = new TH1F("npbkg","N Planes;N planes",37,-0.5,36.5);
    TH1F* fpcon = new TH1F("fpcon","First Plane;planes",37,-0.5,36.5);
    TH1F* fpbkg = new TH1F("fpbkg","First Plane;planes",37,-0.5,36.5);
    TH1F* lpcon = new TH1F("lpcon","Last Plane;planes",37,-0.5,36.5);
    TH1F* lpbkg = new TH1F("lpbkg","Last Plane;planes",37,-0.5,36.5);
    TH1F* pgapcon = new TH1F("pgapcon","Plane Gap;planes",37,-0.5,36.5);
    TH1F* pgapbkg = new TH1F("pgapbkg","Plane Gap;planes",37,-0.5,36.5);
    TH1F* pspancon = new TH1F("pspancon","Plane Span;planes",37,-0.5,36.5);
    TH1F* pspanbkg = new TH1F("pspanbkg","Plane Span;planes",37,-0.5,36.5);
    TH1F* mfcon = new TH1F("mfcon","Plane Fraction;N planes/N Expected",50,0.0,1.0);
    TH1F* mfbkg = new TH1F("mfbkg","Plane Fraction;N planes/N Expected",50,0.0,1.0);

    npmisscon->SetLineColor(kRed);
    npmissbkg->SetLineColor(kBlue);
    npcon->SetLineColor(kRed);
    npbkg->SetLineColor(kBlue);
    fpcon->SetLineColor(kRed);
    fpbkg->SetLineColor(kBlue);
    lpcon->SetLineColor(kRed);
    lpbkg->SetLineColor(kBlue);
    pgapcon->SetLineColor(kRed);
    pgapbkg->SetLineColor(kBlue);
    pspancon->SetLineColor(kRed);
    pspanbkg->SetLineColor(kBlue);
    mfcon->SetLineColor(kRed);
    mfbkg->SetLineColor(kBlue);

    npmisscon->SetStats(0);
    npmissbkg->SetStats(0);
    npcon->SetStats(0);
    npbkg->SetStats(0);
    fpcon->SetStats(0);
    fpbkg->SetStats(0);
    lpcon->SetStats(0);
    lpbkg->SetStats(0);
    pspancon->SetStats(0);
    pspanbkg->SetStats(0);
    pgapcon->SetStats(0);
    pgapbkg->SetStats(0);
    mfcon->SetStats(0);
    mfbkg->SetStats(0);

    bdiag->Project("npmisscon","lp-fp+1-np",pri+cluster);
    bdiag->Project("npmissbkg","lp-fp+1-np",ebkg+cluster);
    bdiag->Project("npcon","np",pri+cluster);
    bdiag->Project("npbkg","np",ebkg+cluster);
    bdiag->Project("fpcon","fp",pri+cluster);
    bdiag->Project("fpbkg","fp",ebkg+cluster);
    bdiag->Project("lpcon","lp",pri+cluster);
    bdiag->Project("lpbkg","lp",ebkg+cluster);
    bdiag->Project("pgapcon","pgap",pri+cluster);
    bdiag->Project("pgapbkg","pgap",ebkg+cluster);
    bdiag->Project("pspancon","lp-fp+1",pri+cluster);
    bdiag->Project("pspanbkg","lp-fp+1",ebkg+cluster);
    bdiag->Project("mfcon","np/(lp-fp+1)",pri+cluster);
    bdiag->Project("mfbkg","np/(lp-fp+1)",ebkg+cluster);

    npmisscon->Scale(10);
    npcon->Scale(10);
    fpcon->Scale(10);
    lpcon->Scale(10);
    pgapcon->Scale(10);
    pspancon->Scale(10);
    mfcon->Scale(10);

    TLegend* pleg = new TLegend(0.2,0.7,0.8,0.9);
    pleg->AddEntry(npbkg,"Background clusters","L");
    pleg->AddEntry(npcon,"Conversion clusters (X10)","L");

    TCanvas* pcan = new TCanvas("pcan","pcan",1200,500);
    pcan->Divide(4,2);
    pcan->cd(1);
    npcon->Draw();
    npbkg->Draw("same");
    pleg->Draw();
    pcan->cd(2);
    npmissbkg->Draw();
    npmisscon->Draw("same");
    pcan->cd(3);
    mfcon->Draw();
    mfbkg->Draw("same");
    pcan->cd(4);
    pgapcon->Draw();
    pgapbkg->Draw("same");
    pcan->cd(5);
    fpcon->Draw();
    fpbkg->Draw("same");
    pcan->cd(6);
    lpcon->Draw();
    lpbkg->Draw("same");
    pcan->cd(7);
    pgapcon->Draw();
    pgapbkg->Draw("same");
    pcan->cd(8);
    pspancon->Draw();
    pspanbkg->Draw("same");


  } else if(spage == "nhits"){

    TH1F* nshcon = new TH1F("nshcon","N Straw Hits",100,0.5,99.5);
    TH1F* nshbkg = new TH1F("nshbkg","N Straw Hits",100,0.5,99.5);
    TH1F* nchcon = new TH1F("nchcon","N ComboHits",60,0.5,59.5);
    TH1F* nchbkg = new TH1F("nchbkg","N ComboHits",60,0.5,59.5);

    nshcon->SetLineColor(kRed);
    nshbkg->SetLineColor(kBlue);
    nchcon->SetLineColor(kRed);
    nchbkg->SetLineColor(kBlue);

    nshcon->SetStats(0);
    nshbkg->SetStats(0);
    nchcon->SetStats(0);
    nchbkg->SetStats(0);

    bdiag->Project("nshcon","nsh",pri+cluster);
    bdiag->Project("nshbkg","nsh",ebkg+cluster);
    bdiag->Project("nchcon","nch",pri+cluster);
    bdiag->Project("nchbkg","nch",ebkg+cluster);

    nshcon->Scale(2);
    nchcon->Scale(2);

    TLegend* nhleg = new TLegend(0.2,0.7,0.8,0.9);
    nhleg->AddEntry(nshbkg,"Background clusters","L");
    nhleg->AddEntry(nshcon,"Conversion clusters (X2)","L");

    TCanvas* nhcan = new TCanvas("nhcan","nhcan",800,400);
    nhcan->Divide(2,1);
    nhcan->cd(1);
    nshbkg->Draw();
    nshcon->Draw("same");
    nhleg->Draw();
    nhcan->cd(2);
    nchbkg->Draw();
    nchcon->Draw("same");
  } else if (spage=="mva") {
    TH1F* mvacon = new TH1F("mvacon","Background MVA output;MVA Output",200,-0.05,1.05);
    TH1F* mvabkg = new TH1F("mvabkg","Background MVA output;MVA Output",200,-0.05,1.05);
    mvacon->SetLineColor(kRed);
    mvabkg->SetLineColor(kBlue);
    mvacon->SetStats(0);
    mvabkg->SetStats(0);


    bdiag->Project("mvacon","mvaout",pri+cluster);
    bdiag->Project("mvabkg","mvaout",ebkg+cluster);
    Double_t factor(10.0);
    mvacon->Scale(factor);
    TCanvas* mvacan = new TCanvas("mvacan","MVA output",800,400);
    mvacan->Divide(1,1);
    mvacan->cd(1);
    gPad->SetLogy();
    mvabkg->Draw();
    mvacon->Draw("same");
    TBox* sel = new TBox(0.5,mvabkg->GetMinimum(),1.0,mvabkg->GetMaximum());
    sel->SetFillColor(kYellow);
    sel->SetFillStyle(3004);
    sel->Draw();
    TLegend* mcleg = new TLegend(0.5,0.7,0.8,0.9);
    mcleg->AddEntry(mvabkg,"Background electron","L");
    mcleg->AddEntry(mvacon,"Conversion electron (X10)","L");
    mcleg->AddEntry(sel,"Selection","F");
    mcleg->Draw();

  }  else if(spage=="bhits") {
    TH1F* drhobkgp = new TH1F("drhobkgp","Bkg Hit #rho difference;#Delta #rho (mm)",100,-100,100);
    TH1F* drhobkgu = new TH1F("drhobkgu","Bkg Hit #rho difference;#Delta #rho (mm)",100,-100,100);
    TH1F* cdbkgp = new TH1F("cdbkgp","Bkg Hit cluster distance;distance",100,0,120);
    TH1F* cdbkgu = new TH1F("cdbkgu","Bkg Hit cluster distance;distance",100,0,120);
    TH1F* chibkgp = new TH1F("chibkgp","Bkg Hit cluster chi;chi",100,0,20.0);
    TH1F* chibkgu = new TH1F("chibkgu","Bkg Hit cluster chi;chi",100,0,20.0);
    TH1F* dtbkgp = new TH1F("dtbkgp","Bkg Hit time difference;#Delta t (nsec)",100,0,40);
    TH1F* dtbkgu = new TH1F("dtbkgu","Bkg Hit time difference;#Delta t (nsec)",100,0,40);
    TH1F* hrhocon = new TH1F("hrhocon","Cluster Hit #rho;#rho (mm)",100,350,700);
    TH1F* hrhobkg = new TH1F("hrhobkg","Cluster Hit #rho;#rho (mm)",100,350,700);

    hrhocon->SetLineColor(kRed);
    hrhobkg->SetLineColor(kBlue);
    drhobkgp->SetLineColor(kBlue);
    cdbkgp->SetLineColor(kBlue);
    chibkgp->SetLineColor(kBlue);
    dtbkgp->SetLineColor(kBlue);
    drhobkgu->SetLineColor(kRed);
    cdbkgu->SetLineColor(kRed);
    chibkgu->SetLineColor(kRed);
    dtbkgu->SetLineColor(kRed);

    bdiag->Project("drhobkgp","bkghinfo._pos.Rho()-cpos.Rho()",ebkg+cluster+main);
    bdiag->Project("cdbkgp","sqrt((bkghinfo._pos.X()-cpos.X())^2+(bkghinfo._pos.Y()-cpos.Y())^2)",ebkg+cluster+main);
    bdiag->Project("chibkgp","(bkghinfo._pos.Rho()-cpos.Rho())/_rerr",ebkg+cluster+main);
    bdiag->Project("dtbkgp","_time-ctime",ebkg+cluster+main);
    bdiag->Project("drhobkgu","bkghinfo._pos.Rho()-cpos.Rho()",ebkg+cluster+prihit);
    bdiag->Project("dtbkgu","_time-ctime",ebkg+cluster+prihit);
    bdiag->Project("cdbkgu","sqrt((bkghinfo._pos.X()-cpos.X())^2+(bkghinfo._pos.Y()-cpos.Y())^2)",ebkg+cluster+prihit);
    bdiag->Project("chibkgu","(bkghinfo._pos.Rho()-cpos.Rho())/_rerr",ebkg+cluster+prihit);

    drhobkgu->Scale(100);
    cdbkgu->Scale(100);
    chibkgu->Scale(100);
    dtbkgu->Scale(100);


    hrhocon->SetStats(0);
    hrhobkg->SetStats(0);
    bdiag->Project("hrhocon","bkghinfo._pos.Rho()",cluster+pri);
    bdiag->Project("hrhobkg","bkghinfo._pos.Rho()",cluster+ebkg);

    TCanvas* dhcan = new TCanvas("dhcan","Delta hits",800,800);
    dhcan->Divide(2,2);
    dhcan->cd(1);
    drhobkgp->Draw();
    drhobkgu->Draw("same");
    TLegend* hleg = new TLegend(0.5,0.7,0.9,0.9);
    hleg->AddEntry(drhobkgp,"Bkg Cluster Primary hit","L");
    hleg->AddEntry(drhobkgu,"Bgk Cluster Ce hit (X100)","L");
    hleg->Draw();
    dhcan->cd(2);
    dtbkgp->Draw();
    dtbkgu->Draw("same");
    dhcan->cd(3);
    cdbkgp->Draw();
    cdbkgu->Draw("same");
    dhcan->cd(4);
    gPad->SetLogy();
    chibkgp->Draw();
    chibkgu->Draw("same");
  }
}
