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
  TCut pri("prel==0&&mmom.R()>100");
  TCut ebkg("prel<0&&abs(mpdg)==11&&mmom.R()<10.0");
  // hit cuts
  TCut main("_mrel>=0");
  TCut norel("_mrel<0");
  TCut prihit("_prel==0");

  if(spage == "rho"){

    TH1F* crhopri = new TH1F("crhopri","Cluster Transverse Radius;#rho (mm)",100,350,700);
    TH1F* crhobkg = new TH1F("crhobkg","Cluster Transverse Radius;#rho (mm)",100,350,700);
    crhopri->SetLineColor(kRed);
    crhobkg->SetLineColor(kBlue);
    crhopri->SetStats(0);
    crhobkg->SetStats(0);
    bdiag->Project("crhopri","cpos.Rho()",cluster+pri);
    bdiag->Project("crhobkg","cpos.Rho()",cluster+ebkg);

    crhopri->Scale(10);

    TLegend* rleg = new TLegend(0.6,0.7,0.9,0.9);
    rleg->AddEntry(crhobkg,"Background clusters","L");
    rleg->AddEntry(crhopri,"Conversion clusters (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",400,400);
    crhobkg->Draw();
    crhopri->Draw("same");
    rleg->Draw();

  } else if(spage == "rms"){
    TH1F* rhospri = new TH1F("rhospri","Cluster #rho RMS;#rho RMS (mm)",100,0,100);
    TH1F* rhosbkg = new TH1F("rhosbkg","Cluster #rho RMS;#rho RMS(mm)",100,0,100);
    rhospri->SetLineColor(kRed);
    rhosbkg->SetLineColor(kBlue);
    rhospri->SetStats(0);
    rhosbkg->SetStats(0);
    bdiag->Project("rhospri","sqrt(rmscposx^2+rmscposy^2)",cluster+pri);
    bdiag->Project("rhosbkg","sqrt(rmscposx^2+rmscposy^2)",cluster+ebkg);

    TH1F* timespri = new TH1F("timespri","Cluster time RMS;time RMS (ns)",100,0,50);
    TH1F* timesbkg = new TH1F("timesbkg","Cluster time RMS;time RMS(ns)",100,0,50);
    timespri->SetLineColor(kRed);
    timesbkg->SetLineColor(kBlue);
    timespri->SetStats(0);
    timesbkg->SetStats(0);
    bdiag->Project("timespri","rmsctime",cluster+pri);
    bdiag->Project("timesbkg","rmsctime",cluster+ebkg);

    rhospri->Scale(10);
    timespri->Scale(10);

    TLegend* rmsleg = new TLegend(0.6,0.7,0.9,0.9);
    rmsleg->AddEntry(rhosbkg,"Background clusters","L");
    rmsleg->AddEntry(rhospri,"Conversion clusters (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",800,400);
    rhocan->Divide(2,1);
    rhocan->cd(1);
    rhosbkg->Draw();
    rhospri->Draw("same");
    rmsleg->Draw();
    rhocan->cd(2);
    timesbkg->Draw();
    timespri->Draw("same");



  } else if(spage == "z"){

    TH2F* zminmaxpri = new TH2F("zminmaxpri","peak zmax vs zmin;zmin (mm);zmax (mm)",100,-1600,1600,100,-1600,1600);
    TH2F* zminmaxbkg = new TH2F("zminmaxbkg","peak zmax vs zmin;zmin (mm);zmax (mm)",100,-1600,1600,100,-1600,1600);
    TH1F* zgappri = new TH1F("zgappri","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zgapbkg = new TH1F("zgapbkg","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zfracpri = new TH1F("zfracpri","peak zfrac;zfrac (mm)",100,0,1.0);
    TH1F* zfracbkg = new TH1F("zfracbkg","peak zfrac;zfrac (mm)",100,0,1.0);

    zminmaxpri->SetLineColor(kRed);
    zminmaxbkg->SetLineColor(kBlue);
    zgappri->SetLineColor(kRed);
    zgapbkg->SetLineColor(kBlue);
    zfracpri->SetLineColor(kRed);
    zfracbkg->SetLineColor(kBlue);

    zminmaxpri->SetStats(0);
    zminmaxbkg->SetStats(0);
    zgappri->SetStats(0);
    zgapbkg->SetStats(0);
    zfracpri->SetStats(0);
    zfracbkg->SetStats(0);

    bdiag->Project("zminmaxpri","zmax:zmin",cluster+pri);
    bdiag->Project("zminmaxbkg","zmax:zmin",cluster+ebkg);
    bdiag->Project("zgappri","zgap",cluster+pri);
    bdiag->Project("zgapbkg","zgap",cluster+ebkg);
    bdiag->Project("zfracpri","zgap/(zmax-zmin)",cluster+pri);
    bdiag->Project("zfracbkg","zgap/(zmax-zmin)",cluster+ebkg);

    zminmaxpri->Scale(10);
    zgappri->Scale(10);
    zfracpri->Scale(10);

    TLegend* zleg = new TLegend(0.2,0.7,0.8,0.9);
    zleg->AddEntry(zgapbkg,"Background clusters","L");
    zleg->AddEntry(zgappri,"Conversion clusters (X10)","L");

    TCanvas* zcan = new TCanvas("zcan","zcan",1000,600);
    zcan->Divide(2,2);
    zcan->cd(1);
    gPad->SetLogz();
    zminmaxpri->Draw("colorz");
    zcan->cd(2);
    gPad->SetLogz();
    zminmaxbkg->Draw("colorz");
    zcan->cd(3);
    zgapbkg->Draw();
    zgappri->Draw("same");
    zleg->Draw();
    zcan->cd(4);
    zfracbkg->SetMaximum(max(zfracbkg->GetMaximum(),zfracpri->GetMaximum()));
    zfracbkg->Draw();
    zfracpri->Draw("same");

  } else if(spage == "planes"){

    TH1F* npmisspri = new TH1F("npmisspri","Missing Planes;N Missing",36,-0.5,36.5);
    TH1F* npmissbkg = new TH1F("npmissbkg","Missing Planes;N Missing",36,-0.5,36.5);
    TH1F* nppri = new TH1F("nppri","N Planes;N planes",37,-0.5,36.5);
    TH1F* npbkg = new TH1F("npbkg","N Planes;N planes",37,-0.5,36.5);
    TH2F* flppri = new TH2F("flppri","Last vs First Plane;First Plane;Last Plane",37,-0.5,36.5,37,-0.5,36.5);
    TH2F* flpbkg = new TH2F("flpbkg","Last vs First Plane;First Plane;Last Plane",37,-0.5,36.5,37,-0.5,36.5);
    TH1F* pgappri = new TH1F("pgappri","Plane Gap;planes",37,-0.5,36.5);
    TH1F* pgapbkg = new TH1F("pgapbkg","Plane Gap;planes",37,-0.5,36.5);
    TH1F* mfpri = new TH1F("mfpri","Plane Fraction;N planes/N Expected",50,0.0,1.001);
    TH1F* mfbkg = new TH1F("mfbkg","Plane Fraction;N planes/N Expected",50,0.0,1.001);

    npmisspri->SetLineColor(kRed);
    npmissbkg->SetLineColor(kBlue);
    nppri->SetLineColor(kRed);
    npbkg->SetLineColor(kBlue);
    flppri->SetLineColor(kRed);
    flpbkg->SetLineColor(kBlue);
    pgappri->SetLineColor(kRed);
    pgapbkg->SetLineColor(kBlue);
    mfpri->SetLineColor(kRed);
    mfbkg->SetLineColor(kBlue);

    npmisspri->SetStats(0);
    npmissbkg->SetStats(0);
    nppri->SetStats(0);
    npbkg->SetStats(0);
    flppri->SetStats(0);
    flpbkg->SetStats(0);
    pgappri->SetStats(0);
    pgapbkg->SetStats(0);
    mfpri->SetStats(0);
    mfbkg->SetStats(0);

    bdiag->Project("npmisspri","lp-fp+1-np",pri+cluster);
    bdiag->Project("npmissbkg","lp-fp+1-np",ebkg+cluster);
    bdiag->Project("nppri","np",pri+cluster);
    bdiag->Project("npbkg","np",ebkg+cluster);
    bdiag->Project("flppri","lp:fp",pri+cluster);
    bdiag->Project("flpbkg","lp:fp",ebkg+cluster);
    bdiag->Project("pgappri","pgap",pri+cluster);
    bdiag->Project("pgapbkg","pgap",ebkg+cluster);
    bdiag->Project("mfpri","np/(lp-fp+1)",pri+cluster);
    bdiag->Project("mfbkg","np/(lp-fp+1)",ebkg+cluster);

    npmisspri->Scale(10);
    nppri->Scale(10);
    flppri->Scale(10);
    pgappri->Scale(10);
    mfpri->Scale(10);

    TLegend* pleg = new TLegend(0.2,0.7,0.8,0.9);
    pleg->AddEntry(npbkg,"Background clusters","L");
    pleg->AddEntry(nppri,"Conversion clusters (X10)","L");

    TCanvas* pcan = new TCanvas("pcan","pcan",1200,500);
    pcan->Divide(3,2);
    pcan->cd(1);
    nppri->Draw();
    npbkg->Draw("same");
    pleg->Draw();
    pcan->cd(2);
    npmissbkg->Draw();
    npmisspri->Draw("same");
    pcan->cd(3);
    mfpri->Draw();
    mfbkg->Draw("same");
    pcan->cd(4);
    pgappri->Draw();
    pgapbkg->Draw("same");
    pcan->cd(5);
    gPad->SetLogz();
    flppri->Draw("colorz");
    pcan->cd(6);
    gPad->SetLogz();
    flpbkg->Draw("colorz");


  } else if(spage == "nhits"){

    TH1F* nshpri = new TH1F("nshpri","N Straw Hits",100,0.5,99.5);
    TH1F* nshbkg = new TH1F("nshbkg","N Straw Hits",100,0.5,99.5);
    TH1F* nchpri = new TH1F("nchpri","N ComboHits",60,0.5,59.5);
    TH1F* nchbkg = new TH1F("nchbkg","N ComboHits",60,0.5,59.5);

    nshpri->SetLineColor(kRed);
    nshbkg->SetLineColor(kBlue);
    nchpri->SetLineColor(kRed);
    nchbkg->SetLineColor(kBlue);

    nshpri->SetStats(0);
    nshbkg->SetStats(0);
    nchpri->SetStats(0);
    nchbkg->SetStats(0);

    bdiag->Project("nshpri","nsh",pri+cluster);
    bdiag->Project("nshbkg","nsh",ebkg+cluster);
    bdiag->Project("nchpri","nch",pri+cluster);
    bdiag->Project("nchbkg","nch",ebkg+cluster);

    nshpri->Scale(2);
    nchpri->Scale(2);

    TLegend* nhleg = new TLegend(0.2,0.7,0.8,0.9);
    nhleg->AddEntry(nshbkg,"Background clusters","L");
    nhleg->AddEntry(nshpri,"Conversion clusters (X2)","L");

    TCanvas* nhcan = new TCanvas("nhcan","nhcan",800,400);
    nhcan->Divide(2,1);
    nhcan->cd(1);
    nshbkg->Draw();
    nshpri->Draw("same");
    nhleg->Draw();
    nhcan->cd(2);
    nchbkg->Draw();
    nchpri->Draw("same");
  } else if (spage=="mva") {
    TH1F* mvapri = new TH1F("mvapri","Background MVA output;MVA Output",200,-0.05,1.05);
    TH1F* mvabkg = new TH1F("mvabkg","Background MVA output;MVA Output",200,-0.05,1.05);
    mvapri->SetLineColor(kRed);
    mvabkg->SetLineColor(kBlue);
    mvapri->SetStats(0);
    mvabkg->SetStats(0);


    bdiag->Project("mvapri","kQ",pri+cluster);
    bdiag->Project("mvabkg","kQ",ebkg+cluster);
    Double_t factor(10.0);
    mvapri->Scale(factor);
    TCanvas* mvacan = new TCanvas("mvacan","MVA output",800,400);
    mvacan->Divide(1,1);
    mvacan->cd(1);
    gPad->SetLogy();
    mvabkg->Draw();
    mvapri->Draw("same");
    TBox* sel = new TBox(0.5,mvabkg->GetMinimum(),1.0,mvabkg->GetMaximum());
    sel->SetFillColor(kYellow);
    sel->SetFillStyle(3004);
    sel->Draw();
    TLegend* mcleg = new TLegend(0.5,0.7,0.8,0.9);
    mcleg->AddEntry(mvabkg,"Background electron","L");
    mcleg->AddEntry(mvapri,"Conversion electron (X10)","L");
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
    TH1F* hrhopri = new TH1F("hrhopri","Cluster Hit #rho;#rho (mm)",100,350,700);
    TH1F* hrhobkg = new TH1F("hrhobkg","Cluster Hit #rho;#rho (mm)",100,350,700);

    hrhopri->SetLineColor(kRed);
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


    hrhopri->SetStats(0);
    hrhobkg->SetStats(0);
    bdiag->Project("hrhopri","bkghinfo._pos.Rho()",cluster+pri);
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
