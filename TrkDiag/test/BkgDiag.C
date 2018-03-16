#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBox.h"

void BDiag(TTree* bdiag, const char* page="rho",bool train=false) {
  TString spage(page);
  TCut cluster("mvastat>0");
  if(train)cluster += TCut("ppdg==11&&nprimary/nshits>0.8");
  TCut con("pgen==2&&pmom>90");
  TCut bkg("pproc<20");
  // hit cuts
  TCut primary("_relation>=0");
  TCut norel("_relation<0");
  TCut stht("_stereo");
  TCut nstht("!_stereo");

  if(spage == "rho"){
    TH1F* rhocon = new TH1F("rhocon","Cluster Hit #rho;#rho (mm)",100,0,50);
    TH1F* rhobkg = new TH1F("rhobkg","Cluster Hit #rho;#rho (mm)",100,0,50);
    rhocon->SetLineColor(kRed);
    rhobkg->SetLineColor(kBlue);
    rhocon->SetStats(0);
    rhobkg->SetStats(0);
    bdiag->Project("rhocon","HitRho",cluster+con);
    bdiag->Project("rhobkg","HitRho",cluster+bkg);

    TH1F* radiuscon = new TH1F("radiuscon","Cluster Transverse Radius;#rho (mm)",100,350,700);
    TH1F* radiusbkg = new TH1F("radiusbkg","Cluster Transverse Radius;#rho (mm)",100,350,700);
    radiuscon->SetLineColor(kRed);
    radiusbkg->SetLineColor(kBlue);
    radiuscon->SetStats(0);
    radiusbkg->SetStats(0);
    bdiag->Project("radiuscon","ClusterRho",cluster+con);
    bdiag->Project("radiusbkg","ClusterRho",cluster+bkg);

    TH1F* rhoscon = new TH1F("rhoscon","Cluster Hit #rho RMS;#rho RMS (mm)",100,0,50);
    TH1F* rhosbkg = new TH1F("rhosbkg","Cluster Hit #rho RMS;#rho RMS(mm)",100,0,50);
    rhoscon->SetLineColor(kRed);
    rhosbkg->SetLineColor(kBlue);
    rhoscon->SetStats(0);
    rhosbkg->SetStats(0);
    bdiag->Project("rhoscon","HitRhoSpread",cluster+con);
    bdiag->Project("rhosbkg","HitRhoSpread",cluster+bkg);

    TH1F* timescon = new TH1F("timescon","Cluster time RMS;time RMS (ns)",100,0,20);
    TH1F* timesbkg = new TH1F("timesbkg","Cluster time RMS;time RMS(ns)",100,0,20);
    timescon->SetLineColor(kRed);
    timesbkg->SetLineColor(kBlue);
    timescon->SetStats(0);
    timesbkg->SetStats(0);
    bdiag->Project("timescon","TimeSpread",cluster+con);
    bdiag->Project("timesbkg","TimeSpread",cluster+bkg);

    rhocon->Scale(10);
    radiuscon->Scale(10);
    rhoscon->Scale(10);
    timescon->Scale(10);

    TLegend* rleg = new TLegend(0.6,0.7,0.9,0.9);
    rleg->AddEntry(rhobkg,"Background clusters","L");
    rleg->AddEntry(rhocon,"Conversion clusters (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",800,600);
    rhocan->Divide(2,2);
    rhocan->cd(1);
    rhobkg->Draw();
    rhocon->Draw("same");
    rleg->Draw();
    rhocan->cd(2);
    rhosbkg->Draw();
    rhoscon->Draw("same");
    rhocan->cd(3);
    timesbkg->Draw();
    timescon->Draw("same");
    rhocan->cd(4);
    radiusbkg->Draw();
    radiuscon->Draw("same");
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

    bdiag->Project("zmincon","ZMin",cluster+con);
    bdiag->Project("zminbkg","ZMin",cluster+bkg);
    bdiag->Project("zmaxcon","ZMax",cluster+con);
    bdiag->Project("zmaxbkg","ZMax",cluster+bkg);
    bdiag->Project("zgapcon","ZGap",cluster+con);
    bdiag->Project("zgapbkg","ZGap",cluster+bkg);
    bdiag->Project("zfraccon","ZGap/(ZMax-ZMin)",cluster+con);
    bdiag->Project("zfracbkg","ZGap/(ZMax-ZMin)",cluster+bkg);

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
    zfracbkg->Draw();
    zfraccon->Draw("same");

  } else if(spage == "planes"){

    TH1F* npmisscon = new TH1F("npmisscon","Missing Planes;N miss ",36,-0.5,36.5);
    TH1F* npmissbkg = new TH1F("npmissbkg","Missing Planes;N miss ",36,-0.5,36.5);
    TH1F* npcon = new TH1F("npcon","N Planes;N planes",37,-0.5,36.5);
    TH1F* npbkg = new TH1F("npbkg","N Planes;N planes",37,-0.5,36.5);
    TH1F* mfcon = new TH1F("mfcon","Plane Fraction;N planes",50,0.0,1.0);
    TH1F* mfbkg = new TH1F("mfbkg","Plane Fraction;N planes",50,0.0,1.0);

    npmisscon->SetLineColor(kRed);
    npmissbkg->SetLineColor(kBlue);
    npcon->SetLineColor(kRed);
    npbkg->SetLineColor(kBlue);
    mfcon->SetLineColor(kRed);
    mfbkg->SetLineColor(kBlue);

    npmisscon->SetStats(0);
    npmissbkg->SetStats(0);
    npcon->SetStats(0);
    npbkg->SetStats(0);
    mfcon->SetStats(0);
    mfbkg->SetStats(0);

    bdiag->Project("npmisscon","NExpectedPlanes-NPlanes",con+cluster);
    bdiag->Project("npmissbkg","NExpectedPlanes-NPlanes",bkg+cluster);
    bdiag->Project("npcon","NPlanes",con+cluster);
    bdiag->Project("npbkg","NPlanes",bkg+cluster);
    bdiag->Project("mfcon","PlaneFraction",con+cluster);
    bdiag->Project("mfbkg","PlaneFraction",bkg+cluster);

    npmisscon->Scale(10);
    npcon->Scale(10);
    mfcon->Scale(10);

    TLegend* stleg = new TLegend(0.2,0.7,0.8,0.9);
    stleg->AddEntry(npbkg,"Background clusters","L");
    stleg->AddEntry(npcon,"Conversion clusters (X10)","L");

    TCanvas* pcan = new TCanvas("pcan","pcan",1000,600);
    pcan->Divide(2,2);
    pcan->cd(1);
    npcon->Draw();
    npbkg->Draw("same");
    pcan->cd(2);
    npmissbkg->Draw();
    npmisscon->Draw("same");
    pcan->cd(3);
    mfcon->Draw();
    mfbkg->Draw("same");

  } else if(spage == "nhits"){

    TH1F* nshitscon = new TH1F("nshitscon","N Straw Hits",100,0.5,99.5);
    TH1F* nshitsbkg = new TH1F("nshitsbkg","N Straw Hits",100,0.5,99.5);
    TH1F* nchitscon = new TH1F("nchitscon","N ComboHits",100,0.5,99.5);
    TH1F* nchitsbkg = new TH1F("nchitsbkg","N ComboHits",100,0.5,99.5);
    TH1F* nphitscon = new TH1F("nphitscon","N Plane Hits",100,0.0,5.0);
    TH1F* nphitsbkg = new TH1F("nphitsbkg","N Plane Hits",100,0.0,5.0);
    TH1F* sfraccon = new TH1F("sfraccon","Stereo Fraction",100,0.0,1.0);
    TH1F* sfracbkg = new TH1F("sfracbkg","Stereo Fraction",100,0.0,1.0);

    nshitscon->SetLineColor(kRed);
    nshitsbkg->SetLineColor(kBlue);
    nchitscon->SetLineColor(kRed);
    nchitsbkg->SetLineColor(kBlue);
    nphitscon->SetLineColor(kRed);
    nphitsbkg->SetLineColor(kBlue);
    sfraccon->SetLineColor(kRed);
    sfracbkg->SetLineColor(kBlue);

    nshitscon->SetStats(0);
    nshitsbkg->SetStats(0);
    nphitscon->SetStats(0);
    nphitsbkg->SetStats(0);
    sfraccon->SetStats(0);
    sfracbkg->SetStats(0);
    
    bdiag->Project("nshitscon","nshits",con+cluster);
    bdiag->Project("nshitsbkg","nshits",bkg+cluster);
    bdiag->Project("nchitscon","nchits",con+cluster);
    bdiag->Project("nchitsbkg","nchits",bkg+cluster);
    bdiag->Project("nphitscon","NPlaneHits",con+cluster);
    bdiag->Project("nphitsbkg","NPlaneHits",bkg+cluster);
    bdiag->Project("sfraccon","StereoFraction",con+cluster);
    bdiag->Project("sfracbkg","StereoFraction",bkg+cluster);

    nshitscon->Scale(10);
    nchitscon->Scale(10);
    nphitscon->Scale(10);
    sfraccon->Scale(10);

    TLegend* nhleg = new TLegend(0.2,0.7,0.8,0.9);
    nhleg->AddEntry(nshitsbkg,"Background clusters","L");
    nhleg->AddEntry(nshitscon,"Conversion clusters (X10)","L");

    TCanvas* nhcan = new TCanvas("nhcan","nhcan",800,800);
    nhcan->Divide(2,2);
    nhcan->cd(1);
    nshitsbkg->Draw();
    nshitscon->Draw("same");
    nhleg->Draw();
    nhcan->cd(2);
    nchitsbkg->Draw();
    nchitscon->Draw("same");
    nhcan->cd(3);
    nphitsbkg->Draw();
    nphitscon->Draw("same");
    nhcan->cd(4);
    sfracbkg->Draw();
    sfraccon->Draw("same");
 } else if (spage=="mva") {
    TH1F* mvacon = new TH1F("mvacon","Background MVA output;MVA Output",200,-0.05,1.05);
    TH1F* mvabkg = new TH1F("mvabkg","Background MVA output;MVA Output",200,-0.05,1.05);
    mvacon->SetLineColor(kRed);
    mvabkg->SetLineColor(kBlue);
    mvacon->SetStats(0);
    mvabkg->SetStats(0);
    bdiag->Project("mvacon","mvaout",(con+cluster)*"nconv");
    bdiag->Project("mvabkg","mvaout",(bkg+cluster)*"nebkg");
    Double_t factor(50.0);
    mvacon->Scale(factor);
    TCanvas* can = new TCanvas("mvacan","MVA output",800,400);
    can->Divide(1,1);
    can->cd(1);
    gPad->SetLogy();
    mvabkg->Draw();
    mvacon->Draw("same");
    TBox* sel = new TBox(0.5,mvabkg->GetMinimum(),1.0,mvabkg->GetMaximum());
    sel->SetFillColor(kYellow);
    sel->SetFillStyle(3004);
    sel->Draw();
    TLegend* mcleg = new TLegend(0.5,0.7,0.8,0.9);
    mcleg->AddEntry(mvabkg,"Background electron","L");
    mcleg->AddEntry(mvacon,"Conversion electron (X50)","L");
    mcleg->AddEntry(sel,"Selection","F");
    mcleg->Draw();
 
  }  else if(spage=="hits") {
    TH1F* drhobkgp = new TH1F("drhobkgp","Bkg Hit #rho difference;#Delta #rho (mm)",100,-50.0,100.0);
    TH1F* drhobkgu = new TH1F("drhobkgu","Bkg Hit #rho difference;#Delta #rho (mm)",100,-50.0,100.0);
    TH1F* cdbkgp = new TH1F("cdbkgp","Bkg Hit cluster distance;distance",100,0,5.0);
    TH1F* cdbkgu = new TH1F("cdbkgu","Bkg Hit cluster distance;distance",100,0,5.0);
    TH1F* chibkgp = new TH1F("chibkgp","Bkg Hit cluster chi;chi",100,0,10.0);
    TH1F* chibkgu = new TH1F("chibkgu","Bkg Hit cluster chi;chi",100,0,10.0);
    TH1F* dtbkgp = new TH1F("dtbkgp","Bkg Hit time difference;#Delta t (nsec)",100,0,40);
    TH1F* dtbkgu = new TH1F("dtbkgu","Bkg Hit time difference;#Delta t (nsec)",100,0,40);
    drhobkgp->SetLineColor(kBlue);
    cdbkgp->SetLineColor(kBlue);
    chibkgp->SetLineColor(kBlue);
    dtbkgp->SetLineColor(kBlue);
    drhobkgu->SetLineColor(kRed);
    cdbkgu->SetLineColor(kRed);
    chibkgu->SetLineColor(kRed);
    dtbkgu->SetLineColor(kRed);

    TCut conhit("_mcgen==2&&_mcmom>90");

    bdiag->Project("drhobkgp","_rrho-HitRho",bkg+cluster+primary);
    bdiag->Project("cdbkgp","_gdist",bkg+cluster+primary);
    bdiag->Project("chibkgp","(_rrho-HitRho)/_rerr",bkg+cluster+primary);
    bdiag->Project("dtbkgp","_time-ctime",bkg+cluster+primary);
    bdiag->Project("drhobkgu","_rrho-HitRho",bkg+cluster+conhit);
    bdiag->Project("dtbkgu","_time-ctime",bkg+cluster+conhit);
    bdiag->Project("cdbkgu","_gdist",bkg+cluster+conhit);
    bdiag->Project("chibkgu","(_rrho-HitRho)/_rerr",bkg+cluster+conhit);

    drhobkgu->Scale(100);
    cdbkgu->Scale(100);
    chibkgu->Scale(100);
    dtbkgu->Scale(100);

    TCanvas* dhcan = new TCanvas("dhcan","Delta hits",800,800);
    dhcan->Divide(2,2);
    dhcan->cd(1);
    drhobkgp->Draw();
    drhobkgu->Draw("same");
    TLegend* hleg = new TLegend(0.5,0.7,0.9,0.9);
    hleg->AddEntry(drhobkgp,"Bkg Primary hit","L");
    hleg->AddEntry(drhobkgu,"Bgk unrelated hit (X100)","L");
    hleg->Draw();
    dhcan->cd(2);
    dtbkgp->Draw();
    dtbkgu->Draw("same");
    dhcan->cd(3);
    cdbkgp->Draw();
    cdbkgu->Draw("same");
    dhcan->cd(4);
    chibkgp->Draw();
    chibkgu->Draw("same");
  }    
}
