#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBox.h"

void BkgDiag(TTree* bdiag, const char* page="rho",bool train=false) {
  TString spage(page);
  TCut cluster("mvastat>0");
  if(train)cluster += TCut("ppdg==11&&nprimary/nhits>0.8");
  TCut con("pgen==2&&pmom>100");
  TCut bkg("pproc<20");
  // hit cuts
  TCut primary("_relation>=0");
  TCut norel("_relation<0");
  TCut stht("_stereo");
  TCut nstht("!_stereo");

  if(spage == "rho"){
    TH1F* rhocon = new TH1F("rhocon","Cluster #rho;#rho (mm)",100,0,35);
    TH1F* rhobkg = new TH1F("rhobkg","Cluster #rho;#rho (mm)",100,0,35);
    rhocon->SetLineColor(kRed);
    rhobkg->SetLineColor(kBlue);
    rhocon->SetStats(0);
    rhobkg->SetStats(0);
    bdiag->Project("rhocon","Rho",cluster+con);
    bdiag->Project("rhobkg","Rho",cluster+bkg);

    TH1F* rhoscon = new TH1F("rhoscon","Cluster #rho RMS;#rho RMS (mm)",100,0,25);
    TH1F* rhosbkg = new TH1F("rhosbkg","Cluster #rho RMS;#rho RMS(mm)",100,0,25);
    rhoscon->SetLineColor(kRed);
    rhosbkg->SetLineColor(kBlue);
    rhoscon->SetStats(0);
    rhosbkg->SetStats(0);
    bdiag->Project("rhoscon","RhoSpread",cluster+con);
    bdiag->Project("rhosbkg","RhoSpread",cluster+bkg);

    TH1F* timescon = new TH1F("timescon","Cluster time RMS;time RMS (ns)",100,0,30);
    TH1F* timesbkg = new TH1F("timesbkg","Cluster time RMS;time RMS(ns)",100,0,30);
    timescon->SetLineColor(kRed);
    timesbkg->SetLineColor(kBlue);
    timescon->SetStats(0);
    timesbkg->SetStats(0);
    bdiag->Project("timescon","TimeSpread",cluster+con);
    bdiag->Project("timesbkg","TimeSpread",cluster+bkg);

    rhocon->Scale(10);
    rhoscon->Scale(10);
    timescon->Scale(10);

    TLegend* rleg = new TLegend(0.2,0.7,0.6,0.9);
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
 } else if(spage == "z"){

    TH1F* zmincon = new TH1F("zmincon","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zminbkg = new TH1F("zminbkg","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zmaxcon = new TH1F("zmaxcon","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zmaxbkg = new TH1F("zmaxbkg","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zgapcon = new TH1F("zgapcon","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zgapbkg = new TH1F("zgapbkg","peak zgap;zgap (mm)",100,0,1600);

    zmincon->SetLineColor(kRed);
    zminbkg->SetLineColor(kBlue);
    zmaxcon->SetLineColor(kRed);
    zmaxbkg->SetLineColor(kBlue);
    zgapcon->SetLineColor(kRed);
    zgapbkg->SetLineColor(kBlue);

    zmincon->SetStats(0);
    zminbkg->SetStats(0);
    zmaxcon->SetStats(0);
    zmaxbkg->SetStats(0);
    zgapcon->SetStats(0);
    zgapbkg->SetStats(0);

    bdiag->Project("zmincon","MinZ",cluster+con);
    bdiag->Project("zminbkg","MinZ",cluster+bkg);
    bdiag->Project("zmaxcon","MaxZ",cluster+con);
    bdiag->Project("zmaxbkg","MaxZ",cluster+bkg);
    bdiag->Project("zgapcon","ZGap",cluster+con);
    bdiag->Project("zgapbkg","ZGap",cluster+bkg);

    zmincon->Scale(10);
    zmaxcon->Scale(10);
    zgapcon->Scale(10);

    TLegend* zleg = new TLegend(0.2,0.7,0.8,0.9);
    zleg->AddEntry(zminbkg,"Background clusters","L");
    zleg->AddEntry(zmincon,"Conversion clusters (X10)","L");

    TCanvas* zcan = new TCanvas("zcan","zcan",1200,800);
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

  } else if(spage == "planes"){

    TH1F* npmisscon = new TH1F("npmisscon","Missing Planes;N miss ",36,-0.5,36.5);
    TH1F* npmissbkg = new TH1F("npmissbkg","Missing Planes;N miss ",36,-0.5,36.5);
    TH1F* npcon = new TH1F("npcon","N Planes;N planes",37,-0.5,36.5);
    TH1F* npbkg = new TH1F("npbkg","N Planes;N planes",37,-0.5,36.5);
    TH1F* mfcon = new TH1F("mfcon","Missing Plane Fraction;N planes",50,-0.1,10.0);
    TH1F* mfbkg = new TH1F("mfbkg","Missing Plane Fraction;N planes",50,-0.1,10.0);

    npmisscon->SetLineColor(kRed);
    npmissbkg->SetLineColor(kBlue);
    npcon->SetLineColor(kRed);
    npbkg->SetLineColor(kBlue);

    npmisscon->SetStats(0);
    npmissbkg->SetStats(0);
    npcon->SetStats(0);
    npbkg->SetStats(0);

    bdiag->Project("npmisscon","MissingPlanes",con+cluster);
    bdiag->Project("npmissbkg","MissingPlanes",bkg+cluster);
    bdiag->Project("npcon","NPlanes",con+cluster);
    bdiag->Project("npbkg","NPlanes",bkg+cluster);

    npmisscon->Scale(10);
    npcon->Scale(10);

    TLegend* stleg = new TLegend(0.2,0.7,0.8,0.9);
    stleg->AddEntry(npbkg,"Background clusters","L");
    stleg->AddEntry(npcon,"Conversion clusters (X10)","L");

    TCanvas* pcan = new TCanvas("pcan","pcan",1200,800);
    pcan->Divide(1,2);
    pcan->cd(1);
    npbkg->Draw();
    npcon->Draw("same");
    pcan->cd(2);
    npmissbkg->Draw();
    npmisscon->Draw("same");

  } else if(spage == "nhits"){

    TH1F* nhitscon = new TH1F("nhitscon","N Hits",100,0.5,99.5);
    TH1F* nhitsbkg = new TH1F("nhitsbkg","N Hits",100,0.5,99.5);
    TH1F* nphitscon = new TH1F("nphitscon","N Plane Hits",100,0.0,5.0);
    TH1F* nphitsbkg = new TH1F("nphitsbkg","N Plane Hits",100,0.0,5.0);
    TH1F* sfraccon = new TH1F("sfraccon","Stereo Fraction",100,0.0,1.0);
    TH1F* sfracbkg = new TH1F("sfracbkg","Stereo Fraction",100,0.0,1.0);

    nhitscon->SetLineColor(kRed);
    nhitsbkg->SetLineColor(kBlue);
    nphitscon->SetLineColor(kRed);
    nphitsbkg->SetLineColor(kBlue);
    sfraccon->SetLineColor(kRed);
    sfracbkg->SetLineColor(kBlue);

    nhitscon->SetStats(0);
    nhitsbkg->SetStats(0);
    nphitscon->SetStats(0);
    nphitsbkg->SetStats(0);
    sfraccon->SetStats(0);
    sfracbkg->SetStats(0);
    
    bdiag->Project("nhitscon","NHits",con+cluster);
    bdiag->Project("nhitsbkg","NHits",bkg+cluster);
    bdiag->Project("nphitscon","NPlaneHits",con+cluster);
    bdiag->Project("nphitsbkg","NPlaneHits",bkg+cluster);
    bdiag->Project("sfraccon","StereoFraction",con+cluster);
    bdiag->Project("sfracbkg","StereoFraction",bkg+cluster);

    nhitscon->Scale(10);
    sfraccon->Scale(10);

    TLegend* nhleg = new TLegend(0.2,0.7,0.8,0.9);
    nhleg->AddEntry(nhitsbkg,"Background clusters","L");
    nhleg->AddEntry(nhitscon,"Conversion clusters (X10)","L");

    TCanvas* nhcan = new TCanvas("nhcan","nhcan",1200,800);
    nhcan->Divide(2,2);
    nhcan->cd(1);
    nhitsbkg->Draw();
    nhitscon->Draw("same");
    nhleg->Draw();
    nhcan->cd(2);
    nphitsbkg->Draw();
    nphitscon->Draw("same");
    nhcan->cd(3);
    sfracbkg->Draw();
    sfraccon->Draw("same");
 } else if (spage=="clustermva") {
    TH1F* stclmvacon = new TH1F("stclmvacon","Stereo Cluster MVA output;MVA Output",200,-0.2,1.2);
    TH1F* stclmvabkg = new TH1F("stclmvabkg","Stereo Cluster MVA output;MVA Output",200,-0.2,1.2);
    stclmvacon->SetLineColor(kRed);
    stclmvabkg->SetLineColor(kBlue);
    TH1F* nstclmvacon = new TH1F("nstclmvacon","Non-Stereo Cluster MVA output;MVA Output",200,-0.2,1.2);
    TH1F* nstclmvabkg = new TH1F("nstclmvabkg","Non-Stereo Cluster MVA output;MVA Output",200,-0.2,1.2);
    nstclmvacon->SetLineColor(kRed);
    nstclmvabkg->SetLineColor(kBlue);
    TCut stcl("ngdstereo/ngdhits>0.49");
    TCut nstcl("ngdstereo/ngdhits<0.49");
    bdiag->Project("stclmvacon","pmvaout",(con+cluster+stcl)*"nchits");
    bdiag->Project("stclmvabkg","pmvaout",(bkg+cluster+stcl)*"nchits");
    bdiag->Project("nstclmvacon","pmvaout",(con+cluster+nstcl)*"nchits");
    bdiag->Project("nstclmvabkg","pmvaout",(bkg+cluster+nstcl)*"nchits");
    Double_t factor(50.0);
    stclmvacon->Scale(factor);
    nstclmvacon->Scale(factor);
    TCanvas* can = new TCanvas("cmvacan","Cluster MVA output",800,400);
    can->Divide(2,1);
    can->cd(1);
    stclmvabkg->Draw();
    stclmvacon->Draw("same");
    TBox* ssel = new TBox(0.8,stclmvabkg->GetMinimum(),1.2,stclmvabkg->GetMaximum());
    ssel->SetFillColor(kYellow);
    ssel->SetFillStyle(3004);
    ssel->Draw();
    TLegend* mcleg = new TLegend(0.2,0.7,0.8,0.9);
    mcleg->AddEntry(stclmvabkg,"Background electron","L");
    mcleg->AddEntry(stclmvacon,"Conversion electron (X50)","L");
    mcleg->AddEntry(ssel,"Rejected","F");
    mcleg->Draw();
 
    can->cd(2);
    nstclmvabkg->Draw();
    nstclmvacon->Draw("same");
    TBox* nssel = new TBox(0.8,nstclmvabkg->GetMinimum(),1.2,nstclmvabkg->GetMaximum());
    nssel->SetFillColor(kYellow);
    nssel->SetFillStyle(3004);
    nssel->Draw();


  } else if(spage=="hitmva") {
    TH1F* sthtmvap = new TH1F("sthtmvap","Stereo Hit MVA output",200,-0.2,1.2);
    TH1F* sthtmvau = new TH1F("sthtmvau","Stereo Hit MVA output",200,-0.2,1.2);
    sthtmvau->SetLineColor(kRed);
    sthtmvap->SetLineColor(kBlue);
    TH1F* nsthtmvap = new TH1F("nsthtmvap","Non-Stereo Hit MVA output",200,-0.2,1.2);
    TH1F* nsthtmvau = new TH1F("nsthtmvau","Non-Stereo Hit MVA output",200,-0.2,1.2);
    nsthtmvap->SetLineColor(kBlue);
    nsthtmvau->SetLineColor(kRed);

    bdiag->Project("sthtmvap","_hgd",bkg+cluster+stht+primary);
    bdiag->Project("sthtmvau","_hgd",bkg+cluster+stht+norel);
    bdiag->Project("nsthtmvap","_hgd",bkg+cluster+nstht+primary);
    bdiag->Project("nsthtmvau","_hgd",bkg+cluster+nstht+norel);

    TLegend* mhleg = new TLegend(0.2,0.7,0.6,0.9);
    mhleg->AddEntry(sthtmvap,"#delta primary hits","L");
    mhleg->AddEntry(sthtmvau,"#delta unrelated (X10)","L");

    Double_t factor(10.0);
    sthtmvau->Scale(factor);
    nsthtmvau->Scale(factor);

    TCanvas* can = new TCanvas("hmvacan","Hit MVA output",800,400);
    can->Divide(2,1);
    can->cd(1);
    sthtmvap->Draw();
    sthtmvau->Draw("same");
    can->cd(2);
    nsthtmvap->Draw();
    nsthtmvau->Draw("same");
    mhleg->Draw();
  }  else if(spage=="hits") {

    TH1F* drhobkgp = new TH1F("drhobkgp","Bkg Hit #rho difference;#Delta #rho (mm)",100,-25,50.0);
    TH1F* cdbkgp = new TH1F("cdbkgp","Bkg Hit cluster distance;distance",100,0,2.0);
    TH1F* dtbkgp = new TH1F("dtbkgp","Bkg Hit time difference;#Delta t (nsec)",100,0,50);
     TH1F* drhobkgu = new TH1F("drhobkgu","Bkg Hit #rho difference;#Delta #rho (mm)",100,-25,50.0);
    TH1F* cdbkgu = new TH1F("cdbkgu","Bkg Hit cluster distance;distance",100,0,2.0);
    TH1F* dtbkgu = new TH1F("dtbkgu","Bkg Hit time difference;#Delta t (nsec)",100,0,50);
    drhobkgp->SetLineColor(kBlue);
    cdbkgp->SetLineColor(kBlue);
    dtbkgp->SetLineColor(kBlue);
    drhobkgu->SetLineColor(kRed);
    cdbkgu->SetLineColor(kRed);
    dtbkgu->SetLineColor(kRed);

    bdiag->Project("drhobkgp","_rho-Rho",bkg+cluster+primary);
    bdiag->Project("cdbkgp","_dist",bkg+cluster+primary);
    bdiag->Project("dtbkgp","_time-ctime",bkg+cluster+primary);
    bdiag->Project("drhobkgu","_rho-Rho",bkg+cluster+norel);
    bdiag->Project("dtbkgu","_time-ctime",bkg+cluster+norel);
    bdiag->Project("cdbkgu","_dist",bkg+cluster+norel);

    TCanvas* dhcan = new TCanvas("dhcan","Delta hits",1200,800);
    dhcan->Divide(2,2);
    dhcan->cd(1);
    drhobkgp->Draw();
    drhobkgu->Draw("same");
    TLegend* hleg = new TLegend(0.2,0.7,0.6,0.9);
    hleg->AddEntry(drhobkgp,"Bkg Primary hit","L");
    hleg->AddEntry(drhobkgu,"Bgk unrelated hit","L");
    hleg->Draw();
 
    dhcan->cd(2);
    drhobkgp->Draw();
    drhobkgu->Draw("same");
    dhcan->cd(3);
    dtbkgp->Draw();
    dtbkgu->Draw("same");
    dhcan->cd(4);
    cdbkgp->Draw();
    cdbkgu->Draw("same");
  }    
}
