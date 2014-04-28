#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
void Delta(TTree* ddiag, const char* page="rho",bool train=false) {
  TString spage(page);
  TCut cluster("nchits>4&&ns>1");
  if(train)cluster += TCut("nprimary/nchits>0.8");
  TCut con("pgen==2&&mcmom>100");
  TCut bkg("pproc<20");
  TCut primary("_relation>=0");
  TCut norel("_relation<0");
  TCut stht("_stereo");
  TCut nstht("!_stereo");

  if(spage == "rho"){
    TH1F* mrhocon = new TH1F("mrhocon","Cluster #rho;#rho (mm)",100,330,800);
    TH1F* mrhobkg = new TH1F("mrhobkg","Cluster #rho;#rho (mm)",100,330,800);
    mrhocon->SetLineColor(kRed);
    mrhobkg->SetLineColor(kBlue);

    mrhocon->SetStats(0);
    mrhobkg->SetStats(0);

    ddiag->Project("mrhocon","rmean",cluster+con);
    ddiag->Project("mrhobkg","rmean",cluster+bkg);

    mrhocon->Scale(10);

    TLegend* rleg = new TLegend(0.2,0.7,0.6,0.9);
    rleg->AddEntry(mrhobkg,"Background peaks","L");
    rleg->AddEntry(mrhocon,"Conversion peaks (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",800,600);
    rhocan->Divide(1,1);
    rhocan->cd(1);
    mrhobkg->Draw();
    mrhocon->Draw("same");
    rleg->Draw();

  } else if(spage == "spread"){

    TH1F* srhocon = new TH1F("srhocon","Sigma of hit #rho distribution;#sigma #rho (mm)",100,0,50);
    TH1F* srhobkg = new TH1F("srhobkg","Sigma of hit #rho distribution;#sigma #rho (mm)",100,0,50);
    TH1F* sphicon = new TH1F("sphicon","Sigma of hit #phi distribution;#sigma #phi",100,0,0.20);
    TH1F* sphibkg = new TH1F("sphibkg","Sigma of hit #phi distribution;#sigma #phi",100,0,0.20);

    srhocon->SetLineColor(kRed);
    srhobkg->SetLineColor(kBlue);
    sphicon->SetLineColor(kRed);
    sphibkg->SetLineColor(kBlue);

    srhocon->SetStats(0);
    srhobkg->SetStats(0);
    sphicon->SetStats(0);
    sphibkg->SetStats(0);

    ddiag->Project("srhocon","srho",cluster+con);
    ddiag->Project("srhobkg","srho",cluster+bkg);
    ddiag->Project("sphicon","sphi",cluster+con);
    ddiag->Project("sphibkg","sphi",cluster+bkg);

    srhocon->Scale(10);
    sphicon->Scale(10);

    TLegend* sleg = new TLegend(0.2,0.7,0.8,0.9);
    sleg->AddEntry(srhobkg,"Background peaks","L");
    sleg->AddEntry(srhocon,"Conversion peaks (X10)","L");

    TCanvas* scan = new TCanvas("scan","scan",800,600);
    scan->Divide(2,1);
    scan->cd(1);
    srhobkg->Draw();
    srhocon->Draw("same");
    sleg->Draw();
    scan->cd(2);
    sphibkg->Draw();
    sphicon->Draw("same");
  
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

    ddiag->Project("zmincon","zmin",cluster+con);
    ddiag->Project("zminbkg","zmin",cluster+bkg);
    ddiag->Project("zmaxcon","zmax",cluster+con);
    ddiag->Project("zmaxbkg","zmax",cluster+bkg);
    ddiag->Project("zgapcon","zgap",cluster+con);
    ddiag->Project("zgapbkg","zgap",cluster+bkg);

    zmincon->Scale(10);
    zmaxcon->Scale(10);
    zgapcon->Scale(10);

    TLegend* zleg = new TLegend(0.2,0.7,0.8,0.9);
    zleg->AddEntry(zminbkg,"Background peaks","L");
    zleg->AddEntry(zmincon,"Conversion peaks (X10)","L");

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

  } else if(spage == "stations"){

    TH1F* smincon = new TH1F("smincon","peak smin;smin ",23,-0.5,22.5);
    TH1F* sminbkg = new TH1F("sminbkg","peak smin;smin ",23,-0.5,22.5);
    TH1F* smaxcon = new TH1F("smaxcon","peak smax;smax ",23,-0.5,22.5);
    TH1F* smaxbkg = new TH1F("smaxbkg","peak smax;smax ",23,-0.5,22.5);
    TH1F* nsmisscon = new TH1F("nsmisscon","peak nsmiss;nsmiss ",23,-0.5,22.5);
    TH1F* nsmissbkg = new TH1F("nsmissbkg","peak nsmiss;nsmiss ",23,-0.5,22.5);
    TH1F* nscon = new TH1F("nscon","peak ns;ns (mm)",23,-0.5,22.5);
    TH1F* nsbkg = new TH1F("nsbkg","peak ns;ns (mm)",23,-0.5,22.5);

    smincon->SetLineColor(kRed);
    sminbkg->SetLineColor(kBlue);
    smaxcon->SetLineColor(kRed);
    smaxbkg->SetLineColor(kBlue);
    nsmisscon->SetLineColor(kRed);
    nsmissbkg->SetLineColor(kBlue);
    nscon->SetLineColor(kRed);
    nsbkg->SetLineColor(kBlue);

    smincon->SetStats(0);
    sminbkg->SetStats(0);
    smaxcon->SetStats(0);
    smaxbkg->SetStats(0);
    nsmisscon->SetStats(0);
    nsmissbkg->SetStats(0);
    nscon->SetStats(0);
    nsbkg->SetStats(0);

    ddiag->Project("smincon","smin",con+cluster);
    ddiag->Project("sminbkg","smin",bkg+cluster);
    ddiag->Project("smaxcon","smax",con+cluster);
    ddiag->Project("smaxbkg","smax",bkg+cluster);
    ddiag->Project("nsmisscon","nsmiss",con+cluster);
    ddiag->Project("nsmissbkg","nsmiss",bkg+cluster);
    ddiag->Project("nscon","ns",con+cluster);
    ddiag->Project("nsbkg","ns",bkg+cluster);

    smincon->Scale(10);
    smaxcon->Scale(10);
    nsmisscon->Scale(10);
    nscon->Scale(10);

    TLegend* stleg = new TLegend(0.2,0.7,0.8,0.9);
    stleg->AddEntry(sminbkg,"Background peaks","L");
    stleg->AddEntry(smincon,"Conversion peaks (X10)","L");

    TCanvas* stcan = new TCanvas("stcan","stcan",1200,800);
    stcan->Divide(2,2);
    stcan->cd(1);
    sminbkg->Draw();
    smincon->Draw("same");
    stleg->Draw();
    stcan->cd(2);
    smaxbkg->Draw();
    smaxcon->Draw("same");
    stcan->cd(3);
    nsmissbkg->Draw();
    nsmisscon->Draw("same");
    stcan->cd(4);
    nsbkg->Draw();
    nscon->Draw("same");

  } else if(spage == "nhits"){

    TH1F* ngdhitscon = new TH1F("ngdhitscon","ngdhits",100,0,200);
    TH1F* ngdhitsbkg = new TH1F("ngdhitsbkg","ngdhits",100,0,200);
    TH1F* nchitscon = new TH1F("nchitscon","nchits",100,0,200);
    TH1F* nchitsbkg = new TH1F("nchitsbkg","nchits",100,0,200);
    TH1F* cwratiocon = new TH1F("cwratiocon","Ratio core/wide",100,-0.01,1.01);
    TH1F* cwratiobkg = new TH1F("cwratiobkg","Ratio core/wide",100,-0.01,1.01);

    ngdhitscon->SetLineColor(kRed);
    ngdhitsbkg->SetLineColor(kBlue);
    nchitscon->SetLineColor(kRed);
    nchitsbkg->SetLineColor(kBlue);
    cwratiocon->SetLineColor(kRed);
    cwratiobkg->SetLineColor(kBlue);

    ngdhitscon->SetStats(0);
    ngdhitsbkg->SetStats(0);
    nchitscon->SetStats(0);
    nchitsbkg->SetStats(0);
    cwratiocon->SetStats(0);
    cwratiobkg->SetStats(0);

    ddiag->Project("ngdhitscon","ngdhits",con+cluster);
    ddiag->Project("ngdhitsbkg","ngdhits",bkg+cluster);
    ddiag->Project("nchitscon","nchits",con+cluster);
    ddiag->Project("nchitsbkg","nchits",bkg+cluster);
    ddiag->Project("cwratiocon","ngdhits/nchits",con+cluster);
    ddiag->Project("cwratiobkg","ngdhits/nchits",bkg+cluster);

    ngdhitscon->Scale(10);
    nchitscon->Scale(10);
    cwratiocon->Scale(10);

    TLegend* nhleg = new TLegend(0.2,0.7,0.8,0.9);
    nhleg->AddEntry(ngdhitsbkg,"Background peaks","L");
    nhleg->AddEntry(ngdhitscon,"Conversion peaks (X10)","L");

    TCanvas* nhcan = new TCanvas("nhcan","nhcan",1200,800);
    nhcan->Divide(2,2);
    nhcan->cd(1);
    ngdhitsbkg->Draw();
    ngdhitscon->Draw("same");
    nhleg->Draw();
    nhcan->cd(2);
    nchitsbkg->Draw();
    nchitscon->Draw("same");
    nhcan->cd(3);
    cwratiobkg->Draw();
    cwratiocon->Draw("same");
  } else if(spage=="size"){
    TH2F* csize = new TH2F("csize","#delta-ray Cluster Size;#Delta x(mm);#Delta y(mm)",50,-30,30,50,-30,30);
    TH1F* tsize = new TH1F("tsize","#delta-ray Time Width;#Delta t(nsec)",100,-80,80);
    csize->SetStats(0);
    tsize->SetStats(0);
    ddiag->Project("csize","rmed*sin(pmed)-_mcpos.dy:rmed*cos(pmed)-_mcpos.dx",bkg+cluster);
    ddiag->Project("tsize","_time-tmed",bkg+cluster);
    TCanvas* scan = new TCanvas("scan","size",800,400);
    scan->Divide(2,1);
    scan->cd(1);
    csize->Draw("colorz");
    scan->cd(2);
    tsize->Draw();
  } else if (spage=="clustermva") {
    TH1F* stclmvacon = new TH1F("stclmvacon","Stereo Cluster MVA output",200,-0.2,1.2);
    TH1F* stclmvabkg = new TH1F("stclmvabkg","Stereo Cluster MVA output",200,-0.2,1.2);
    stclmvacon->SetLineColor(kRed);
    stclmvabkg->SetLineColor(kBlue);
    TH1F* nstclmvacon = new TH1F("nstclmvacon","Non-Stereo Cluster MVA output",200,-0.2,1.2);
    TH1F* nstclmvabkg = new TH1F("nstclmvabkg","Non-Stereo Cluster MVA output",200,-0.2,1.2);
    nstclmvacon->SetLineColor(kRed);
    nstclmvabkg->SetLineColor(kBlue);
    TCut stcl("ngdstereo/ngdhits>0.5");
    TCut nstcl("ngdstereo/ngdhits<0.5");
    ddiag->Project("stclmvacon","pmvaout",con+cluster+stcl);
    ddiag->Project("stclmvabkg","pmvaout",bkg+cluster+stcl);
    ddiag->Project("nstclmvacon","pmvaout",con+cluster+nstcl);
    ddiag->Project("nstclmvabkg","pmvaout",bkg+cluster+nstcl);
    TLegend* mcleg = new TLegend(0.2,0.7,0.6,0.9);
    mcleg->AddEntry(stclmvabkg,"#delta Background","L");
    mcleg->AddEntry(stclmvacon,"Conversion (X10)","L");
    Double_t factor(20.0);
    stclmvacon->Scale(factor);
    nstclmvacon->Scale(factor);
    TCanvas* can = new TCanvas("cmvacan","Cluster MVA output",800,400);
    can->Divide(2,1);
    can->cd(1);
    stclmvabkg->Draw();
    stclmvacon->Draw("same");
    mcleg->Draw();
    can->cd(2);
    nstclmvabkg->Draw();
    nstclmvacon->Draw("same");


  } else if(spage=="hitmva") {
    TH1F* sthtmvap = new TH1F("sthtmvap","Stereo Hit MVA output",200,-0.2,1.2);
    TH1F* sthtmvau = new TH1F("sthtmvau","Stereo Hit MVA output",200,-0.2,1.2);
    sthtmvau->SetLineColor(kRed);
    sthtmvap->SetLineColor(kBlue);
    TH1F* nsthtmvap = new TH1F("nsthtmvap","Non-Stereo Hit MVA output",200,-0.2,1.2);
    TH1F* nsthtmvau = new TH1F("nsthtmvau","Non-Stereo Hit MVA output",200,-0.2,1.2);
    nsthtmvap->SetLineColor(kBlue);
    nsthtmvau->SetLineColor(kRed);

    ddiag->Project("sthtmvap","_hgd",bkg+cluster+stht+primary);
    ddiag->Project("sthtmvau","_hgd",bkg+cluster+stht+norel);
    ddiag->Project("nsthtmvap","_hgd",bkg+cluster+nstht+primary);
    ddiag->Project("nsthtmvau","_hgd",bkg+cluster+nstht+norel);

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
//    TH1F* drhoconp = new TH1F("drhoconp","CE Hit #rho difference;#Delta #rho (mm)",100,0,150.0);
//    TH1F* dphiconp = new TH1F("dphiconp","CE Hit #phi difference;#Delta #phi",100,0,0.3);
//    TH1F* dtconp = new TH1F("dtconp","CE Hit time difference;#Delta t (nsec)",100,0,50);
//    TH1F* drhoconn = new TH1F("drhoconn","CE Hit #rho difference;#Delta #rho (mm)",100,0,150.0);
//    TH1F* dphiconn = new TH1F("dphiconn","CE Hit #phi difference;#Delta #phi",100,0,0.3);
//    TH1F* dtconn = new TH1F("dtconn","CE Hit time difference;#Delta t (nsec)",100,0,50);

    TH1F* drhobkgps = new TH1F("drhobkgps","Bkg Hit #rho difference;#Delta #rho (mm)",100,0,80.0);
    TH1F* dphibkgps = new TH1F("dphibkgps","Bkg Hit #phi difference;#Delta #phi",100,0,80.0);
    TH1F* drhobkgpn = new TH1F("drhobkgpn","Bkg Hit #rho difference;#Delta #rho (mm)",100,0,80.0);
    TH1F* dphibkgpn = new TH1F("dphibkgpn","Bkg Hit #phi difference;#Delta #phi",100,0,80.0);
    TH1F* cdbkgps = new TH1F("cdbkgps","Bkg Hit cluster distance;distance",100,0,8.0);
    TH1F* cdbkgpn = new TH1F("cdbkgpn","Bkg Hit cluster distance;distance",100,0,8.0);
    TH1F* dtbkgp = new TH1F("dtbkgp","Bkg Hit time difference;#Delta t (nsec)",100,0,50);
    TH1F* drhobkgn = new TH1F("drhobkgn","Unrelated Hit #rho difference;#Delta #rho (mm)",100,0,80.0);
    TH1F* dphibkgn = new TH1F("dphibkgn","Unrelated Hit #phi difference;#Delta #phi",100,0,80.0);
    TH1F* dtbkgn = new TH1F("dtbkgn","Unrelated Hit time difference;#Delta t (nsec)",100,0,50);
    TH1F* cdbkgn = new TH1F("cdbkgn","Unrelated Hit cluster distance;distance)",100,0,8.0);
    drhobkgps->SetLineColor(kBlue);
    dphibkgps->SetLineColor(kBlue);
    drhobkgpn->SetLineColor(kGreen);
    dphibkgpn->SetLineColor(kGreen);
    cdbkgps->SetLineColor(kBlue);
    cdbkgpn->SetLineColor(kGreen);
    dtbkgp->SetLineColor(kBlue);
    drhobkgn->SetLineColor(kRed);
    dphibkgn->SetLineColor(kRed);
    dtbkgn->SetLineColor(kRed);
    cdbkgn->SetLineColor(kRed);

    ddiag->Project("drhobkgps","_drho",bkg+cluster+primary+stht);
    ddiag->Project("dphibkgps","_dphi*rmean",bkg+cluster+primary+stht);
    ddiag->Project("cdbkgps","_cdist",bkg+cluster+primary+stht);
    ddiag->Project("drhobkgpn","_drho",bkg+cluster+primary+nstht);
    ddiag->Project("dphibkgpn","_dphi*rmean",bkg+cluster+primary+nstht);
    ddiag->Project("cdbkgpn","_cdist",bkg+cluster+primary+nstht);
    ddiag->Project("dtbkgp","_dt",bkg+cluster+primary);
    ddiag->Project("drhobkgn","_drho",bkg+cluster+norel);
    ddiag->Project("dphibkgn","_dphi*rmean",bkg+cluster+norel);
    ddiag->Project("dtbkgn","_dt",bkg+cluster+norel);
    ddiag->Project("cdbkgn","_cdist",bkg+cluster+norel);

    TCanvas* dhcan = new TCanvas("dhcan","Delta hits",1200,800);
    dhcan->Divide(2,2);
    dhcan->cd(1);
    drhobkgps->Draw();
    drhobkgpn->Draw("same");
    drhobkgn->Draw("same");


    TLegend* hleg = new TLegend(0.2,0.7,0.6,0.9);
    hleg->AddEntry(drhobkgps,"#delta Primary stereo hit","L");
    hleg->AddEntry(drhobkgpn,"#delta Primary non-stereo hit","L");
    hleg->AddEntry(drhobkgn,"#delta unrelated hit","L");
    hleg->Draw();

    dhcan->cd(2);
    dphibkgps->Draw();
    dphibkgpn->Draw("same");
    dphibkgn->Draw("same");
    dhcan->cd(3);
    dtbkgp->Draw();
    dtbkgn->Draw("same");
    dhcan->cd(4);
    cdbkgps->Draw();
    cdbkgpn->Draw("same");
    cdbkgn->Draw("same");
  }    
}
