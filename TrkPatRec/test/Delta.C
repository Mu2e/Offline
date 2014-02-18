#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
void Delta(TTree* ddiag, const char* page="rho",const char* addcut="") {
  TString spage(page);
  TCut add(addcut);
  TCut con("pgen==2&&nprimary/nchits>0.8");
  con += add;
  TCut bkg("pproc<20&&nprimary/nchits>0.8");
  bkg += add;
  TCut cluster("nchits>=50&&nprimary/nchits>0.99");

  if(spage == "rho"){
    TH1F* mrhocon = new TH1F("mrhocon","Cluster #rho;#rho (mm)",100,330,800);
    TH1F* mrhobkg = new TH1F("mrhobkg","Cluster #rho;#rho (mm)",100,330,800);
    mrhocon->SetLineColor(kRed);
    mrhobkg->SetLineColor(kBlue);

    mrhocon->SetStats(0);
    mrhobkg->SetStats(0);

    ddiag->Project("mrhocon","rmed",con);
    ddiag->Project("mrhobkg","rmed",bkg);

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

    ddiag->Project("srhocon","srho",con);
    ddiag->Project("srhobkg","srho",bkg);
    ddiag->Project("sphicon","sphi",con);
    ddiag->Project("sphibkg","sphi",bkg);

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

    ddiag->Project("zmincon","zmin",con);
    ddiag->Project("zminbkg","zmin",bkg);
    ddiag->Project("zmaxcon","zmax",con);
    ddiag->Project("zmaxbkg","zmax",bkg);
    ddiag->Project("zgapcon","zgap",con);
    ddiag->Project("zgapbkg","zgap",bkg);

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

    ddiag->Project("smincon","smin",con);
    ddiag->Project("sminbkg","smin",bkg);
    ddiag->Project("smaxcon","smax",con);
    ddiag->Project("smaxbkg","smax",bkg);
    ddiag->Project("nsmisscon","nsmiss",con);
    ddiag->Project("nsmissbkg","nsmiss",bkg);
    ddiag->Project("nscon","ns",con);
    ddiag->Project("nsbkg","ns",bkg);

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

    ddiag->Project("ngdhitscon","ngdhits",con);
    ddiag->Project("ngdhitsbkg","ngdhits",bkg);
    ddiag->Project("nchitscon","nchits",con);
    ddiag->Project("nchitsbkg","nchits",bkg);
    ddiag->Project("cwratiocon","ngdhits/nchits",con);
    ddiag->Project("cwratiobkg","ngdhits/nchits",bkg);

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
  }
}
