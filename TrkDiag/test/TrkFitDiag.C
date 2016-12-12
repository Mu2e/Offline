//
// $Id: TrkFitDiag.C,v 1.5 2013/02/14 22:22:07 genser Exp $
// $Author: genser $
// $Date: 2013/02/14 22:22:07 $
//
//
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TFile.h"
#include <vector>
#include <iostream>

// Track Fit diagonstics

void TrkFitDiag(TFile* tfile,std::vector<TH1*>& plots) {

  TString tdirn("RKFDownstreameMinus");
  TDirectory* tdir = static_cast<TDirectory*>(tfile->Get(tdirn));
  if(tdir == 0){
    std::cout << "TrkFitDiag: TDirectory " << tdirn << " not found" << std::endl;
    return;
  }

  TString ttreen("trkdiag");
  TTree* trks = static_cast<TTree*>((tdir->Get(ttreen)));
  if(trks == 0){
    std::cout << "TrkFitDiag: " << ttreen << " TTree not found!" << std::endl;
    return;
  }

// setup cuts

  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(720);
  double momlow(103.4);
  double momhigh(104.8);
  unsigned minnhits(0);
  unsigned minnactive(25);
  unsigned maxndif(10);
  double maxt0err(0.9);
  double maxmomerr(0.3);
  double minfitcon(1e-3);
  TCut reco,goodfit,cosmic,rmom,rpitch,livegate;
  TCut goodmc, tpitch, tmom, nmch;
  TCut nacut, t0errcut, momerrcut,fitconcut;
  char cutstring[100];
  snprintf(cutstring,100,"nactive>=%i&&nhits-nactive<=%i",minnactive,maxndif);
  nacut = TCut(cutstring);
  snprintf(cutstring,100,"t0err<%f",maxt0err);
  t0errcut = TCut(cutstring);
  snprintf(cutstring,100,"fit.momerr<%f",maxmomerr);
  momerrcut = TCut(cutstring);
  snprintf(cutstring,100,"fit.con>%f",minfitcon);
  fitconcut = TCut(cutstring);

  snprintf(cutstring,100,"td>%4.3f&&td<%4.3f",tdlow,tdhigh);
  rpitch = TCut(cutstring);
  snprintf(cutstring,100,"t0>%f",t0min);
  livegate = TCut(cutstring);
  snprintf(cutstring,100,"mcent.td>%4.3f&&mcent.td<%4.3f",tdlow-0.02,tdhigh+0.02);
  tpitch = TCut(cutstring);
  tmom = TCut("mcent.mom>100");
  snprintf(cutstring,100,"mc.ndigigood>=%i",minnhits);
  nmch = TCut(cutstring);
  reco = TCut("fit.status>0");
  cosmic = TCut("abs(d0)<105 && d0+2/om>450 && d0+2/om<680");
  snprintf(cutstring,100,"fit.mom>%f&&fit.mom<%f",momlow,momhigh);
  rmom = TCut(cutstring);

  goodmc = tpitch+tmom+nmch;
  goodfit = reco+nacut+t0errcut+momerrcut+fitconcut;

  //cout << "cuts done " << endl;

  TH1F* fitstatus = new TH1F("fitstatus","Kalman fit status",11,-5.5,5.5);
  TH1F* chisq = new TH1F("chisq","Chisq/NDof",100,0,10);
  TH1F* fitcon = new TH1F("fitcon","fit consistency",101,-0.01,1.01);
  TH1F* fitconl = new TH1F("fitconl","log_{10} fit consistency",101,-8,0);
  TH1F* t0 = new TH1F("t0","Track Fit t_{0};t_{0} (ns)",100,0,1800);
  TH1F* t0err = new TH1F("t0err","Fit t_{0} error; t_{0} error (ns)",100,0,2.0);
  TH1F* nh = new TH1F("nh","Fit N hits;N hits",101,-0.5,100.5);
  TH1F* na = new TH1F("na","Fit N active hits;N hits",101,-0.5,100.5);
  TH1F* nmc = new TH1F("nmc","N CE hits;N hits",101,-0.5,100.5);

  //  cout << "first plots booked " << endl;

  plots.push_back(fitstatus);
  plots.push_back(chisq);
  plots.push_back(fitcon);
  plots.push_back(fitconl);
  plots.push_back(t0);
  plots.push_back(t0err);
  plots.push_back(nh);
  plots.push_back(na);
  plots.push_back(nmc);

  //  cout << "first plots pushed " << endl;

  trks->Project("fitstatus","fit.status");
  trks->Project("chisq","chisq/ndof");
  trks->Project("fitcon","fit.con",goodmc);
  trks->Project("fitconl","log10(fit.con)",goodmc);
  trks->Project("t0","t0",goodmc+goodfit);
  trks->Project("t0err","t0err",goodmc);
  trks->Project("nh","nhits",goodmc);
  trks->Project("na","nactive",goodmc);
  trks->Project("nmc","mc.ndigigood");

  //  cout << "first plots done" << endl;

  TH1F* d0 = new TH1F("d0","Track fit d_{0};d_{0} (mm)",100,-150,150);
  TH1F* p0 = new TH1F("p0","Track Fit #phi_{0};#phi_{0}",100,-3.15,3.15);
  TH1F* om = new TH1F("om","Track Fit #omega;#omega (mm^{-1})",100,-0.006,0.006);
  TH1F* z0 = new TH1F("z0","Track Fit z_{0};z_{0} (mm)",100,-1000,1000);
  TH1F* td = new TH1F("td","Track Fit tan(#lambda);tan(#lambda)",100,0.5,1.5);
  TH1F* rmax = new TH1F("rmax","Track fit rmax;d_{0}+2/#omega (mm)",100,300,900);
  plots.push_back(d0);
  plots.push_back(p0);
  plots.push_back(om);
  plots.push_back(z0);
  plots.push_back(td);
  plots.push_back(rmax);
  trks->Project("d0","d0",goodmc+reco);
  trks->Project("p0","p0",goodmc+reco);
  trks->Project("om","om",goodmc+reco);
  trks->Project("z0","z0",goodmc+reco);
  trks->Project("td","td",goodmc+reco);
  trks->Project("rmax","d0+2.0/om",goodmc+reco);

  TH1F* d0pull = new TH1F("d0pull","d0 pull",100,-10,10);
  TH1F* p0pull = new TH1F("p0pull","#phi0 pull",100,-10,10);
  TH1F* ompull = new TH1F("ompull","#omega pull",100,-10,10);
  TH1F* z0pull = new TH1F("z0pull","z0 pull",100,-10,10);
  TH1F* tdpull = new TH1F("tdpull","tan(#lambda) pull",100,-10,10);
  plots.push_back(d0pull);
  plots.push_back(p0pull);
  plots.push_back(ompull);
  plots.push_back(z0pull);
  plots.push_back(tdpull);
  trks->Project("d0pull","(d0-mcent.d0)/d0err",goodmc+goodfit);
  trks->Project("p0pull","(p0-mcent.p0)/p0err",goodmc+goodfit);
  trks->Project("ompull","(om-mcent.om)/omerr",goodmc+goodfit);
  trks->Project("z0pull","(z0-mcent.z0)/z0err",goodmc+goodfit);
  trks->Project("tdpull","(td-mcent.td)/tderr",goodmc+goodfit);
//  d0pull->Fit("gaus","Q0");
//  p0pull->Fit("gaus","Q0");
//  ompull->Fit("gaus","Q0");
//  z0pull->Fit("gaus","Q0");
//  tdpull->Fit("gaus","Q0");

  TH1F* fitmom = new TH1F("fitmom","Track fit momentum at tracker entrance;fit momentum (MeV)",100,90,107);
  TH1F* mcmom = new TH1F("mcmom","True CE momentum at tracker entrance;CE momentum (MeV)",100,90,107);
  TH1F* momerr = new TH1F("momerr","Fit momentum error;momentum error (MeV)",100,0,0.5);
  TH1F* mres = new TH1F("mres","momentum resolution at tracker entrance;MeV",100,-2,2);
  TH1F* mpull = new TH1F("mpull","momentum pull at tracker entrance",100,-10,10);
  plots.push_back(fitmom);
  plots.push_back(mcmom);
  plots.push_back(momerr);
  plots.push_back(mres);
  plots.push_back(mpull);
  trks->Project("fitmom","fit.mom",goodmc+goodfit);
  trks->Project("mcmom","mcent.mom",tpitch+nmch);
  trks->Project("momerr","fit.momerr",goodmc+goodfit);
  trks->Project("mres","fit.mom-mcent.mom",goodmc+goodfit);
  trks->Project("mpull","(fit.mom-mcent.mom)/fit.momerr",goodmc+goodfit);
//  mpull->Fit("gaus","Q0");

  unsigned nbins(10);
  double bmax = nbins-0.5;
  TH1F* acc = new TH1F("acc","Cummulative Acceptance;;cummulative acceptance",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","Relative Acceptance;;relative acceptance",nbins,-0.5,bmax);
  plots.push_back(acc);
  plots.push_back(racc);
  unsigned ibin(1);
  acc->GetXaxis()->SetBinLabel(ibin++,"All");
  acc->GetXaxis()->SetBinLabel(ibin++,">=20 SH");
  acc->GetXaxis()->SetBinLabel(ibin++,"p>100 MeV/c");
  acc->GetXaxis()->SetBinLabel(ibin++,"pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  acc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  acc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");

  ibin = 1;
  racc->GetXaxis()->SetBinLabel(ibin++,"All");
  racc->GetXaxis()->SetBinLabel(ibin++,">=20 SH");
  racc->GetXaxis()->SetBinLabel(ibin++,"p>100 MeV/c");
  racc->GetXaxis()->SetBinLabel(ibin++,"pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  racc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  racc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");
  
  ibin = 0;
  const char* binnames[11] ={"0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0","8.0","9.0","10.0"};
  trks->Project("acc",binnames[ibin++]);
  trks->Project("+acc",binnames[ibin++],nmch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+fitconcut);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+fitconcut+livegate);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+fitconcut+livegate+rpitch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+fitconcut+livegate+rpitch+cosmic);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+fitconcut+livegate+rpitch+cosmic+rmom);

  double all = acc->GetBinContent(1);
  double prev = all;
  for(ibin=1;ibin<=nbins;ibin++){
    racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    prev = acc->GetBinContent(ibin);
  }
  racc->SetMaximum(1.1);
  acc->Scale(1.0/all);
  acc->SetMaximum(1.1);
  acc->SetStats(0);
  racc->SetStats(0);
  acc->GetXaxis()->SetLabelSize(0.06);
  racc->GetXaxis()->SetLabelSize(0.06);
  acc->SetMarkerSize(2.0);
  racc->SetMarkerSize(2.0);
  acc->GetYaxis()->SetTitleSize(0.05);
  racc->GetYaxis()->SetTitleSize(0.05);

}
