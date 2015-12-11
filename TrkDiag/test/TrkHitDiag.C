//
// $Id: TrkHitDiag.C,v 1.3 2013/02/14 22:19:12 genser Exp $
// $Author: genser $
// $Date: 2013/02/14 22:19:12 $
//
//
#include "TCanvas.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include <vector>

// Track Hit diagonstics

void TrkHitDiag(TFile* tfile,std::vector<TH1*>& plots) {

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

  TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
  TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-10,10);
  TH1F* rpull = new TH1F("rpull","residual pull",100,-10,10);
  plots.push_back(dres);
  plots.push_back(dpull);
  plots.push_back(rpull);
  trks->Project("dres","_rdrift-_mcdist","_active");
  trks->Project("dpull","(_rdrift-_mcdist)/_rdrifterr","_active");
  trks->Project("rpull","_resid/_residerr","_active");
  dres->Fit("gaus","Q0");
  dpull->Fit("gaus","Q0");
  rpull->Fit("gaus","Q0");

  TH1F* t0res = new TH1F("t0res","hit t0 resolution;nsec",100,-10,10);
  TH1F* t0pull = new TH1F("t0pull","hit t0 pull",100,-10,10);
  plots.push_back(t0res);
  plots.push_back(t0pull);
  trks->Project("t0res","_ht-_mcht","_active");
  trks->Project("t0pull","(_ht-_mcht)/_t0err","_active");
  t0res->Fit("gaus","Q0");
  t0pull->Fit("gaus","Q0");

  TH1F* tdivres = new TH1F("tdivres","#Deltat V resolution;mm",100,-200,200);
  TH1F* tdivpull = new TH1F("tdivpull","#Deltat V pull",100,-10,10);
  plots.push_back(tdivres);
  plots.push_back(tdivpull);
  trks->Project("tdivres","_tddist-_mcdist","_active");
  trks->Project("tdivpull","(_tddist-_mcdist)/_tdderr","_active");
  tdivres->Fit("gaus","Q0");
  tdivpull->Fit("gaus","Q0");

}
