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
#include <vector>

// Track Hit diagonstics

void TrkHitDiag(TDirectory& tdir,std::vector<TObject*> plots) {
  TTree* trks = tdir.Get("trkdiag");

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
  trks->Project("t0res","hitt0-mchitt0","_active");
  trks->Project("t0pull","(hitt0-mchitt0)/hitt0err","_active");
  t0res->Fit("gaus","Q0");
  t0pull->Fit("gaus","Q0");

  TH1F* tdivres = new TH1F("tdivres","#Deltat V resolution;mm",100,-200,200);
  TH1F* tdivpull = new TH1F("tdivpull","#Deltat V pull",100,-10,10);
  plots.push_back(tdivres);
  plots.push_back(tdivpull);
  trks->Project("tdivres","dmid-mcdmid","_active");
  trks->Project("tdivpull","(dmid-mcdmid)/dmiderr","_active");
  tdivres->Fit("gaus","Q0");
  tdivpull->Fit("gaus","Q0");

}
