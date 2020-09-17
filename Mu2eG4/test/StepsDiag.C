//
//
//
// Original author KLG 
//
// to Retrieve/prepare histograms based on histograms/ntuples from the
// files that were created by g4validate_01.fcl
//
// to be used with g4validate_01.C which handles the plot creation/overlay
//

#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include <TNtuple.h>
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

// Steps & (original) StrawHits diagonstic plot retrieval/preparation

bool getPushNextHist(TString hname, TDirectory* tdir, std::vector<TH1*>& plots) {
  TH1F*    _tmp = 0x0;
  _tmp = static_cast<TH1F*>(tdir->Get(hname));
  if (_tmp == 0x0 ) {cerr << "missing histogram" <<endl; return false;} 
  plots.push_back(_tmp);
  return true;
}

void StepsDiag(TFile* tfile,std::vector<TH1*>& plots) {

  TString sdirn("checkhits");
  TDirectory* sdir = static_cast<TDirectory*>(tfile->Get(sdirn));
  if(sdir == 0x0){
    std::cout << "StepsDiag: TDirectory " << sdirn << " not found" << std::endl;
    return;
  }

  TString shdirn("readStrawHits");
  TDirectory* shdir = static_cast<TDirectory*>(tfile->Get(shdirn));
  if(sdir == 0x0){
    std::cout << "StepsDiag: TDirectory " << sdirn << " not found" << std::endl;
    return;
  }

  TString tntuplen("ntup");
  TNtuple* hits = static_cast<TNtuple*>((shdir->Get(tntuplen)));
  if(hits == 0x0){
    std::cout << "StepsDiag: " << tntuplen << " TNtuple not found!" << std::endl;
    return;
  }

  // checkhits plots

  if (!getPushNextHist("hMultiplicity",   sdir, plots)) return;
  if (!getPushNextHist("hHitNeighbours",  sdir, plots)) return;
  if (!getPushNextHist("hEnergyDep",      sdir, plots)) return;
  if (!getPushNextHist("hStepLength",     sdir, plots)) return;
  if (!getPushNextHist("hRadius",         sdir, plots)) return;
  if (!getPushNextHist("hxHit",           sdir, plots)) return;
  if (!getPushNextHist("hyHit",           sdir, plots)) return;
  if (!getPushNextHist("hzHit",           sdir, plots)) return;
					  
  // readStrawHits plots
					  
  if (!getPushNextHist("hHitTime",        shdir, plots)) return; 
  if (!getPushNextHist("hHitDeltaTime",   shdir, plots)) return; 
  if (!getPushNextHist("hHitEnergy",      shdir, plots)) return; 
  if (!getPushNextHist("hNHits",          shdir, plots)) return; 
  if (!getPushNextHist("hNHitsPerWire",   shdir, plots)) return; 
  if (!getPushNextHist("hDriftTime",      shdir, plots)) return; 
  if (!getPushNextHist("hDriftDistance",  shdir, plots)) return; 
  if (!getPushNextHist("hDistanceToMid",  shdir, plots)) return; 
  if (!getPushNextHist("hNG4Steps",       shdir, plots)) return; 
  if (!getPushNextHist("hG4StepLength",   shdir, plots)) return; 
  if (!getPushNextHist("hG4StepRelTimes", shdir, plots)) return; 
  if (!getPushNextHist("hG4StepEdep",     shdir, plots)) return; 

  // ntuple based plot

  TH1F* hitR = new TH1F("hHitR","sqrt(hity*hity+hitx*hitx);mm",500,300.,800.);
  hits->Project("hHitR","sqrt(hity*hity+hitx*hitx)");
  plots.push_back(hitR);

}
