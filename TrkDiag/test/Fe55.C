#define Fe55_cxx
#include "Fe55.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void Fe55::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Int_t oldevtid(-1);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
       Long64_t nb = fChain->GetEntry(jentry);
     // look for entries in the same event
      if(eventid == oldevtid){
//	cout << "found match eventid" << eventid << endl;
	std::array<HitComp,2> hc;
	fillHC(hc[0]);
//	cout << "before straw = " << straw << endl;
	Long64_t kentry = LoadTree(jentry-1);
	Long64_t knb = fChain->GetEntry(jentry-1);
//	cout << "after straw = " << straw << endl;
	fillHC(hc[1]);
	_eid = eventid;
	plotDiff(hc);
      }
      oldevtid = eventid;
   }
}

void Fe55::fillHC(HitComp& hc) {
  hc._straw = straw;
  hc._plane = plane;
  hc._panel = panel;
  hc._mcsptime = mcsptime;
  hc._mcgtime = mcgt;
  hc._mcpoe = mcpoe;
  hc._mcedep = mcedep;
  hc._mcshpos.dx = mcshpos_dx;
  hc._mcshpos.dy = mcshpos_dy;
  hc._mcshpos.dz = mcshpos_dz;
}

void Fe55::plotDiff(std::array<HitComp,2>& hc) {
  size_t ifirst(0), ilast(1);
  if(hc[0]._mcsptime > hc[1]._mcsptime){
    ifirst = 1;
    ilast = 0;
  }
  _dstraw->Fill(hc[ilast]._straw-hc[ifirst]._straw);
  _dplane->Fill(hc[ilast]._plane-hc[ifirst]._plane);
  _dpanel->Fill(hc[ilast]._panel-hc[ifirst]._panel);
  _dmcsptime->Fill(hc[ilast]._mcsptime-hc[ifirst]._mcsptime);
  _dmcgtime->Fill(hc[ilast]._mcgtime-hc[ifirst]._mcgtime);
  double dx = hc[ilast]._mcshpos.dx-hc[ifirst]._mcshpos.dx;
  double dy = hc[ilast]._mcshpos.dy-hc[ifirst]._mcshpos.dy;
  double dz = hc[ilast]._mcshpos.dz-hc[ifirst]._mcshpos.dz;
  _ddx->Fill(dx);
  _ddy->Fill(dy);
  _ddz->Fill(dz);
  _dist->Fill(sqrt(dx*dx+dy*dy+dz*dz));
  _earlyoe->Fill(1000.0*hc[ifirst]._mcpoe);
  _lateoe->Fill(1000.0*hc[ilast]._mcpoe);
  _earlyde->Fill(1000.0*hc[ifirst]._mcedep);
  _latede->Fill(1000.0*hc[ilast]._mcedep);
  _hc1 = hc[ifirst];
  _hc2 = hc[ilast];
  _hcomp->Fill();
}

void Fe55::bookPlots() {
  _dstraw = new TH1F("dstraw","#Delta straw number",201,-100,100);
  _dpanel = new TH1F("dpanel","#Delta panel number",11,-5,5);
  _dplane = new TH1F("dplane","#Delta plane number",41,-20,20);
  _dmcsptime = new TH1F("dmcsptime","#Delta mcsptime",200,-0.001,10.0);
  _dmcgtime = new TH1F("dmcgtime","#Delta mcgtime",200,-0.001,10.0);
  _ddx = new TH1F("ddx","#Delta dx",200,-0.001,100.0);
  _ddy = new TH1F("ddy","#Delta dy",200,-0.001,100.0);
  _ddz = new TH1F("ddz","#Delta dz",200,-0.001,100.0);
  _dist = new TH1F("dist","Distance",200,-0.001,100.0);
  _earlyoe = new TH1F("earlyoe","Early Origin Energy;KeV",100,0.0,10.0);
  _lateoe = new TH1F("lateoe","Late Origin Energy;KeV",100,0.0,10.0);
  _earlyde = new TH1F("earlyde","Early Deposited Energy;KeV",100,0.0,10.0);
  _latede = new TH1F("latede","Late Deposited Energy;KeV",100,0.0,10.0);

  _hcomp = new TTree("hcomp","HitComp");
  _hcomp->Branch("eid",&_eid,"eid/I");
  _hcomp->Branch("hc1.",&_hc1);
  _hcomp->Branch("hc2.",&_hc2);
}
