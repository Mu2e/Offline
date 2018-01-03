#define SHCombo_cxx
#include "SHCombo.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <array>

void SHCombo::MakeHists(){
  _ds = new TH1F("ds","StrawId difference",101,-0.5,100.5);
  _dst = new TH1F("dst","StrawId difference",101,-0.5,100.5);
  _dt = new TH1F("dt","Time difference",100,-100.0,100.0);
  _dtt = new TH1F("dtt","Time difference",100,-100.0,100.0);
  _dwd = new TH1F("dwd","Wire length difference",100,-400.0,400.0);
  _dwdt = new TH1F("dwdt","Wire length difference",100,-400.0,400.0);
  _dwdp = new TH1F("dwdp","Wire length difference pull",100,-10.0,10.0);
  _dwdpt = new TH1F("dwdpt","Wire length difference pull",100,-10.0,10.0);
  _np = new TH1F("np","N hits/panel",100,-0.5,99.5);
  _mcdwd = new TH1F("mcdwd","MC Wire length difference",100,-200.0,200.0);
  _wres = new TH1F("wres","Wire Resolution",100,-400.0,400.0);
  _awres = new TH1F("awres","Average Wire Resolution",100,-400.0,400.0);
  _wpull = new TH1F("wpull","Wire Pull",100,-10.0,10.0);
  _awpull = new TH1F("awpull","Average Wire Pull",100,-10.0,10.0);
  _nmatch = new TH1F("nmatch","N Matching",10,-0.5,9.5);
  _nmatcht = new TH1F("nmatcht","N Matching",10,-0.5,9.5);
  _nmatchst = new TH1F("nmatchstt","N Matching",10,-0.5,9.5);
  _dst->SetLineColor(kBlue);
  _dtt->SetLineColor(kBlue);
  _dwdt->SetLineColor(kBlue);
  _dwdpt->SetLineColor(kBlue);
  _mcdwd->SetLineColor(kRed);
  _nmatch->SetLineColor(kGreen);
  _nmatcht->SetLineColor(kBlue);
  _nmatchst->SetLineColor(kRed);
  _wres->SetLineColor(kRed);
  _awres->SetLineColor(kCyan);
  _wpull->SetLineColor(kRed);
  _awpull->SetLineColor(kCyan);

}

void SHCombo::SaveCans(const char* suffix) {
  std::string ssuf(suffix);
  for(auto can : _cans){
    std::string fname(can->GetName());
    fname += ssuf;
    can->SaveAs(fname.c_str());
  }
}

void SHCombo::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  // get ID of 1st event
  Long64_t ientry = LoadTree(0);
  fChain->GetEntry(0);
  int run = runid;
  int subrun = subrunid;
  int event = eventid;
  for(Long64_t jentry=0;jentry < nentries; ++jentry){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if(event != eventid || subrun != subrunid || run != runid){
      // end of this event: process the data
      processEvent();
      // advance to the next event
      event = eventid;
      subrun = subrunid;
      run = runid;
      _plane.clear(); _panel.clear(); _straw.clear(); _mcid.clear(); _wd.clear(); _wderr.clear(); _time.clear(); _mcwd.clear();
    }
    // store data for this straw hit
    _plane.push_back(plane);
    _panel.push_back(panel);
    _straw.push_back(straw);
    _mcid.push_back(mcid);
    _mcwd.push_back(mcshlen);
    _wd.push_back(shlen);
    _wderr.push_back(wres);
    _time.push_back(std::min(time_tcal,time_thv));
  }
  // process last event
  processEvent();
  // draw
  Draw();
}

void SHCombo::processEvent() {
    // loop over hit pairs
  std::array<unsigned,240> np{0};
  for(size_t ihit = 0;ihit < _plane.size(); ++ihit){
    _wres->Fill(_wd[ihit]-_mcwd[ihit]);
    _wpull->Fill((_wd[ihit]-_mcwd[ihit])/_wderr[ihit]);
    unsigned pid = 6*_plane[ihit]+_panel[ihit];	
    ++np[pid];
    unsigned nmatch(0), nmatchst(0), nmatcht(0);
    for(size_t jhit=ihit+1;jhit < _plane.size(); ++jhit){
      if(_plane[ihit] == _plane[jhit] &&
	  _panel[ihit] == _panel[jhit]){
	bool mcmatch = _mcid[ihit] == _mcid[jhit];
	if(mcmatch){
	  _mcdwd->Fill(_mcwd[ihit]-_mcwd[jhit]);
	  nmatcht++;
	}
	int ds = abs( (int)_straw[ihit]-(int)_straw[jhit]);
	_ds->Fill(ds);
	if(mcmatch) _dst->Fill(ds);
	if(ds !=0 && ds <= _maxds){
	  float dt = _time[ihit]-_time[jhit];
	  _dt->Fill(dt);
	  if(mcmatch) _dtt->Fill(dt);
	  if(fabs(dt) < _maxdt){
	    float dw = _wd[ihit]-_wd[jhit];
	    float dwe1 = _wderr[ihit]*_wderr[ihit];
	    float dwe2 = _wderr[jhit]*_wderr[jhit];
	    float dwerr = sqrt(dwe1+dwe2);
	    float awerr = 1.0/sqrt(1.0/dwe1+1.0/dwe2);
	    float dwp = dw/dwerr;
	    float wa = (_wd[ihit]/dwe1 + _wd[jhit]/dwe2)/(1.0/dwe1+1.0/dwe2);
	    float mcwa = 0.5*(_mcwd[ihit]+_mcwd[jhit]);
	    _dwdp->Fill(dwp);
	    if(mcmatch)_dwdpt->Fill(dwp);
	    if(fabs(dwp) < _maxdwp) {
	      _dwd->Fill(dw);
	      ++nmatch;
	      _awres->Fill(wa-mcwa);
	      _awpull->Fill((wa-mcwa)/awerr);
	      if(mcmatch){
		_dwdt->Fill(dw);
		++nmatchst;
	      }
	    }
	  }
	}
      }
    }
    _nmatch->Fill(nmatch);
    _nmatcht->Fill(nmatcht);
    _nmatchst->Fill(nmatchst);
  }
  for(auto nh : np)
    _np->Fill(nh);
}  
void SHCombo::Draw() {
  // draw results
  TCanvas* selcan = new TCanvas("selcan","selcan",700,700);
  _cans.push_back(selcan);
  selcan->Divide(2,2);
  selcan->cd(1);
  _ds->Draw();
  _dst->Draw("same");
  selcan->cd(2);
  _dt->Draw();
  _dtt->Draw("same");
  selcan->cd(3);
  _dwdp->Draw();
  _dwdpt->Draw("same");
  selcan->cd(4);
  _dwd->Draw();
  _dwdt->Draw("same");

  TCanvas* rcan = new TCanvas("rcan","rcan",700,700);
  _cans.push_back(rcan);
  rcan->Divide(2,2);
  rcan->cd(1);
  _wres->Draw();
  _awres->Draw("same");
  rcan->cd(2);
  _wpull->Draw();
  _awpull->Draw("same");
  rcan->cd(3);
  _mcdwd->Draw();

  TCanvas* ncan = new TCanvas("ncan","ncan",700,700);
  _cans.push_back(ncan);
  ncan->Divide(2,1);
  ncan->cd(1);
  _nmatcht->Draw();
  _nmatch->Draw("same");
  _nmatchst->Draw("same");
  ncan->cd(2);
  _np->Draw();
}

