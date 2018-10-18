#define SHCombo_cxx
#include "SHCombo.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <array>

void SHCombo::MakeHists(){
  _dsa = new TH1F("dsa","StrawId difference",101,-0.5,100.5);
  _dst = new TH1F("dst","StrawId difference",101,-0.5,100.5);
  _dso = new TH1F("dso","StrawId difference",101,-0.5,100.5);
  _dta = new TH1F("dta","Time difference",100,-100.0,100.0);
  _dtt = new TH1F("dtt","Time difference",100,-100.0,100.0);
  _dto = new TH1F("dto","Time difference",100,-100.0,100.0);
  _dwda = new TH1F("dwda","Wire length difference",100,-400.0,400.0);
  _dwdt = new TH1F("dwdt","Wire length difference",100,-400.0,400.0);
  _dwdo = new TH1F("dwdo","Wire length difference",100,-400.0,400.0);
  _dedepa = new TH1F("dedepa","EDep difference",100,-10.0,10.0);
  _dedept = new TH1F("dedept","EDep difference",100,-10.0,10.0);
  _dedepo = new TH1F("dedepo","EDep difference",100,-10.0,10.0);
  _dwdpa = new TH1F("dwdpa","Wire length difference pull",100,-10.0,10.0);
  _dwdpt = new TH1F("dwdpt","Wire length difference pull",100,-10.0,10.0);
  _dwdpo = new TH1F("dwdpo","Wire length difference pull",100,-10.0,10.0);
  _np = new TH1F("np","N hits/panel",100,-0.5,99.5);
  _mcdwd = new TH1F("mcdwd","MC Wire length difference",100,-200.0,200.0);
  _wres = new TH1F("wres","Wire Resolution",100,-400.0,400.0);
  _awres = new TH1F("awres","Average Wire Resolution",100,-400.0,400.0);
  _wpull = new TH1F("wpull","Wire Pull",100,-10.0,10.0);
  _awpull = new TH1F("awpull","Average Wire Pull",100,-10.0,10.0);
  _nmatch = new TH1F("nmatch","N Matching",10,-0.5,9.5);
  _nmatcht = new TH1F("nmatcht","N Matching",10,-0.5,9.5);
  _nmatchst = new TH1F("nmatchstt","N Matching",10,-0.5,9.5);

  _dsa->SetStats(0);
  _dta->SetStats(0);
  _dwda->SetStats(0);
  _dwdpa->SetStats(0);
  _dedepa->SetStats(0);

  _dsa->SetLineColor(kRed);
  _dta->SetLineColor(kRed);
  _dwda->SetLineColor(kRed);
  _dwdpa->SetLineColor(kRed);
  _dedepa->SetLineColor(kRed);
  _dso->SetLineColor(kGreen);
  _dto->SetLineColor(kGreen);
  _dwdo->SetLineColor(kGreen);
  _dwdpo->SetLineColor(kGreen);
  _dedepo->SetLineColor(kGreen);
  _dst->SetLineColor(kBlue);
  _dtt->SetLineColor(kBlue);
  _dwdt->SetLineColor(kBlue);
  _dwdpt->SetLineColor(kBlue);
  _dedept->SetLineColor(kBlue);

  _nmatch->SetLineColor(kGreen);
  _nmatcht->SetLineColor(kBlue);
  _nmatchst->SetLineColor(kOrange);

  _mcdwd->SetLineColor(kRed);
  _wres->SetLineColor(kOrange);
  _wpull->SetLineColor(kOrange);
  _awres->SetLineColor(kCyan);
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
      _plane.clear(); _panel.clear(); _straw.clear(); _mcid.clear(); _wd.clear(); _wderr.clear(); _time.clear(); _mcwd.clear(); _edep.clear();
    }
    // store data for this straw hit
    _plane.push_back(plane);
    _panel.push_back(panel);
    _straw.push_back(straw);
    _mcid.push_back(mcid);
    _mcwd.push_back(mcshlen);
    _edep.push_back(1000.0*edep);
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
	int ds = abs( (int)_straw[ihit]-(int)_straw[jhit]);
	float dt = _time[ihit]-_time[jhit];
	float dw = _wd[ihit]-_wd[jhit];
	float dwe1 = _wderr[ihit]*_wderr[ihit];
	float dwe2 = _wderr[jhit]*_wderr[jhit];
	float dwerr = sqrt(dwe1+dwe2);
	float awerr = 1.0/sqrt(1.0/dwe1+1.0/dwe2);
	float dwp = dw/dwerr;
	float wa = (_wd[ihit]/dwe1 + _wd[jhit]/dwe2)/(1.0/dwe1+1.0/dwe2);
	float mcwa = 0.5*(_mcwd[ihit]+_mcwd[jhit]);
	float dedep = _edep[ihit]-_edep[jhit];
	_dsa->Fill(ds);
	_dta->Fill(dt);
	_dwdpa->Fill(dwp);
	_dwda->Fill(dw);
	_dedepa->Fill(dedep);
	if(mcmatch){
	  _mcdwd->Fill(_mcwd[ihit]-_mcwd[jhit]);
	  nmatcht++;
	  _dst->Fill(ds);
	  _dtt->Fill(dt);
	  _dwdpt->Fill(dwp);
	  _dwdt->Fill(dw);
	  _dedept->Fill(dedep);
	}
	bool goodds = ds !=0 && ds <= _maxds;
	bool gooddt = fabs(dt) < _maxdt;
	bool gooddwp = abs(dwp) < _maxdwp;
	if(goodds&&gooddt) _dwdpo->Fill(dwp);
	if(goodds&&gooddwp) _dto->Fill(dt);
	if(gooddt&&gooddwp) _dso->Fill(ds);
	if(gooddt&&gooddwp&&goodds){
	  _dwdo->Fill(dw);
	  _dedepo->Fill(dedep);
	  ++nmatch;
	  _awpull->Fill((wa-mcwa)/awerr);
	  _awres->Fill(wa-mcwa);
	  if(mcmatch)++nmatchst;
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
  TCanvas* selcan = new TCanvas("selcan","selcan",1000,700);
  _cans.push_back(selcan);

  TLegend* sleg = new TLegend(0.4,0.6,0.9,0.9);
  sleg->AddEntry(_dsa,"All Panel Hit Pairs","L");
  sleg->AddEntry(_dso,"Selected Hit Pairs","L");
  sleg->AddEntry(_dst,"True Panel Hit Pairs","L");
  selcan->Divide(3,2);
  selcan->cd(1);
  _dsa->Draw();
  _dst->Draw("same");
  _dso->Draw("same");
  sleg->Draw();
  selcan->cd(2);
  _dta->Draw();
  _dtt->Draw("same");
  _dto->Draw("same");
  selcan->cd(3);
  _dwdpa->Draw();
  _dwdpt->Draw("same");
  _dwdpo->Draw("same");
  selcan->cd(4);
  _dwda->Draw();
  _dwdt->Draw("same");
  _dwdo->Draw("same");
  selcan->cd(5);
  _dedepa->Draw();
  _dedept->Draw("same");
  _dedepo->Draw("same");


  TCanvas* rcan = new TCanvas("rcan","rcan",700,700);
  _cans.push_back(rcan);
  rcan->Divide(2,2);
  rcan->cd(1);
  _wres->Draw();
  _awres->Draw("sames");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)_wres->FindObject("stats");
  Coord_t dx = st->GetX2NDC()-st->GetX1NDC();
  st->SetX2NDC(st->GetX1NDC());
  st->SetX1NDC(st->GetX1NDC()-dx);
  st->Draw();
  TLegend* rleg = new TLegend(0.1,0.6,0.4,0.9);
  rleg->AddEntry(_wres,"Single Hit","L");
  rleg->AddEntry(_awres,"Pair Average","L");
  rleg->Draw();
  rcan->cd(2);
  _wpull->Draw();
  _awpull->Draw("sames");
  gPad->Update();
  st = (TPaveStats*)_wpull->FindObject("stats");
  dx = st->GetX2NDC()-st->GetX1NDC();
  st->SetX2NDC(st->GetX1NDC());
  st->SetX1NDC(st->GetX1NDC()-dx);
  st->Draw();
  rcan->cd(3);
  _mcdwd->Draw();

  TCanvas* ncan = new TCanvas("ncan","ncan",700,700);
  _cans.push_back(ncan);
  ncan->Divide(2,1);
  ncan->cd(1);
  _nmatcht->Draw();
  _nmatch->Draw("same");
  _nmatchst->Draw("same");
  TLegend* nleg = new TLegend(0.3,0.6,0.8,0.9);
  nleg->AddEntry(_nmatcht,"MC True","L"); 
  nleg->AddEntry(_nmatch,"Reco","L"); 
  nleg->AddEntry(_nmatchst,"True Reco","L");
  nleg->Draw();
  ncan->cd(2);
  _np->Draw();
}

