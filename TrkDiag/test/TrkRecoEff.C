#define TrkRecoEff_cxx
#include "TrkRecoEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

TrkRecoEff::TrkRecoEff(TTree *tree, double norm) : fChain(0) ,_tffval(32,0), _norm(norm),
_eff(0), _rej(0), _cuts(0)
{
// build branches 
   Init(tree);
   // setup flag value vector: this should come from TrkFitFlag FIXME!!
   _tffval[hitsOK] = 1<<hitsOK;
   _tffval[circleOK] = 1<<circleOK;
   _tffval[phizOK] = 1<<phizOK;
   _tffval[helixOK] = 1<<helixOK;
   _tffval[seedOK] = 1<<seedOK;
   _tffval[kalmanOK] = 1<<kalmanOK;
   createHistos();
}


void TrkRecoEff::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::vector<bool> effcuts, rejcuts;
      effcuts.push_back(true); rejcuts.push_back(effcuts.back());// 0 bin counts number of events
      effcuts.push_back(ndigi >= 15); // events with at least 15 digis from the primary particle
      // take time modulo MB period
      float mctime = fmod(mcmidt0,1695.0);
      effcuts.push_back(mctime > 500.0 ); // primary track passed the tracker during 
      float  mccost = mcmidpz/mcmidmom;
      static float highcost = 1.0/sqrt(2.0)+0.05;
      effcuts.push_back(mccost > 0.45 && mccost < highcost ); // track direction (with buffer)
      effcuts.push_back(tcn > 10 ); rejcuts.push_back(effcuts.back());// reconstructed time cluster associated with primary track
      effcuts.push_back((hsf__value&_tffval[circleOK])>0);
      effcuts.push_back((hsf__value&_tffval[helixOK])>0);rejcuts.push_back(effcuts.back());
      effcuts.push_back((ksf__value&_tffval[seedOK])>0);rejcuts.push_back(effcuts.back());
      effcuts.push_back((kff__value&_tffval[kalmanOK])>0);rejcuts.push_back(effcuts.back());

      // effcuts are cumulative
      for(unsigned icut=0;icut< effcuts.size(); ++icut){
	Float_t fcut = (Float_t)icut;
	if(effcuts[icut])
	  _eff->Fill(fcut);
	else
	  break;
      }
      // as are rejection cuts
      for(unsigned icut=0;icut< rejcuts.size(); ++icut){
	Float_t fcut = (Float_t)icut;
	if(rejcuts[icut])
	  _rej->Fill(fcut);
	else
	  break;
      } 
      // histogram values that pass final Kalman fit
      bool helixfit = (hsf__value&_tffval[helixOK])>0;
      bool seedfit = (ksf__value&_tffval[seedOK])>0;
      bool kalfit = (kff__value&_tffval[kalmanOK])>0;
      if(helixfit){
	_hn->Fill(hsn);
	_hna->Fill(hsna);
	if(seedfit){
	  _shn->Fill(hsn);
	  _shna->Fill(hsna);
	  if(kalfit){
	    _fhn->Fill(hsn);
	    _fhna->Fill(hsna);
	  }
	}
      }
   }
   // normalize
   if(_norm <= 0.0 )
     _norm = _eff->GetBinContent(1);
   cout << "Normalizing efficiency to " << _norm << endl;
   _eff->Scale(1.0/_norm);
   cout << "Normalizing rejection to " << nentries << endl;
   _rej->Scale(1.0/nentries);
}

void TrkRecoEff::drawHistos(){
  if(_cuts == 0)
    _cuts =new TCanvas("cuts","Cut Values",1200,800);
  _cuts->Clear();
  _cuts->Divide(2,2);
  _cuts->cd(1);
  _hn->Draw();
  _shn->Draw("same");
  _fhn->Draw("same");
  _cuts->cd(2);
  _hna->Draw();
  _shna->Draw("same");
  _fhna->Draw("same");


}

void TrkRecoEff::createHistos() {
  cout << "in createHistos" << endl;
  // define the bin labels for efficiency plot
  std::vector<string> effbins, rejbins;
  effbins.push_back("All Events");rejbins.push_back(effbins.back());
  effbins.push_back("MC NDigis");
  effbins.push_back("MC t0");
  effbins.push_back("MC Pitch");
  effbins.push_back("Time Cluster");rejbins.push_back(effbins.back());
  effbins.push_back("Circle Fit");
  effbins.push_back("Helix Fit");rejbins.push_back(effbins.back());
  effbins.push_back("Seed Fit");rejbins.push_back(effbins.back());
  effbins.push_back("Kalman Fit");rejbins.push_back(effbins.back());

  _eff = new TH1F("eff","Efficiency",effbins.size(),-0.5,effbins.size()-0.5);
  for(unsigned ibin=0; ibin< effbins.size(); ++ibin){
    _eff->GetXaxis()->SetBinLabel(ibin+1,effbins[ibin].c_str());
  }

  _rej = new TH1F("rej","Rejection",rejbins.size(),-0.5,rejbins.size()-0.5);
  for(unsigned ibin=0; ibin< rejbins.size(); ++ibin){
    _rej->GetXaxis()->SetBinLabel(ibin+1,rejbins[ibin].c_str());
  }

  _eff->SetStats(0);
  _rej->SetStats(0);
  _eff->SetMarkerStyle(20);
  _rej->SetMarkerStyle(20);
  _eff->Sumw2();
  _rej->Sumw2();

  _hn = new TH1F("hn","N Helix Hits",101,-0.5,100.5);
  _shn = new TH1F("shn","N Helix Hits",101,-0.5,100.5);
  _fhn = new TH1F("fhn","N Helix Hits",101,-0.5,100.5);

  _hna = new TH1F("hna","N Helix Active Hits",101,-0.5,100.5);
  _shna = new TH1F("shna","N Helix Active Hits",101,-0.5,100.5);
  _fhna = new TH1F("fhna","N Helix Active Hits",101,-0.5,100.5);

  _hn->SetLineColor(kBlack);
  _hna->SetLineColor(kBlack);

  _shn->SetLineColor(kBlue);
  _shna->SetLineColor(kBlue);

  _fhn->SetLineColor(kRed);
  _fhna->SetLineColor(kRed);
}
