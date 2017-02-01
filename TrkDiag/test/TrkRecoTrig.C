#define TrkRecoTrig_cxx
#include "TrkRecoTrig.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

TrkRecoTrig::TrkRecoTrig(TTree *tree) : fChain(0) ,_tffval(32,0), _rej(0), _trigcan(0)
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


void TrkRecoTrig::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::vector<bool> rejcuts;
      rejcuts.push_back(true); // 0 bin counts number of events
      rejcuts.push_back(tcn > 10 );
      rejcuts.push_back((hsf__value&_tffval[helixOK])>0);
      rejcuts.push_back((ksf__value&_tffval[seedOK])>0);
      rejcuts.push_back((kff__value&_tffval[kalmanOK])>0);

      // rejection cuts are cumulative
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
   cout << "Normalizing rejection to " << nentries << endl;
   _rej->Scale(1.0/nentries);
}

void TrkRecoTrig::drawHistos(){
  if(_trigcan == 0)
    _trigcan =new TCanvas("trigger","Trigger Rejection",1200,800);
  _trigcan->Clear();
  _trigcan->Divide(2,2);
  _trigcan->cd(1);
  _hn->Draw();
  _shn->Draw("same");
  _fhn->Draw("same");
  _trigcan->cd(2);
  _hna->Draw();
  _shna->Draw("same");
  _fhna->Draw("same");


}

void TrkRecoTrig::createHistos() {
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

  _rej = new TH1F("rej","Rejection",rejbins.size(),-0.5,rejbins.size()-0.5);
  for(unsigned ibin=0; ibin< rejbins.size(); ++ibin){
    _rej->GetXaxis()->SetBinLabel(ibin+1,rejbins[ibin].c_str());
  }

  _rej->SetStats(0);
  _rej->SetMarkerStyle(21);
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
