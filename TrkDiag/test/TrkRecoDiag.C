#define TrkRecoDiag_cxx
#include "TrkRecoDiag.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

TrkRecoDiag::TrkRecoDiag(TTree *tree, double norm) : fChain(0) ,_tffval(32,0), _norm(norm),
  _eff(0), _acc(0), _effcan(0), _tcseln(6), 
  _hseln(6), _hselminm(280.0),_hselmaxm(380.0),
  _sselminm(95.0), _sselmaxm(110.0), _sselmerr(1.2), _sselchi2(9.0), 
  _pseltq(0.4), _pselminm(95.0), _pselmaxm(110.0),
  _mincost(0.45), _maxcost(0.7571),
  _mint0(700.0),
  _tdlow(0.57735027), _tdhigh(1.0)
 
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


void TrkRecoDiag::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::vector<bool> effcuts, acccuts;
      effcuts.push_back(true); acccuts.push_back(effcuts.back());// 0 bin counts number of events
      effcuts.push_back(ndigi >= 15); // events with at least 15 digis from the primary particle
      // take time modulo MB period
      float mctime = fmod(mcmidt0,1695.0);
      effcuts.push_back(mctime > 500.0 ); // primary track passed the tracker during 
      float  mccost = mcmidpz/mcmidmom;
      float hd0 = hsh__rcent - hsh__radius;
      float hrmax  = hsh__rcent + hsh__radius;
      float hmom = sqrt(hsh__radius*hsh__radius + hsh__lambda*hsh__lambda);
      effcuts.push_back(mccost > _mincost && mccost < _maxcost ); // track direction (with buffer)
      effcuts.push_back(tcn>_tcseln); acccuts.push_back(effcuts.back());// reconstructed time cluster associated with primary track
      effcuts.push_back((hsf__value&_tffval[circleOK])>0);
      effcuts.push_back((hsf__value&_tffval[helixOK])>0);acccuts.push_back(effcuts.back());
      effcuts.push_back(hsna>_hseln && hmom < _hselmaxm && hmom > _hselminm);acccuts.push_back(effcuts.back());
      effcuts.push_back((ksf__value&_tffval[seedOK])>0);acccuts.push_back(effcuts.back());
      effcuts.push_back(ksm > _sselminm && ksm < _sselmaxm && kschisq/(ksna-5)<_sselchi2 && ksmerr < _sselmerr);acccuts.push_back(effcuts.back());
      effcuts.push_back((kff__value&_tffval[kalmanOK])>0);acccuts.push_back(effcuts.back());
      // physics cuts
      bool ptqual = kftq > _pseltq;
      bool ppitch = kfh__pars[3] > _tdlow && kfh__pars[3] < _tdhigh;
      bool ptime = kft0> _mint0;
      bool pmom = kfm > _pselminm && kfm < _pselmaxm;
      bool pd0 = kfh__pars[0]<105&&kfh__pars[0]>-80;
      bool prmax =  (kfh__pars[0]+2/kfh__pars[2])>450 && (kfh__pars[0]+2/kfh__pars[2])<680;
      effcuts.push_back(ptqual && ppitch && ptime && pmom && pd0 && prmax );acccuts.push_back(effcuts.back());
      // efficiency 
      if((int)effcuts.size() != _eff->GetNbinsX())
	cout << "Error: Efficiency cuts size " << effcuts.size() << " != " << _eff->GetNbinsX() << endl;
      for(unsigned icut=0;icut< effcuts.size(); ++icut){
	if(effcuts[icut])
	  _eff->Fill((Float_t)icut);
	else
	  break;
      }
      // acceptance
      if((int)acccuts.size() != _acc->GetNbinsX())
	cout << "Error: Acceptance cuts size " << acccuts.size() << " != " << _acc->GetNbinsX() << endl;
      for(unsigned icut=0;icut< acccuts.size(); ++icut){
	if(acccuts[icut])
	  _acc->Fill((Float_t)icut);
	else
	  break;
      }

      // histogram values that pass final Kalman fit
      bool helixfit = (hsf__value&_tffval[helixOK])>0;
      bool seedfit = (ksf__value&_tffval[seedOK])>0;
      bool kalfit = (kff__value&_tffval[kalmanOK])>0;
      double kndf = kschisq/(ksna-5);
      bool psel = kfm > _pselminm && kfm < _pselmaxm && kftq > _pseltq;
      if(helixfit){
	_hn->Fill(hsn);
	_hna->Fill(hsna);
	_hd0->Fill(hd0);
	_hrmax->Fill(hrmax);
	_hmom->Fill(hmom);
	_hrad->Fill(hsh__radius);
	_hlam->Fill(hsh__lambda);
	if(seedfit){
	  _shn->Fill(hsn);
	  _shna->Fill(hsna);
	  _shd0->Fill(hd0);
	  _shrmax->Fill(hrmax);
	  _shmom->Fill(hmom);
	  _shrad->Fill(hsh__radius);
	  _shlam->Fill(hsh__lambda);
	  _ssna->Fill(ksna);
	  _ssmom->Fill(ksm);
	  _ssmerr->Fill(ksmerr);
	  _sschisq->Fill(kndf);
	  if(kalfit){
	    _fhn->Fill(hsn);
	    _fhna->Fill(hsna);
	    _fhd0->Fill(hd0);
	    _fhrmax->Fill(hrmax);
	    _fhmom->Fill(hmom);
	    _fhrad->Fill(hsh__radius);
	    _fhlam->Fill(hsh__lambda);
	    _fsna->Fill(ksna);
	    _fsmom->Fill(ksm);
	    _fsmerr->Fill(ksmerr);
	    _fschisq->Fill(kndf);
	    if(psel){
	      _phn->Fill(hsn);
	      _phna->Fill(hsna);
	      _phd0->Fill(hd0);
	      _phrmax->Fill(hrmax);
	      _phmom->Fill(hmom);
	      _phrad->Fill(hsh__radius);
	      _phlam->Fill(hsh__lambda);
	      _psna->Fill(ksna);
	      _psmom->Fill(ksm);
	      _psmerr->Fill(ksmerr);
	      _pschisq->Fill(kndf);
	    }
	  }
	}
      }
   }
   // normalize
   cout << "Normalizing efficiency to " << _norm << endl;
   _eff->Scale(1.0/_norm);
   _acc->Scale(1.0/_norm);
}

void TrkRecoDiag::createHistos() {
  // define the bin labels for efficiency plot
  std::vector<string> effbins, accbins;
  effbins.push_back("All Events");accbins.push_back(effbins.back());
  effbins.push_back("MC NDigis");
  effbins.push_back("MC t0");
  effbins.push_back("MC Pitch");
  effbins.push_back("Time Cluster");accbins.push_back(effbins.back());
  effbins.push_back("Circle Fit");
  effbins.push_back("Helix Fit");accbins.push_back(effbins.back());
  effbins.push_back("Helix Cuts");accbins.push_back(effbins.back());
  effbins.push_back("Seed Fit");accbins.push_back(effbins.back());
  effbins.push_back("Seed Cuts");accbins.push_back(effbins.back());
  effbins.push_back("Kalman Fit");accbins.push_back(effbins.back());
  effbins.push_back("TrkQual Sel");accbins.push_back(effbins.back());

  _eff = new TH1F("eff",title("Efficiency"),effbins.size(),-0.5,effbins.size()-0.5);
  _acc = new TH1F("acc",title("Acceptance"),accbins.size(),-0.5,accbins.size()-0.5);
  for(unsigned ibin=0; ibin< effbins.size(); ++ibin){
    _eff->GetXaxis()->SetBinLabel(ibin+1,effbins[ibin].c_str());
  }
  _eff->SetStats(0);
  _eff->SetMarkerStyle(20);
  _eff->Sumw2();
  _eff->SetMaximum(1.0);
  _eff->SetLineColor(fChain->GetLineColor());
  _eff->SetMarkerColor(fChain->GetLineColor());

  for(unsigned ibin=0; ibin< accbins.size(); ++ibin){
    _acc->GetXaxis()->SetBinLabel(ibin+1,accbins[ibin].c_str());
  }
  _acc->SetStats(0);
  _acc->SetMarkerStyle(20);
  _acc->Sumw2();
  _acc->SetMaximum(1.0);
  _acc->SetLineColor(fChain->GetLineColor());
  _acc->SetMarkerColor(fChain->GetLineColor());

  _hn = new TH1F("hn",title("Helix NHits"),81,-0.5,80.5);
  _shn = new TH1F("shn",title("Helix NHits"),81,-0.5,80.5);
  _fhn = new TH1F("fhn",title("Helix NHits"),81,-0.5,80.5);
  _phn = new TH1F("phn",title("Helix NHits"),81,-0.5,80.5);

  _hna = new TH1F("hna",title("Helix NActive"),61,-0.5,60.5);
  _shna = new TH1F("shna",title("Helix NActive"),61,-0.5,60.5);
  _fhna = new TH1F("fhna",title("Helix NActive"),61,-0.5,60.5);
  _phna = new TH1F("phna",title("Helix NActive"),61,-0.5,60.5);

  _hd0 = new TH1F("hd0",title("Helix d0;mm"),101,-150.0,250.0);
  _shd0 = new TH1F("shd0",title("Helix d0;mm"),101,-150.0,250.0);
  _fhd0 = new TH1F("fhd0",title("Helix d0;mm"),101,-150.0,250.0);
  _phd0 = new TH1F("phd0",title("Helix d0;mm"),101,-150.0,250.0);

  _hrmax = new TH1F("hrmax",title("Helix rmax;mm"),101,400.0,700.0);
  _shrmax = new TH1F("shrmax",title("Helix rmax;mm"),101,400.0,700.0);
  _fhrmax = new TH1F("fhrmax",title("Helix rmax;mm"),101,400.0,700.0);
  _phrmax = new TH1F("phrmax",title("Helix rmax;mm"),101,400.0,700.0);

  _hrad = new TH1F("hrad",title("Helix rad;mm"),101,100.0,400.0);
  _shrad = new TH1F("shrad",title("Helix rad;mm"),101,100.0,400.0);
  _fhrad = new TH1F("fhrad",title("Helix rad;mm"),101,100.0,400.0);
  _phrad = new TH1F("phrad",title("Helix rad;mm"),101,100.0,400.0);

  _hlam = new TH1F("hlambda",title("Helix lambda;mm"),101,100.0,400.0);
  _shlam = new TH1F("shlambda",title("Helix lambda;mm"),101,100.0,400.0);
  _fhlam = new TH1F("fhlambda",title("Helix lambda;mm"),101,100.0,400.0);
  _phlam = new TH1F("phlambda",title("Helix lambda;mm"),101,100.0,400.0);

  _hmom = new TH1F("hmom",title("Helix mom;mm"),101,200.0,500.0);
  _shmom = new TH1F("shmom",title("Helix mom;mm"),101,200.0,500.0);
  _fhmom = new TH1F("fhmom",title("Helix mom;mm"),101,200.0,500.0);
  _phmom = new TH1F("phmom",title("Helix mom;mm"),101,200.0,500.0);

  _hn->SetLineColor(kBlack);
  _hna->SetLineColor(kBlack);
  _hd0->SetLineColor(kBlack);
  _hrmax->SetLineColor(kBlack);
  _hrad->SetLineColor(kBlack);
  _hlam->SetLineColor(kBlack);
  _hmom->SetLineColor(kBlack);

  _shn->SetLineColor(kBlue);
  _shna->SetLineColor(kBlue);
  _shd0->SetLineColor(kBlue);
  _shrmax->SetLineColor(kBlue);
  _shrad->SetLineColor(kBlue);
  _shlam->SetLineColor(kBlue);
  _shmom->SetLineColor(kBlue);

  _fhn->SetLineColor(kRed);
  _fhna->SetLineColor(kRed);
  _fhd0->SetLineColor(kRed);
  _fhrmax->SetLineColor(kRed);
  _fhrad->SetLineColor(kRed);
  _fhlam->SetLineColor(kRed);
  _fhmom->SetLineColor(kRed);

  _phn->SetLineColor(kGreen);
  _phna->SetLineColor(kGreen);
  _phd0->SetLineColor(kGreen);
  _phrmax->SetLineColor(kGreen);
  _phrad->SetLineColor(kGreen);
  _phlam->SetLineColor(kGreen);
  _phmom->SetLineColor(kGreen);

  _ssna = new TH1F("ssna",title("Seed NActive"),61,-0.5,60.5);
  _fsna = new TH1F("fsna",title("Seed NActive"),61,-0.5,60.5);
  _psna = new TH1F("psna",title("Seed NActive"),61,-0.5,60.5);

  _ssmom = new TH1F("ssmom",title("Seed mom;MeV/c"),101,85.0,120.0);
  _fsmom = new TH1F("fsmom",title("Seed mom;MeV/c"),101,85.0,120.0);
  _psmom = new TH1F("psmom",title("Seed mom;MeV/c"),101,85.0,120.0);

  _ssmerr = new TH1F("ssmerr",title("Seed mom err;MeV/c"),101,0.0,2.0);
  _fsmerr = new TH1F("fsmerr",title("Seed mom err;MeV/c"),101,0.0,2.0);
  _psmerr = new TH1F("psmerr",title("Seed mom err;MeV/c"),101,0.0,2.0);

  _sschisq = new TH1F("sschisq",title("Seed chisq/dof"),101,0.0,15.0);
  _fschisq = new TH1F("fschisq",title("Seed chisq/dof"),101,0.0,15.0);
  _pschisq = new TH1F("pschisq",title("Seed chisq/dop"),101,0.0,15.0);
  
  _ssna->SetLineColor(kBlue);
  _ssmom->SetLineColor(kBlue);
  _ssmerr->SetLineColor(kBlue);
  _sschisq->SetLineColor(kBlue);

  _fsna->SetLineColor(kRed);
  _fsmom->SetLineColor(kRed);
  _fsmerr->SetLineColor(kRed);
  _fschisq->SetLineColor(kRed);

  _psna->SetLineColor(kGreen);
  _psmom->SetLineColor(kGreen);
  _psmerr->SetLineColor(kGreen);
  _pschisq->SetLineColor(kGreen);

}

const char* TrkRecoDiag::title(const char* titl) {
  static string stitle;
  stitle = string(fChain->GetTitle());
  stitle += string(titl);
  return stitle.c_str();
}
