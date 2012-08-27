//
// Root c++ function to test MakeStrawHit_module
// 
// $Id: strawHits.C,v 1.4 2012/08/27 22:31:40 genser Exp $
// $Author: genser $
// $Date: 2012/08/27 22:31:40 $
// 
// Original author KLG somewat based on Rob Kutschke's example
//
// 1) Retrieve histograms and ntuples from the file that was created
//    by mixExample01.fcl
//
// 2) Draw the histograms to the screen (called a canvas).
//
// 3) Split a canvas into multiple pads and draw a different histogram
//    in each pad.
//
// 4) Save the canvas in a format suitable for printing ( postscript )
//    or in a format suitable for inclusion in other documents
//    ( png, jpg, gif ).
//

// run it in root e.g. like .x strawHits.C++ (or .x strawHits.C++g)

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>

#include <TH1.h>
#include <TStyle.h>
#include <TString.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TTree.h>

#include <sstream>
#include <iostream>
#include <vector>

void strawHits()
{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  gStyle->SetOptStat("emruo");

  // Base name of input file and of all plot files.
  TString basename("strawHits_01");

  // Open the input file that contains histograms and ntuples
  // made by ReadStrawHits_module

  std::vector<TFile*>  file;
  std::vector<TPaveLabel*> fileLabel;

  file.push_back(new TFile("mixExample01_2000_before_xtalk_off.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"before","NDC"));
  file.push_back(new TFile("mixExample01_2000_after_xtalk_off_e2a_off_save.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"after xtalk off e2a off","NDC"));
  file.push_back(new TFile("mixExample01_2000_after_xtalk_on2p_e2a_off.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"after xtalk on 2\% e2a off","NDC"));
  file.push_back(new TFile("mixExample01_2000_after_xtalk_off_e2a_on.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"after xtalk off e2a on","NDC"));
  file.push_back(new TFile("mixExample01_2000_after_xtalk_on2p_e2a_on.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"after xtalk on 2\% e2a on","NDC"));
  file.push_back(new TFile("mixExample01_2000_after_xtalk_on1p_e2a_on.root"));
  fileLabel.push_back(new TPaveLabel(0.505,0.905,0.78,0.955,"after xtalk on 1\% e2a on","NDC"));

  const int nfiles = file.size();

  for (int ii=0; ii!=nfiles; ++ii) {
    fileLabel[ii]->SetBorderSize(0);
    fileLabel[ii]->SetFillColor(0);
  }

  // Name of the output pdf file.
  // Postscript is the only graphics format for which root supports multi-page output files.
  TString psfile( basename + ".pdf");

  std::vector<TH1F*> _hHitTime(nfiles);      
  std::vector<TH1F*> _hHitDeltaTime(nfiles); 
  std::vector<TH1F*> _hHitEnergy(nfiles);    
  std::vector<TH1F*> _hNHits(nfiles);        
  std::vector<TH1F*> _hNHitsPerWire(nfiles); 
  std::vector<TH1F*> _hDriftTime(nfiles);    
  std::vector<TH1F*> _hDriftDistance(nfiles);
  std::vector<TH1F*> _hDistanceToMid(nfiles);
  std::vector<TH1F*> _hNG4Steps(nfiles);     
  std::vector<TH1F*> _hG4StepLength(nfiles); 
  std::vector<TH1F*> _hG4StepEdep(nfiles);   
  TH1F* _tmp = 0;

  for (int ii=0; ii!=nfiles; ++ii) {
    file[ii]->GetObject("readStrawHits/hHitTime",       _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hHitDeltaTime",  _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitDeltaTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hHitEnergy",     _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitEnergy[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNHits",         _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHits[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNHitsPerWire",  _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHitsPerWire[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDriftTime",     _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDriftDistance", _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftDistance[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDistanceToMid", _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDistanceToMid[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNG4Steps",      _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNG4Steps[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hG4StepLength",  _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepLength[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hG4StepEdep",    _tmp); 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepEdep[ii]=_tmp; 
  }

  const int nhistograms = 11;
  std::vector<std::vector<TH1F*>*> _histograms(nhistograms);

  _histograms[ 0] = &_hHitTime;      // *_histograms[ 0] is _hHitTime;
  _histograms[ 1] = &_hHitDeltaTime;
  _histograms[ 2] = &_hHitEnergy;
  _histograms[ 3] = &_hNHits;
  _histograms[ 4] = &_hNHitsPerWire;
  _histograms[ 5] = &_hDriftTime;
  _histograms[ 6] = &_hDriftDistance;
  _histograms[ 7] = &_hDistanceToMid;
  _histograms[ 8] = &_hNG4Steps;
  _histograms[ 9] = &_hG4StepLength;
  _histograms[10] = &_hG4StepEdep;

  // Get a pointer to the ntuple.
  // TNtuple* nt; file]->GetObject("readStrawHit/ntup",nt);

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + basename, 900, 900 );

  // Open a multi-page output pdf file .
  canvas->Print( psfile+"[");

  for (int jj=0; jj!=nhistograms; ++jj) {

    // Clear canvas in preparation for next page.
    // Split the canvas into 6 pads.
    canvas->cd(0);
    canvas->Clear();
    canvas->Divide(2,3);

    // Draw some histograms, one per pad.
    // cd(n): move to graphics pad number "n".
    // "H9": draw outline histogram ("H") in high resolution mode (9)

    for (int ii=0; ii!=nfiles; ++ii) {
      canvas->cd(ii+1);
      cout << jj << " " << ii << " Drawing " << ((*_histograms[jj])[ii])->GetTitle() 
           << ", " <<fileLabel[ii]->GetLabel() <<endl;
      ((*_histograms[jj])[ii])->Draw("H9");
      fileLabel[ii]->Draw("9");
    }

    // Flush page to screen
    canvas->Update();

    // Add this canvas to the pdf file.
    canvas->Print(psfile);

    // Prompt and wait for response before continuing.
    cerr << "Double click in the last active pad to continue: " ;
    gPad->WaitPrimitive();
    cerr << endl;

  }

  // Uncomment this line to save this canvas as a png file (slow)
  //canvas->Print( basename + "_1.png" );

  // Close the pdf file.
  canvas->Print(psfile+"]");

}

// This tells emacs to view this file in c++ mode.
// Local Variables:
// mode:c++
// End:
