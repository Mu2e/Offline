//
// Root c++ function to compare plots based on steps, hits, track fits
// 
// $Id: sthktfO.C,v 1.3 2014/04/03 22:26:32 genser Exp $
// $Author: genser $
// $Date: 2014/04/03 22:26:32 $
// 
// Original author KLG somewat based on Rob Kutschke's example
//
// 1) Retrieve histograms, ntuples and tree from the files that was created
//    by g4validate_01.fcl
//
// 2) Draw the histograms to the screen (called a canvas).
//
// 3) Split a canvas into multiple pads and draw a different histogram
//    in each pad.
//
// 4) Save the canvas in a format suitable for printing ( pdf/postscript )
//    or in a format suitable for inclusion in other documents
//    ( png, jpg, gif ).
//
// 5) overlaythe histograms, run Kolmogorov test on them
//

// run it in root e.g. like .x sthktfO.C++ (or .x sthktfO.C++g)

// deprecated
// >>>>>>>>>>>> please use g4validate_01.C with StepsDiag.C instead <<<<<<<<<<<<

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>

#include "TCut.h"

#include <TH1.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TStyle.h>
#include <TString.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLegend.h>

#include <sstream>
#include <iostream>
#include <vector>

#include <iomanip>
#include <fstream>
#include <sstream>

void sthktfO()
{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  // gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);

  // flag controlling the pause after each canvas 
  // bool const interactive = true;
  bool const interactive = false;

  // flag controlling creation of the png files
  //  bool const createpng = true;
  bool const createpng = false;

  TString histtmpName;
  ostringstream histtmpSuffix;
  ostringstream canvasSuffix;
  ostringstream histIdos;

  // Open the input file that contains histograms and ntuples
  // made by ReadStrawHits_module

  std::vector<TFile*>  files;
  std::vector<TPaveLabel*> fileLabel;
  std::vector<TString> fileText;

//    files.push_back(new TFile("g4validate_01_g4961a_20130530135752.root"));
//    fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4961a QGSP_BERT_HP mu2e v3_0_1.34","NDC"));
//    fileText.push_back("g4961aQBHMu2ev3_0_1_34");

//    files.push_back(new TFile("g4validate_01_g4962_20130530121321.root"));
//    fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4962  QGSP_BERT_HP mu2e v3_0_1.34","NDC"));
//    fileText.push_back("g4962QBHMu2ev3_0_1_34");

  files.push_back(new TFile("g4validate_01.g4962.shielding.20130531105554.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4962 Shielding       mu2e v3_0_1.34","NDC"));
  fileText.push_back("g4962Shldv3_0_1_34");

  files.push_back(new TFile("g4validate_01.g4962.shieldingmu2e01.20130531144419.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4962 ShieldingMu2e01 mu2e v3_0_1.34","NDC"));
  fileText.push_back("g4962ShldMu2e01v3_0_1_34");

  files.push_back(new TFile("g4validate_01.g4962.shieldingmu2e00.20130530165103.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4962 ShieldingMu2e00 mu2e v3_0_1.34","NDC"));
  fileText.push_back("g4962ShldMu2e00v3_0_1_34");

  // Base name of input file and of all plot files.
  TString basename("steps_sh_ktf");

  const Int_t nfiles = files.size();

  for (Int_t ii=0; ii!=nfiles; ++ii) {
    fileLabel[ii]->SetBorderSize(0);
    fileLabel[ii]->SetFillColor(0);
    fileLabel[ii]->SetTextColor(602);
    basename = basename + "_" + fileText[ii];
  }

  // Name of the output pdf file.
  // Pdf/Postscript is the only graphics format for which root supports multi-page output files.
  TString pdffile( basename + ".pdf");

  std::vector<TH1F*> _hMultiplicity(nfiles);
  std::vector<TH1F*> _hHitNeighbours(nfiles);
  std::vector<TH1F*> _hEnergyDep(nfiles);
  std::vector<TH1F*> _hStepLength(nfiles);
  std::vector<TH1F*> _hRadius(nfiles);
  std::vector<TH1F*> _hxHit(nfiles);
  std::vector<TH1F*> _hyHit(nfiles);
  std::vector<TH1F*> _hzHit(nfiles);

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
  std::vector<TH1F*> _hG4StepRelTimes(nfiles);   

  std::vector<TNtuple*>  _nt(nfiles);   
  std::vector<TNtuple*> _snt(nfiles);   

  std::vector<TTree*>   _tftree(nfiles);   

  TH1F*    _tmp = 0x0;
  TNtuple* _tmpnt = 0x0;
  TTree*   _tmptree = 0x0;

  Int_t nhistc(0);
  Int_t nntc(0);
  Int_t ntreec(0);

  for (Int_t ii=0; ii!=nfiles; ++ii) {

    // histograms

    files[ii]->GetObject("checkhits/hMultiplicity",  _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hMultiplicity[ii]=_tmp;
    files[ii]->GetObject("checkhits/hHitNeighbours", _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitNeighbours[ii]=_tmp;
    files[ii]->GetObject("checkhits/hEnergyDep",     _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hEnergyDep[ii]=_tmp;
    files[ii]->GetObject("checkhits/hStepLength",    _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hStepLength[ii]=_tmp;
    files[ii]->GetObject("checkhits/hRadius",        _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hRadius[ii]=_tmp;
    files[ii]->GetObject("checkhits/hxHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hxHit[ii]=_tmp;
    files[ii]->GetObject("checkhits/hyHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hyHit[ii]=_tmp;
    files[ii]->GetObject("checkhits/hzHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hzHit[ii]=_tmp;

    files[ii]->GetObject("readStrawHits/hHitTime",       _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitTime[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hHitDeltaTime",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitDeltaTime[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hHitEnergy",     _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitEnergy[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hNHits",         _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHits[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hNHitsPerWire",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHitsPerWire[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hDriftTime",     _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftTime[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hDriftDistance", _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftDistance[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hDistanceToMid", _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDistanceToMid[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hNG4Steps",      _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNG4Steps[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hG4StepLength",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepLength[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hG4StepRelTimes",    _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepRelTimes[ii]=_tmp; 
    files[ii]->GetObject("readStrawHits/hG4StepEdep",    _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepEdep[ii]=_tmp; 

    // ntuples

    files[ii]->GetObject("checkhits/ntup",     _tmpnt); if (ii==0) {++nntc;} 
    if (_tmpnt==0) {cerr << "missing ntuple" <<endl; return;} _nt[ii]=_tmpnt; 

    files[ii]->GetObject("readStrawHits/ntup", _tmpnt); if (ii==0) {++nntc;} 
    if (_tmpnt==0) {cerr << "missing ntuple" <<endl; return;} _snt[ii]=_tmpnt; 

    // trees

    files[ii]->GetObject("RKFDownstreameMinus/trkdiag", _tmptree); if (ii==0) {++ntreec;} 
    if (_tmptree==0) {cerr << "missing tree" <<endl; return;} _tftree[ii]=_tmptree; 

  }

  const Int_t nhistograms(nhistc);
  std::vector<std::vector<TH1F*>*> _histograms(nhistograms);
  std::vector<TH1F*> _histograms_copy(nfiles);
  std::vector<TH1F*> _histograms_orig(nfiles);

  const Int_t nntuples(nntc);
  std::vector<std::vector<TNtuple*>*> _ntuples(nntuples);

  const Int_t ntrees(ntreec);
  std::vector<std::vector<TTree*>*> _trees(ntrees);

  cout << " nhistograms " << nhistc << endl;

  Int_t initialNhist = nhistc;

  nhistc=0;

  _histograms[nhistc++] = &_hMultiplicity;
  _histograms[nhistc++] = &_hHitNeighbours;
  _histograms[nhistc++] = &_hEnergyDep;
  _histograms[nhistc++] = &_hStepLength;
  _histograms[nhistc++] = &_hRadius;
  _histograms[nhistc++] = &_hxHit;
  _histograms[nhistc++] = &_hyHit;
  _histograms[nhistc++] = &_hzHit;

  _histograms[nhistc++] = &_hHitTime;     // *_histograms[...] is _hHitTime;
  _histograms[nhistc++] = &_hHitDeltaTime;
  _histograms[nhistc++] = &_hHitEnergy;
  _histograms[nhistc++] = &_hNHits;
  _histograms[nhistc++] = &_hNHitsPerWire;
  _histograms[nhistc++] = &_hDriftTime;
  _histograms[nhistc++] = &_hDriftDistance;
  _histograms[nhistc++] = &_hDistanceToMid;
  _histograms[nhistc++] = &_hNG4Steps;
  _histograms[nhistc++] = &_hG4StepLength;
  _histograms[nhistc++] = &_hG4StepRelTimes;
  _histograms[nhistc++] = &_hG4StepEdep;

  nntc=0;

  _ntuples[nntc++] = &_nt;
  _ntuples[nntc++] = &_snt;

  ntreec=0;
  _trees[ntreec++] = &_tftree;

  if (initialNhist!=nhistc) {

    cout << " inconsistent number of nhistograms " << nhistc 
         << " expected " << initialNhist 
         << endl;
    return;

  }

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + basename + " for various versions of Geant4", 900, 900 );

  // Open a multi-page output pdf file .
  canvas->Print( pdffile+"[");

  // log scale?
  Int_t islog(0);

  TLegend* leg = new TLegend();

  Int_t canvasCounter(0);

  Float_t vscalemax = 0.;

  for (Int_t jj=0; jj!=nhistograms; ++jj) {

    // Clear canvas in preparation for next page.
    // Split the canvas into n pads.
    canvas->cd(0);
    canvas->Clear();
    canvas->Divide(2,2);

    // Draw some histograms, one per pad.
    // cd(n): move to graphics pad number "n".
    // "H9": draw outline histogram ("H") in high resolution mode (9)
    
    islog = (jj==1000) ? 0 : 1;

    delete leg;
    leg = new TLegend(0.70,0.90,0.90,1.00);
    leg->SetFillColor(0);

    // drawing original hiostograms from each file

    vscalemax = 0;
    for (Int_t ii=0; ii!=nfiles; ++ii) {

      canvas->cd(ii+1);
      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << ((*_histograms[jj])[ii])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      if (ii==0) 
        {((*_histograms[jj])[ii])->SetLineColor(602); }
      else 
        {((*_histograms[jj])[ii])->SetLineColor(ii+1);}
      gStyle->SetOptStat("neMRuo");
      
      // collecting legend info
      leg->AddEntry((*_histograms[jj])[ii],fileText[ii],"L");

      ((*_histograms[jj])[ii])->Draw("H9");
      fileLabel[ii]->Draw("9");
      
      if ( vscalemax < ((*_histograms[jj])[ii])->GetMaximum() ) {
        vscalemax = ((*_histograms[jj])[ii])->GetMaximum();
      }

    }

    cout << "vscalemax = " << vscalemax << endl;
    vscalemax = islog ? vscalemax*2.0 : vscalemax*1.1;

    // now draw the histograms overlayed + legend
    canvas->cd(nfiles+1);
    for (Int_t ii=0; ii!=nfiles; ++ii) {

      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << ((*_histograms[jj])[ii])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      if (ii==0) 
        {((*_histograms[jj])[ii])->SetLineColor(602); }
      else 
        {((*_histograms[jj])[ii])->SetLineColor(ii+1);}
            
      delete _histograms_copy[ii];
      _histograms_copy[ii] = static_cast<TH1F*>(((*_histograms[jj])[ii])->Clone());
      
      histtmpSuffix.str("");
      histtmpSuffix.width(3);
      histtmpSuffix.fill('0');
      histtmpSuffix << ii;
      histtmpName = "hcopy_"+histtmpSuffix.str();
      _histograms_copy[ii]->SetName(histtmpName);
      _histograms_copy[ii]->SetStats(kFALSE);
      _histograms_copy[ii]->SetMaximum(vscalemax);
      if (ii==0) {
        (_histograms_copy[ii])->Draw("H9E");
      } else {
        (_histograms_copy[ii])->Draw("H9ESAME");
        (_histograms_copy[ii])->KolmogorovTest((_histograms_copy[ii-1]),"UNOD");
      }
      if (ii==nfiles-1) {
        leg->Draw("9SAME");
      }

    }

    // Flush page to screen
    canvas->Update();

    if (createpng) {
      // save the canvas to a png file if requested
      canvasSuffix.str("");
      canvasSuffix.width(3);
      canvasSuffix.fill('0');
      canvasSuffix << canvasCounter;
      canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );
      //      canvas->Print( basename + "_" + canvasSuffix.str() + ".root" );
    }

    ++canvasCounter;
    // Add this canvas to the pdf file.
    canvas->Print(pdffile);

    if (interactive) {
      // Prompt and wait for response before continuing.
      cerr << "Double click in the last active pad to continue: " ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

  }

  // a plot based on the ntuples (hist)

  canvas->cd(0);
  canvas->Clear();
  canvas->Divide(2,2);

  vscalemax = 0.;
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+1);
    histIdos.str("");
    histIdos << ii+1;
    TString htmpname = "hitR"+histIdos.str();
    TString drawInputString("sqrt(hity*hity+hitx*hitx)>>"+htmpname+"(500,300.,800.)");

    cout << 0 << " " << ii << " Drawing " << ((*_ntuples[1])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() << ", " << drawInputString
         << ", " << fileText[ii]
         <<endl;
    ((*_ntuples[1])[ii])->Draw(drawInputString ,"","");
    delete _histograms_copy[ii];
    if (ii==0) 
      { (static_cast<TH1F*>(gDirectory->Get(htmpname)))->SetLineColor(602); }
    else
      { (static_cast<TH1F*>(gDirectory->Get(htmpname)))->SetLineColor(ii+1); }
    _histograms_copy[ii] = static_cast<TH1F*>((gDirectory->Get(htmpname))->Clone());
    histtmpSuffix.str("");
    histtmpSuffix.width(3);
    histtmpSuffix.fill('0');
    histtmpSuffix << ii;
    histtmpName = "hcopy_"+histtmpSuffix.str();
    _histograms_copy[ii]->SetName(histtmpName);
    _histograms_copy[ii]->SetStats(kFALSE);
    fileLabel[ii]->Draw("SAME");
    if ( vscalemax < _histograms_copy[ii]->GetMaximum() ) {
      vscalemax = _histograms_copy[ii]->GetMaximum();
    }
  }

  cout << "vscalemax = " << vscalemax << endl;
  vscalemax*=1.1;

  canvas->cd(nfiles+1);
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    cout << 1 << " " << ii << " Drawing " << (_histograms_copy[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel()
         << ", " << fileText[ii]
         <<endl;

    histIdos.str("");
    histIdos << ii+1;
    TString htmpname = "hitR"+histIdos.str();
    static_cast<TH1F*>((gDirectory->Get(htmpname)))->SetMaximum(vscalemax);
    _histograms_copy[ii]->SetMaximum(vscalemax);
    if (ii==0) {
      (_histograms_copy[ii])->Draw("H9E");
    } else {
      (_histograms_copy[ii])->Draw("H9ESAME");
      (_histograms_copy[ii])->KolmogorovTest((_histograms_copy[ii-1]),"UNOD");
    }
    if (ii==nfiles-1) {
      leg->Draw("9SAME");
    }

  }

  // Flush page to screen
  canvas->Update();

  if (createpng) {

    canvasSuffix.str("");
    canvasSuffix.width(3);
    canvasSuffix.fill('0');
    canvasSuffix << canvasCounter;
    canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );
    //    canvas->Print( basename + "_" + canvasSuffix.str() + ".root" );

  }

  // Add this canvas to the pdf file.
  canvas->Print(pdffile);
  ++canvasCounter;

  if (interactive) {
    // Prompt and wait for response before continuing.
    cerr << "Double click in the last active pad to continue: " ;
    gPad->WaitPrimitive();
    cerr << endl;
  }

  // plots based on the track fit tree (based on Dave's code); this will change in future and be "wrapped"

  // setup cuts

  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(720);
  double momlow(103.4);
  double momhigh(104.8);
  unsigned minnhits(20);
  unsigned minnactive(25);
  unsigned maxndif(10);
  double maxt0err(0.9);
  double maxmomerr(0.3);
  double minfitcon(1e-3);
  TCut reco,goodfit,cosmic,rmom,rpitch,livegate;
  TCut goodmc, tpitch, tmom, nmch;
  TCut nacut, t0errcut, momerrcut,fitconcut;
  char cutstring[100];
  snprintf(cutstring,100,"nactive>=%i&&nhits-nactive<=%i",minnactive,maxndif);
  nacut = TCut(cutstring);
  snprintf(cutstring,100,"t0err<%f",maxt0err);
  t0errcut = TCut(cutstring);
  snprintf(cutstring,100,"fit.momerr<%f",maxmomerr);
  momerrcut = TCut(cutstring);
  snprintf(cutstring,100,"fitcon>%f",minfitcon);
  fitconcut = TCut(cutstring);

  snprintf(cutstring,100,"td>%4.3f&&td<%4.3f",tdlow,tdhigh);
  rpitch = TCut(cutstring);
  snprintf(cutstring,100,"t0>%f",t0min);
  livegate = TCut(cutstring);
  snprintf(cutstring,100,"mcent.td>%4.3f&&mcent.td<%4.3f",tdlow-0.02,tdhigh+0.02);
  tpitch = TCut(cutstring);
  tmom = TCut("mcent.mom>100");
  snprintf(cutstring,100,"nchits>=%i",minnhits);
  nmch = TCut(cutstring);
  reco = TCut("fit.status>0");
  cosmic = TCut("abs(d0)<105 && d0+2/om>450 && d0+2/om<680");
  snprintf(cutstring,100,"fit.mom>%f&&fit.mom<%f",momlow,momhigh);
  rmom = TCut(cutstring);

  goodmc = tpitch+tmom+nmch;
  goodfit = reco+nacut+t0errcut+momerrcut+fitconcut;


  canvas->cd(0);
  canvas->Clear();
  canvas->Divide(2,2);


  vscalemax = 0.;
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+1);
    histIdos.str("");
    histIdos << ii+1;
    TString htmpname = "chisq"+histIdos.str();
    delete _histograms_orig[ii];
    _histograms_orig[ii] = new TH1F(htmpname,"KTF Chisq/NDof",100,0.,10.);

    cout << 0 << " " << ii << " Drawing " << ((*_trees[0])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() << ", " << _histograms_orig[ii]->GetTitle()
         << ", " << fileText[ii]
         <<endl;

    ((*_trees[0])[ii])->Project(htmpname,"chisq/ndof");
    _histograms_orig[ii]->Draw();
    cout << 0 << " done"  <<endl;

    delete _histograms_copy[ii];
    if (ii==0) 
      { _histograms_orig[ii]->SetLineColor(602); }
    else
      { _histograms_orig[ii]->SetLineColor(ii+1); }
    _histograms_copy[ii] = static_cast<TH1F*>(_histograms_orig[ii]->Clone());
    histtmpName = "hcopy_"+ htmpname;
    _histograms_copy[ii]->SetName(histtmpName);
    _histograms_copy[ii]->SetStats(kFALSE);
    fileLabel[ii]->Draw("SAME");
    if ( vscalemax < _histograms_copy[ii]->GetMaximum() ) {
      vscalemax = _histograms_copy[ii]->GetMaximum();
    }

  }

  cout << "vscalemax = " << vscalemax << endl;
  vscalemax*=1.1;

  canvas->cd(nfiles+1);
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    cout << 1 << " " << ii << " Drawing " << (_histograms_copy[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel()
         << ", " << fileText[ii]
         <<endl;
    _histograms_copy[ii]->SetMaximum(vscalemax);
    _histograms_orig[ii]->SetMaximum(vscalemax);
    if (ii==0) {
      (_histograms_copy[ii])->Draw("H9E");
    } else {
      (_histograms_copy[ii])->Draw("H9ESAME");
      (_histograms_copy[ii])->KolmogorovTest((_histograms_copy[ii-1]),"UNOD");
    }
    if (ii==nfiles-1) {
      leg->Draw("9SAME");
    }
  }

  // Flush page to screen
  canvas->Update();

  if (createpng) {

    canvasSuffix.str("");
    canvasSuffix.width(3);
    canvasSuffix.fill('0');
    canvasSuffix << canvasCounter;
    canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );
    //    canvas->Print( basename + "_" + canvasSuffix.str() + ".root" );
  }

  ++canvasCounter;
  // Add this canvas to the pdf file.
  canvas->Print(pdffile);

  if (interactive) {
    // Prompt and wait for response before continuing.
    cerr << "Double click in the last active pad to continue: " ;
    gPad->WaitPrimitive();
    cerr << endl;
  }

  canvas->cd(0);
  canvas->Clear();
  canvas->Divide(2,2);

  vscalemax = 0.;
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+1);
    histIdos.str("");
    histIdos << ii+1;
    TString htmpname = "mres"+histIdos.str();
    delete _histograms_orig[ii];
    _histograms_orig[ii] = new TH1F(htmpname,"momentum resolution at tracker entrance;MeV",100,-2.,2.);

    cout << 0 << " " << ii << " Drawing " << ((*_trees[0])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() << ", " << _histograms_orig[ii]->GetTitle()
         << ", " << fileText[ii]
         <<endl;

    ((*_trees[0])[ii])->Project(htmpname,"fit.mom-mcent.mom",goodmc+goodfit);
    _histograms_orig[ii]->Draw();
    cout << 0 << " done"  <<endl;

    delete _histograms_copy[ii];
    if (ii==0) 
      { _histograms_orig[ii]->SetLineColor(602); }
    else
      { _histograms_orig[ii]->SetLineColor(ii+1); }
    _histograms_copy[ii] = static_cast<TH1F*>(_histograms_orig[ii]->Clone());
    histtmpName = "hcopy_"+ htmpname;
    _histograms_copy[ii]->SetName(histtmpName);
    _histograms_copy[ii]->SetStats(kFALSE);
    fileLabel[ii]->Draw("SAME");
    if ( vscalemax < _histograms_copy[ii]->GetMaximum() ) {
      vscalemax = _histograms_copy[ii]->GetMaximum();
    }

  }

  cout << "vscalemax = " << vscalemax << endl;
  vscalemax*=1.1;

  canvas->cd(nfiles+1);
  for (Int_t ii=0; ii!=nfiles; ++ii) {
    cout << 1 << " " << ii << " Drawing " << (_histograms_copy[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel()
         << ", " << fileText[ii]
         <<endl;
    _histograms_orig[ii]->SetMaximum(vscalemax);
    _histograms_copy[ii]->SetMaximum(vscalemax);
    if (ii==0) {
      (_histograms_copy[ii])->Draw("H9E");
    } else {
      (_histograms_copy[ii])->Draw("H9ESAME");
      (_histograms_copy[ii])->KolmogorovTest((_histograms_copy[ii-1]),"UNOD");
    }
    if (ii==nfiles-1) {
      leg->Draw("9SAME");
    }
  }

  // Flush page to screen
  canvas->Update();

  if (createpng) {

    canvasSuffix.str("");
    canvasSuffix.width(3);
    canvasSuffix.fill('0');
    canvasSuffix << canvasCounter;
    canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );
    //    canvas->Print( basename + "_" + canvasSuffix.str() + ".root" );
  }

  ++canvasCounter;
  // Add this canvas to the pdf file.
  canvas->Print(pdffile);

  if (interactive) {
    // Prompt and wait for response before continuing.
    cerr << "Double click in the last active pad to continue: " ;
    gPad->WaitPrimitive();
    cerr << endl;
  }

  cout << "Closing " << pdffile
         <<endl;

  cout << "deprecated" << endl;
  cout << ">>>>>>>>>>>> please use g4validate_01.C with StepsDiag.C instead <<<<<<<<<<<<" << endl;

  // Close the pdf file.
  canvas->Print(pdffile+"]");

}

// This tells emacs to view this file in c++ mode.
// Local Variables:
// mode:c++
// End:
