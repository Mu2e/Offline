//
// Root c++ function to test MakeStrawHit_module
// 
// $Id: strawHitsO.C,v 1.1 2013/01/30 01:13:44 genser Exp $
// $Author: genser $
// $Date: 2013/01/30 01:13:44 $
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
#include <TNtuple.h>
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

void strawHitsO()
{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  gStyle->SetOptStat("emruo");
  // gStyle->SetOptStat(kFALSE);

  // Base name of input file and of all plot files.
  TString basename("steps_strawHits_g4942_g4952");


  TString histtmpName;
  ostringstream histtmpSuffix;
  ostringstream canvasSuffix;
  ostringstream histIdos;

  // Open the input file that contains histograms and ntuples
  // made by ReadStrawHits_module

  std::vector<TFile*>  file;
  std::vector<TPaveLabel*> fileLabel;
  std::vector<char *> fileText;

//   file.push_back(new TFile("g4test_03.g496.ftfpberthp.20121219140908.root"));
//   fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g496 FTFP_BERT_HP 5000","NDC"));

//   file.push_back(new TFile("g4test_03.g4952.20121217183217.root"));
//   fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4952 QGSP_BERT_HP","NDC"));

//   file.push_back(new TFile("g4test_03.g496.20121217182923.root"));
//   fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g496 QGSP_BERT_HP","NDC"));

//   file.push_back(new TFile("g4test_03.g496.ftfpberthp.20121218124143.root"));
//   fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g496 FTFP_BERT_HP","NDC"));

  file.push_back(new TFile("g4test_03.g4942.qgspberthp.20130129144528.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4942 FTFP_BERT_HP","NDC"));
  fileText.push_back("g4942");

  file.push_back(new TFile("g4test_03.g4952.qgspberthp.20130129144405.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4952 FTFP_BERT_HP","NDC"));
  fileText.push_back("g4952");

  const int nfiles = file.size();

  for (int ii=0; ii!=nfiles; ++ii) {
    fileLabel[ii]->SetBorderSize(0);
    fileLabel[ii]->SetFillColor(0);
  }

  // Name of the output pdf file.
  // Postscript is the only graphics format for which root supports multi-page output files.
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

  // Get  pointers to the ntuples
  //  TNtuple* _nt;  file->GetObject("checkhits/ntup",nt);
  //  TNtuple* _snt; file->GetObject("readStrawHits/ntup",snt);

  TH1F*    _tmp = 0;
  TNtuple* _tmpnt = 0;

  int nhistc(0);
  int nntc(0);

  for (int ii=0; ii!=nfiles; ++ii) {

    file[ii]->GetObject("checkhits/hMultiplicity",  _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hMultiplicity[ii]=_tmp;
    file[ii]->GetObject("checkhits/hHitNeighbours", _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitNeighbours[ii]=_tmp;
    file[ii]->GetObject("checkhits/hEnergyDep",     _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hEnergyDep[ii]=_tmp;
    file[ii]->GetObject("checkhits/hStepLength",    _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hStepLength[ii]=_tmp;
    file[ii]->GetObject("checkhits/hRadius",        _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hRadius[ii]=_tmp;
    file[ii]->GetObject("checkhits/hxHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hxHit[ii]=_tmp;
    file[ii]->GetObject("checkhits/hyHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hyHit[ii]=_tmp;
    file[ii]->GetObject("checkhits/hzHit",          _tmp); if (ii==0) {++nhistc;}
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hzHit[ii]=_tmp;

    file[ii]->GetObject("readStrawHits/hHitTime",       _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hHitDeltaTime",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitDeltaTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hHitEnergy",     _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hHitEnergy[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNHits",         _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHits[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNHitsPerWire",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNHitsPerWire[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDriftTime",     _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftTime[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDriftDistance", _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDriftDistance[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hDistanceToMid", _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hDistanceToMid[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hNG4Steps",      _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hNG4Steps[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hG4StepLength",  _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepLength[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hG4StepRelTimes",    _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepRelTimes[ii]=_tmp; 
    file[ii]->GetObject("readStrawHits/hG4StepEdep",    _tmp); if (ii==0) {++nhistc;} 
    if (_tmp==0) {cerr << "missing histogram" <<endl; return;} _hG4StepEdep[ii]=_tmp; 

    file[ii]->GetObject("checkhits/ntup",    _tmpnt); if (ii==0) {++nntc;} 
    if (_tmpnt==0) {cerr << "missing ntuple" <<endl; return;} _nt[ii]=_tmpnt; 

    file[ii]->GetObject("readStrawHits/ntup",    _tmpnt); if (ii==0) {++nntc;} 
    if (_tmpnt==0) {cerr << "missing ntuple" <<endl; return;} _snt[ii]=_tmpnt; 

  }

  const int nhistograms(nhistc);
  std::vector<std::vector<TH1F*>*> _histograms(nhistograms);
  std::vector<TH1F*> _histograms_copy(nfiles);

  const int nntuples(nntc);
  std::vector<std::vector<TNtuple*>*> _ntuples(nntuples);

  cout << " nhistograms " << nhistc << endl;

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

  cout << " nhistograms " << nhistc << endl;

  // Get a pointer to the ntuple.
  // TNtuple* nt; file->GetObject("readStrawHit/ntup",nt);

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + basename + " for various versions of Geant4", 900, 900 );

  // Open a multi-page output pdf file .
  canvas->Print( pdffile+"[");

  int islog =0;

  TLegend* leg = new TLegend();

  int canvasCounter(0);

  //  for (int jj=0; jj!=0; ++jj) {
  for (int jj=0; jj!=nhistograms; ++jj) {

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

    for (int ii=0; ii!=nfiles; ++ii) {

      canvas->cd(ii+1);
      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << ((*_histograms[jj])[ii])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      ((*_histograms[jj])[ii])->SetLineColor(ii+1);
      gStyle->SetOptStat("neMRuo");
      
      // collecting legend info
      leg->AddEntry((*_histograms[jj])[ii],fileText[ii],"L");

      ((*_histograms[jj])[ii])->Draw("H9");
      fileLabel[ii]->Draw("9");

    }

    // now draw the histograms overlayed + legend
    canvas->cd(nfiles+1);
    for (int ii=0; ii!=nfiles; ++ii) {

      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << ((*_histograms[jj])[ii])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      ((*_histograms[jj])[ii])->SetLineColor(ii+1);
      
      // FIXME need to copy the hist to not to draw stats here...
       delete _histograms_copy[ii];
       _histograms_copy[ii] = static_cast<TH1F*>(((*_histograms[jj])[ii])->Clone());
       //       _histograms_copy[ii] = (TH1F*)((*_histograms[jj])[ii])->Clone();
      
      histtmpSuffix.str("");
      histtmpSuffix.width(3);
      histtmpSuffix.fill('0');
      histtmpSuffix << ii;
      histtmpName = "hcopy_"+histtmpSuffix.str();
      (_histograms_copy[ii])->SetName(histtmpName);
      (_histograms_copy[ii])->SetStats(kFALSE);
      if (ii==0) {
        //((*_histograms[jj])[ii])->Draw("H9");
        (_histograms_copy[ii])->Draw("H9");
      } else {
        //((*_histograms[jj])[ii])->Draw("H9SAME");
        (_histograms_copy[ii])->Draw("H9SAME");
         leg->Draw("9SAME");
      }
      if (ii==nfiles-1) {
         leg->Draw("9SAME");
      }

    }

    // Flush page to screen
    canvas->Update();

    // Uncomment this line to save this canvas as a png file (slow)
    canvasSuffix.str("");
    canvasSuffix.width(3);
    canvasSuffix.fill('0');
    canvasSuffix << ++canvasCounter;

    canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );
    
    // Add this canvas to the pdf file.
    canvas->Print(pdffile);

    // Prompt and wait for response before continuing.
    cerr << "Double click in the last active pad to continue: " ;
    gPad->WaitPrimitive();
    cerr << endl;

  }

  // a plot (scatterplot) based on the ntuple

  canvas->cd(0);
  canvas->Clear();
  canvas->Divide(2,2);

  // Draw some histograms, one per pad.
  // cd(n): move to graphics pad number "n".
  // "H9": draw outline histogram ("H") in high resolution mode (9)
    
  for (int ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+1);
    cout << 0 << " " << ii << " Drawing " << ((*_ntuples[0])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() 
         << ", " << fileText[ii]
         <<endl;
    TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
    frame->SetTitle("StepPoint y vs. x (mm)");
    ((*_ntuples[0])[ii])->Draw( "hx:hy","","PSAME");

    fileLabel[ii]->Draw("9SAME");
  }

  for (int ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+3);
    cout << 1 << " " << ii << " Drawing " << ((*_ntuples[0])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() 
         << ", " << fileText[ii]
         <<endl;
    TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
    frame->SetTitle("StrawHit y vs. x (mm)");
    ((*_ntuples[1])[ii])->Draw( "hitx:hity","","PSAME");

    fileLabel[ii]->Draw("9SAME");
  }

  // Flush page to screen
  canvas->Update();

  // Uncomment this line to save this canvas as a png file (slow)
  canvasSuffix.str("");
  canvasSuffix.width(3);
  canvasSuffix.fill('0');
  canvasSuffix << ++canvasCounter;

  canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );

  // Add this canvas to the pdf file.
  canvas->Print(pdffile);

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // another ntuple plot (hist)

  canvas->cd(0);
  canvas->Clear();
  canvas->Divide(2,2);

  for (int ii=0; ii!=nfiles; ++ii) {
    canvas->cd(ii+1);
    histIdos.str("");
    histIdos << ii+1;
    TString htmpname = "mhg"+histIdos.str();
    TString drawInputString("sqrt(hity*hity+hitx*hitx)>>"+htmpname+"(500,300.,800.)");

    cout << 2 << " " << ii << " Drawing " << ((*_ntuples[0])[ii])->GetTitle() 
         << ", " << fileLabel[ii]->GetLabel() << ", " << drawInputString
         << ", " << fileText[ii]
         <<endl;
    ((*_ntuples[1])[ii])->Draw(drawInputString ,"","");
    TH1F *h = (TH1F*)gDirectory->Get(htmpname);
    h->SetLineColor(ii+1);
    fileLabel[ii]->Draw("SAME");
  }

  // Flush page to screen
  canvas->Update();

  // Uncomment this line to save this canvas as a png file (slow)
  canvasSuffix.str("");
  canvasSuffix.width(3);
  canvasSuffix.fill('0');
  canvasSuffix << ++canvasCounter;

  canvas->Print( basename + "_" + canvasSuffix.str() + ".png" );

  // Add this canvas to the pdf file.
  canvas->Print(pdffile);

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Close the pdf file.
  canvas->Print(pdffile+"]");

}

// This tells emacs to view this file in c++ mode.
// Local Variables:
// mode:c++
// End:
