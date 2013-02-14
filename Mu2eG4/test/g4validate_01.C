//
// Root c++ function to compare tracking plots specified by another macro
// 
// $Id: g4validate_01.C,v 1.1 2013/02/14 00:21:52 genser Exp $
// $Author: genser $
// $Date: 2013/02/14 00:21:52 $
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

// run it in root e.g. like .x g4validate_01.C++ (or .x g4validate_01.C++g)

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

#include "KalmanTests/test/TrkFitDiag.C"

void g4validate_01()
{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  gStyle->SetOptStat("emruo");
  // gStyle->SetOptStat(kFALSE);

  // flag controlling the pause after each canvas 
  //  bool const interactive = true;
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

  files.push_back(new TFile("g4validate_01.g4942.qgspberthp.20130211222638.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4942 QGSP_BERT_HP","NDC"));
  fileText.push_back("g4942QBH");

  files.push_back(new TFile("g4validate_01.g4952.qgspberthp.20130211234925.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4952 QGSP_BERT_HP","NDC"));
  fileText.push_back("g4952QBH");

  files.push_back(new TFile("g4validate_01.g4961.ftfpberthp.20130212083557.root"));
  fileLabel.push_back(new TPaveLabel(0.50,0.90,0.78,0.955,"g4961 FTFP_BERT_HP","NDC"));
  fileText.push_back("g4961FBH");

  // Base name of input file and of all plot files.
  TString basename("ktf");

  const unsigned nfiles = files.size();

  for (unsigned ii=0; ii!=nfiles; ++ii) {
    fileLabel[ii]->SetBorderSize(0);
    fileLabel[ii]->SetFillColor(0);
    fileLabel[ii]->SetTextColor(602);
    basename = basename + "_" + fileText[ii];
  }

  // Name of the output pdf file.
  // Pdf/Postscript is the only graphics format for which root supports multi-page output files.
  TString pdffile( basename + ".pdf");

  std::vector<std::vector<TH1F*> > histograms(nfiles);

  std::vector<TH1F*> histograms_copy(nfiles);

  unsigned nhistograms(0);

  for (unsigned ii=0; ii!=nfiles; ++ii) {

    // we fetch histograms from a file

    TrkFitDiag(files[ii],histograms[ii]);

    if ( ii==0 ) {

      nhistograms = histograms[ii].size();

    } else {

      if ( nhistograms!=histograms[ii].size()) {

        cout << "G4validate : number of histograms in file " << files[ii]->GetName() << histograms[ii].size()
             << " is different from initial " << nhistograms << endl;
      }

    }

  }

  cout << " nhistograms " << nhistograms << endl;

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + basename + " for various versions of Geant4", 900, 900 );

  // Open a multi-page output pdf file .
  canvas->Print( pdffile+"[");

  // log scale?
  Int_t islog(0);

  TLegend* leg = new TLegend();

  unsigned canvasCounter(0);

  Float_t vscalemax = 0.;

  for (unsigned jj=0; jj!=nhistograms; ++jj) {

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
    for (unsigned ii=0; ii!=nfiles; ++ii) {

      canvas->cd(ii+1);
      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << (histograms[ii][jj])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      if (ii==0) 
        {(histograms[ii][jj])->SetLineColor(602); }
      else 
        {(histograms[ii][jj])->SetLineColor(ii+1);}
      gStyle->SetOptStat("neMRuo");
      
      // collecting legend info
      leg->AddEntry(histograms[ii][jj],fileText[ii],"L");

      (histograms[ii][jj])->Draw("H9");
      fileLabel[ii]->Draw("9");
      
      if ( vscalemax < (histograms[ii][jj])->GetMaximum() ) {
        vscalemax = (histograms[ii][jj])->GetMaximum();
      }

    }

    cout << "vscalemax = " << vscalemax << endl;
    vscalemax = islog ? vscalemax*2.0 : vscalemax*1.1;

    // now draw the histograms overlayed + legend
    canvas->cd(nfiles+1);
    for (unsigned ii=0; ii!=nfiles; ++ii) {

      gPad->SetLogy(islog);
      cout << jj << " " << ii << " Drawing " << (histograms[ii][jj])->GetTitle() 
           << ", " << fileLabel[ii]->GetLabel() 
           << ", " << fileText[ii]
           <<endl;
      if (ii==0) 
        {(histograms[ii][jj])->SetLineColor(602); }
      else 
        {(histograms[ii][jj])->SetLineColor(ii+1);}
            
      delete histograms_copy[ii];
      histograms_copy[ii] = static_cast<TH1F*>((histograms[ii][jj])->Clone());
      
      histtmpSuffix.str("");
      histtmpSuffix.width(3);
      histtmpSuffix.fill('0');
      histtmpSuffix << ii;
      histtmpName = "hcopy_"+histtmpSuffix.str();
      histograms_copy[ii]->SetName(histtmpName);
      histograms_copy[ii]->SetStats(kFALSE);
      histograms_copy[ii]->SetMaximum(vscalemax);
      if (ii==0) {
        (histograms_copy[ii])->Draw("H9E");
      } else {
        (histograms_copy[ii])->Draw("H9ESAME");
        (histograms_copy[ii])->KolmogorovTest((histograms_copy[ii-1]),"UNOD");
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

  cout << "Closing " << pdffile
         <<endl;

  // Close the pdf file.
  canvas->Print(pdffile+"]");

}

// This tells emacs to view this file in c++ mode.
// Local Variables:
// mode:c++
// End:
