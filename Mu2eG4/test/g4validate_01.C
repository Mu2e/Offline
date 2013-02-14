//
// Root c++ function to compare tracking plots specified by another macro
// 
// $Id: g4validate_01.C,v 1.3 2013/02/14 18:02:12 genser Exp $
// $Author: genser $
// $Date: 2013/02/14 18:02:12 $
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

  TString histTmpName;
  ostringstream histTmpSuffix;
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

  for (unsigned ff=0; ff!=nfiles; ++ff) {
    fileLabel[ff]->SetBorderSize(0);
    fileLabel[ff]->SetFillColor(0);
    fileLabel[ff]->SetTextColor(602);
    basename = basename + "_" + fileText[ff];
  }

  // Name of the output pdf file.
  // Pdf/Postscript is the only graphics format for which root supports multi-page output files.
  TString pdffile( basename + ".pdf");

  std::vector<std::vector<TH1*> > histograms(nfiles);
  std::vector<TH1*> histograms_copy(nfiles);

  unsigned nhistograms(0);

  for (unsigned ff=0; ff!=nfiles; ++ff) {

    // we fetch histograms from a file

    TrkFitDiag(files[ff],histograms[ff]);

    if ( ff==0 ) {

      nhistograms = histograms[ff].size();

    } else {

      if ( nhistograms!=histograms[ff].size()) {

        cout << "G4validate : number of histograms in file " 
             << files[ff]->GetName() << histograms[ff].size()
             << " is different from initial " << nhistograms << endl;
      }

    }

    // rename/append a file related suffix to the histograms so that
    // each has a unique name (not the title)

    for (unsigned hh=0; hh!=nhistograms; ++hh) {

      histTmpSuffix.str("");
      histTmpSuffix << "_";
      histTmpSuffix.width(1); // we assume the number of files is <9
      histTmpSuffix << ff+1;
      histTmpName = histograms[ff][hh]->GetName()+histTmpSuffix.str();
      histograms[ff][hh]->SetName(histTmpName);
    
    }

  }

  cout << " nhistograms " << nhistograms << endl;

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " 
                                + basename + " for various versions of Geant4", 900, 900 );

  // Open a multi-page output pdf file .
  canvas->Print( pdffile+"[");

  // log scale?
  Int_t islog(0);

  TLegend* leg = new TLegend();

  unsigned canvasCounter(0);

  Float_t vscalemax = 0.;

  for (unsigned hh=0; hh!=nhistograms; ++hh) {

    // Clear canvas in preparation for next page.
    // Split the canvas into n pads.
    canvas->cd(0);
    canvas->Clear();
    canvas->Divide(2,2);

    // Draw some histograms, one per pad.
    // cd(n): move to graphics pad number "n".
    // "H9": draw outline histogram ("H") in high resolution mode (9)
    
    islog = (hh==1000) ? 0 : 1; // decide on the log scale; could be per histogram if needed

    delete leg;
    leg = new TLegend(0.70,0.90,0.90,1.00);
    leg->SetFillColor(0);

    // drawing original histograms from each file

    // first get their maxima
    vscalemax = 0;
    for (unsigned ff=0; ff!=nfiles; ++ff) {
      if ( vscalemax < histograms[ff][hh]->GetMaximum() ) {
        vscalemax = histograms[ff][hh]->GetMaximum();
      }
    }

    cout << "vscalemax = " << vscalemax << endl;
    vscalemax = islog ? vscalemax*2.0 : vscalemax*1.1;

    for (unsigned ff=0; ff!=nfiles; ++ff) {

      canvas->cd(ff+1);
      gPad->SetLogy(islog);
      cout << hh << " " << ff << " Drawing " << histograms[ff][hh]->GetTitle() 
           << ", " << fileLabel[ff]->GetLabel() 
           << ", " << fileText[ff]
           <<endl;
      if (ff==0) 
        {histograms[ff][hh]->SetLineColor(602); }
      else 
        {histograms[ff][hh]->SetLineColor(ff+1);}
      gStyle->SetOptStat("neMRuo");
      
      // collecting legend info
      leg->AddEntry(histograms[ff][hh],fileText[ff],"L");
      histograms[ff][hh]->SetMaximum(vscalemax);
      histograms[ff][hh]->Draw("H9");
      fileLabel[ff]->Draw("9");
      
    }

    // now draw the histograms overlayed + legend
    canvas->cd(nfiles+1);
    for (unsigned ff=0; ff!=nfiles; ++ff) {

      gPad->SetLogy(islog);
      cout << hh << " " << ff << " Drawing " << histograms[ff][hh]->GetTitle() 
           << ", " << fileLabel[ff]->GetLabel() 
           << ", " << fileText[ff]
           <<endl;
      if (ff==0) 
        {histograms[ff][hh]->SetLineColor(602); }
      else 
        {histograms[ff][hh]->SetLineColor(ff+1);}
            
      delete histograms_copy[ff];
      histograms_copy[ff] = static_cast<TH1*>(histograms[ff][hh]->Clone());
      
      histTmpSuffix.str("");
      histTmpSuffix.width(3);
      histTmpSuffix.fill('0');
      histTmpSuffix << ff;
      histTmpName = "hcopy_"+histTmpSuffix.str();
      histograms_copy[ff]->SetName(histTmpName);
      histograms_copy[ff]->SetStats(kFALSE);
      histograms_copy[ff]->SetMaximum(vscalemax);

      if (ff==0) {
        (histograms_copy[ff])->Draw("H9E");
      } else {
        (histograms_copy[ff])->Draw("H9ESAME");
        (histograms_copy[ff])->KolmogorovTest((histograms_copy[ff-1]),"UNOD");
      }
      if (ff==nfiles-1) {
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
