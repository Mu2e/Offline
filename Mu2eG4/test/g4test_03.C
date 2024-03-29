//
// Root script to display some of the plots in made by ReadBack_plugin.cc
//
//
// Original author Rob Kutschke
//
// This also serves as a introduction to root scripts for those who
// have not seen root.  It shows how to:
//
// 1) Retrieve histograms and ntuples from the file that was created
//    by g4test_03.py.
//
// 2) Draw the histograms to the screen (called a canvas).
//
// 3) Split a canvas into multiple pads and draw a different histogram
//    in each pad.
//
// 4) Project scatter plots from ntuples.
//
// 5) Save the canvas in a format suitable for printing ( postscript )
//    or in a format suitable for inclusion in other documents
//    ( png, jpg, gif ).
//
// 6) Root documentation is available at: http://root.cern.ch
//
//

{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  gStyle->SetOptStat("emruo");

  // Base name of input file and of all plot files.
  TString filebasename("g4test_03");

  // Open the input file that contains histograms and ntuples
  // made by ReadBack.cc
  TFile* file = new TFile( filebasename + ".root");

  // Name of the output postscript file.
  // Postscript is the only graphics format for which root supports multi-page output files.
  TString psfile( filebasename + ".pdf");

  // Get pointers to some of the histograms.
  TH1F* hMultiplicity;  file->GetObject("checkhits/hMultiplicity",  hMultiplicity);
  TH1F* hHitNeighbours; file->GetObject("checkhits/hHitNeighbours", hHitNeighbours);
  TH1F* hEnergyDep;     file->GetObject("checkhits/hEnergyDep",     hEnergyDep);
  TH1F* hStepLength;    file->GetObject("checkhits/hStepLength",    hStepLength);

  // Get a pointer to the ntuple.
  TNtuple* nt; file->GetObject("checkhits/ntup",nt);

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + filebasename, 900, 900 );

  // Open a multi-page output postscript file .
  canvas->Print( psfile+"[");

  // Split the canvas into 4 pads.
  canvas->Divide(2,2);

  // Draw some histograms, one per pad.
  // cd(n): move to graphics pad number "n".
  // "H9": draw outline histogram ("H") in high resolution mode (9)
  canvas->cd(1); hMultiplicity->Draw("H9");
  canvas->cd(2); hHitNeighbours->Draw("H9");
  canvas->cd(3); hEnergyDep->Draw("H9");
  canvas->cd(4); hStepLength->Draw("H9");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  //  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a png file (slow)
  canvas->Print( filebasename + "_1.png" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for page 2.
  canvas->cd(0);
  canvas->Clear();

  // Draw a y vs x scatter plot of the hit positions.
  // Let root choose all options.
  nt->Draw( "hx:hy","","");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  //  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  canvas->Print( filebasename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for page 3.
  canvas->cd(0);
  canvas->Clear();

  // Make a scatterplot of Energy deposition vs step length.
  // Override the root supplied titles.
  // DrawFrame: arguments are xmin, ymin, xmax, ymax
  // SetTitle:  semi-colon separates "Main title; x-axis title; y-axis title"
  // PSAME:     plot as points (P), overlay on existing plot ("SAME").
  TH1F* frame = canvas->DrawFrame(0., 0., 10., 15.);
  frame->SetTitle("Energy Deposition vs Step Length in Cell;(mm);(keV)");
  nt->Draw( "edep:step","","PSAME");

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment these lines to save this canvas as a gif file and as a pdf file.
  // canvas->Print( filebasename + "_3.gif" );
  // canvas->Print( filebasename + "_3.pdf" );

  // Close the postscript file.
  canvas->Print(psfile+"]");

}
