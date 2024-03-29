//
// Root script to display some of the plots made by ReadBack_module.cc & ReadStrawHit_module.cc
//
//
// Original author Rob Kutschke, extednd by K. Genser for ReadStrawHit_module.cc
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

  TString basename("g4test_03");

  // Open the input file that contains histograms and ntuples
  // made by ReadBack.cc
  TFile* file = new TFile( basename + ".root");

  // Name of the output postscript file.
  // Postscript/pdf is the only graphics format for which root supports multi-page output files.
  TString psfile( basename + ".pdf");

  // Get pointers to some of the histograms.
  TH1F* hMultiplicity;   file->GetObject("checkhits/hMultiplicity",  hMultiplicity);
  TH1F* hHitNeighbours;  file->GetObject("checkhits/hHitNeighbours", hHitNeighbours);
  TH1F* hEnergyDep;      file->GetObject("checkhits/hEnergyDep",     hEnergyDep);
  TH1F* hStepLength;     file->GetObject("checkhits/hStepLength",    hStepLength);

  TH1F* hRadius;         file->GetObject("checkhits/hRadius", hRadius);
  TH1F* hxHit;           file->GetObject("checkhits/hxHit",   hxHit);
  TH1F* hyHit;           file->GetObject("checkhits/hyHit",   hyHit);
  TH1F* hzHit;           file->GetObject("checkhits/hzHit",   hzHit);

  TH1F* hHitDeltaTime;   file->GetObject("readStrawHits/hHitDeltaTime",  hHitDeltaTime);
  TH1F* hHitEnergy;      file->GetObject("readStrawHits/hHitEnergy",     hHitEnergy);
  TH1F* hNHits;          file->GetObject("readStrawHits/hNHits",         hNHits);
  TH1F* hNHitsPerWire;   file->GetObject("readStrawHits/hNHitsPerWire",  hNHitsPerWire);

  TH1F* hNG4Steps;       file->GetObject("readStrawHits/hNG4Steps",      hNG4Steps);
  TH1F* hG4StepLength;   file->GetObject("readStrawHits/hG4StepLength",  hG4StepLength);
  TH1F* hG4StepEdep;     file->GetObject("readStrawHits/hG4StepEdep",    hG4StepEdep);
  TH1F* hG4StepRelTimes; file->GetObject("readStrawHits/hG4StepRelTimes",hG4StepRelTimes);

  // Get  pointers to the ntuples
  TNtuple* nt;  file->GetObject("checkhits/ntup",nt);
  TNtuple* snt; file->GetObject("readStrawHits/ntup",snt);

  // Open a new canvas on the screen.
  // The last two arguments are the size of the window.
  TCanvas *canvas = new TCanvas("c", "Plots from " + basename, 900, 900 );

  // Open a multi-page output postscript file .
  canvas->Print( psfile+"[");

  // Split the canvas into 4 pads.
  canvas->Divide(2,2);

  // Draw some histograms, one per pad.
  // cd(n): move to graphics pad number "n".
  // "H9": draw outline histogram ("H") in high resolution mode (9)
  canvas->cd(1); gPad->SetLogy(1); hMultiplicity->Draw("H9");
  canvas->cd(2); hHitNeighbours->Draw("H9");
  canvas->cd(3); gPad->SetLogy(1); hEnergyDep->Draw("H9");
  canvas->cd(4); gPad->SetLogy(1); hStepLength->Draw("H9");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a png file (slow)
  //canvas->Print( basename + "_1.png" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->Divide(2,2);

  canvas->cd(1); hRadius->Draw("H9");
  canvas->cd(2); hxHit->Draw("H9");
  canvas->cd(3); hyHit->Draw("H9");
  canvas->cd(4); hzHit->Draw("H9");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->Divide(2,2);

  canvas->cd(1); hHitDeltaTime->Draw("H9");
  canvas->cd(2); gPad->SetLogy(1); hHitEnergy->Draw("H9");
  canvas->cd(3); gPad->SetLogy(1); hNHits->Draw("H9");
  canvas->cd(4); gPad->SetLogy(1); hNHitsPerWire->Draw("H9");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->Divide(2,2);

  canvas->cd(1); gPad->SetLogy(1); hNG4Steps->Draw("H9");
  canvas->cd(2); gPad->SetLogy(1); hG4StepLength->Draw("H9");
  canvas->cd(3); gPad->SetLogy(1); hG4StepEdep->Draw("H9");
  canvas->cd(4); gPad->SetLogy(1); hG4StepRelTimes->Draw("H9");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // we specify the histogram names and bin info explicitly

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->Divide(2,2);

  canvas->cd(1); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhg1(500,300.,800.)","","");
  canvas->cd(2); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhg2(500,300.,800.)","genId==2","");
  canvas->cd(3); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhg3(500,300.,800.)","genId==0","");
  canvas->cd(4); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhg4(500,300.,800.)","genId!=0","");
  // canvas->cd(5); nt->Draw( "sqrt(hy*hy+hx*hx)>>mh5(500,300.,800.)","genId==12","");
  // canvas->cd(6); nt->Draw( "sqrt(hy*hy+hx*hx)>>mh6(500,300.,800.)","genId==20","");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->cd(1); snt->Draw( "sqrt(hity*hity+hitx*hitx)>>msh1(500,300.,800.)","","");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  canvas->Divide(2,2);

  canvas->cd(1); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhs1(500,300.,800.)","","");
  //  mh1->SetTitle("sqrt(hy*hy+hx*hx)");
  canvas->cd(2); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhs2(500,300.,800.)","trk==1","");
  //  mh2->SetTitle("sqrt(hy*hy+hx*hx),trk==1");
  canvas->cd(3); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhs3(500,300.,800.)","trk!=1","");
  //  mh3->SetTitle("sqrt(hy*hy+hx*hx),trk!=1");
  canvas->cd(4); nt->Draw( "sqrt(hy*hy+hx*hx)>>mhs4(500,300.,800.)","trk!=1&&abs(pdgId)<100","");
  //  mh4->SetTitle("sqrt(hy*hy+hx*hx),trk!=1&&abs(pdgId)<100");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  // Draw a y vs x scatter plot of the hit positions.
  TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
  frame->SetTitle("StepPoint y vs. x (mm)");
  nt->Draw( "hx:hy","","PSAME");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  // Draw a y vs x scatter plot of the hit positions.
  TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
  frame->SetTitle("StrawHit y vs. x (mm)");
  snt->Draw( "hitx:hity","","PSAME");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

  // Clear canvas in preparation for next page
  canvas->cd(0);
  canvas->Clear();

  // Draw a y vs x scatter plot of the hit positions.
  TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
  frame->SetTitle("StepPoint y vs. x (mm) genId==2");
  nt->Draw( "hx:hy","genId==2","PSAME");

  // Flush page to screen
  canvas->Update();

  // Add this canvas to the postscript file.
  canvas->Print(psfile);

  // Uncomment this line to save this canvas as a jpg file.
  // canvas->Print( basename + "_2.jpg" );

  // Prompt and wait for response before continuing.
  cerr << "Double click in the last active pad to continue: " ;
  gPad->WaitPrimitive();
  cerr << endl;

//   // Clear canvas in preparation for next page
//   canvas->cd(0);
//   canvas->Clear();

//   // Draw a y vs x scatter plot of the hit positions.
//   TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
//   frame->SetTitle("StepPoint y vs. x (mm) edep==0.");
//   nt->Draw( "hx:hy","edep==0.","PSAME");

//   // this works for simple histograms, not ntuple TGraphs ("scatterplots")
//   // TPaveStats *pave = (TPaveStats*)nt->GetListOfFunctions()->FindObject("stats");
//   // pave->SetOptStat(1111111);

//   // Flush page to screen
//   canvas->Update();

//   // Add this canvas to the postscript file.
//   canvas->Print(psfile);

//   // Uncomment this line to save this canvas as a jpg file.
//   // canvas->Print( basename + "_2.jpg" );

//   // Prompt and wait for response before continuing.
//   cerr << "Double click in the last active pad to continue: " ;
//   gPad->WaitPrimitive();
//   cerr << endl;

//   // Clear canvas in preparation for next page
//   canvas->cd(0);
//   canvas->Clear();

//   // Draw a y vs x scatter plot of the hit positions.
//   TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
//   frame->SetTitle("StepPoint y vs. x (mm) abs(pdgId)<100");
//   nt->Draw( "hx:hy","abs(pdgId)<100","PSAME");

//   // Flush page to screen
//   canvas->Update();

//   // Add this canvas to the postscript file.
//   canvas->Print(psfile);

//   // Uncomment this line to save this canvas as a jpg file.
//   // canvas->Print( basename + "_2.jpg" );

//   // Prompt and wait for response before continuing.
//   cerr << "Double click in the last active pad to continue: " ;
//   gPad->WaitPrimitive();
//   cerr << endl;

//   // Clear canvas in preparation for next page
//   canvas->cd(0);
//   canvas->Clear();

//   // Draw a y vs x scatter plot of the hit positions.
//   TH1F* frame = canvas->DrawFrame(-800., -800., 800., 800.);
//   frame->SetTitle("StepPoint y vs. x (mm) abs(pdgId)>100");
//   nt->Draw( "hx:hy","abs(pdgId)>100","PSAME");

//   // Flush page to screen
//   canvas->Update();

//   // Add this canvas to the postscript file.
//   canvas->Print(psfile);

//   // Uncomment this line to save this canvas as a jpg file.
//   // canvas->Print( basename + "_2.jpg" );

//   // Prompt and wait for response before continuing.
//   cerr << "Double click in the last active pad to continue: " ;
//   gPad->WaitPrimitive();
//   cerr << endl;

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
  // canvas->Print( basename + "_3.gif" );
  // canvas->Print( basename + "_3.pdf" );

  // Close the postscript file.
  canvas->Print(psfile+"]");

}

// This tells emacs to view this file in c++ mode.
// Local Variables:
// mode:c++
// End:
