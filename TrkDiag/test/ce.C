// Create two standard pdf files with tracking diagnostic plots
// Version for use with root v5.
//
// This code reads a TTree made by the module:
//    TrkDiag/src/ReadKalFits_module.cc
//
// For an example of using this module:
//    Analyses/test/genReco.fcl
//
// The output is two pdf files and a root file
// acan_ce.pdf  - Acceptance/Efficiency cut tree
//              - top plot is cumulative
//              - bottom plot is differential
// rcan_ce.pdf  - Tracker resolution plots
//              - for 4 different values of the trkqual cut
// save_ce.root - A root file containing the histograms and TFunctions
//                that were displayed in the pdf file.
//
// Instructions:
// 1) Edit the mychain->Add("filename") line to specify the right
//    input filename.
// 2) Add additional mychain->Add lines as needed.
// 3) If desired, edit the names of the output files
// 4) root -l ce>C
//
{

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("emruo");

  TChain* mychain = new TChain("RKFDownstreameMinus/trkdiag");
  mychain->Add("genReco.hist");

  gROOT->LoadMacro("TrkDiag/test/KalFit.C+");

  TFile* outFile = new TFile("save_ce.root","NEW");

  KalFit fit(mychain);
  fit.Acc();
  fit.Res();
  acan->Print("acan_ce.pdf");
  rcan->Print("rcan_ce.pdf");

  outFile->cd();
  outFile->Write();
  outFile->Close();
}
