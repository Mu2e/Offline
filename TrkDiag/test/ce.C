// Create two standard pdf files with tracking diagnostic plots
// Version for use with root v6.
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
// 1) In the c'tor of TChain, the TDirectory element of the string
//    is the module label of the module that ran ReadKalFits;
//    edit if necessary.
// 2) Edit the mychain->Add("filename") line to specify the right
//    input filename.
// 3) Add additional mychain->Add lines as needed.
// 4) If desired, edit the names of the output files
// 5) To run: root -l ce>C
//
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("emruo");

  TChain* mychain = new TChain("RKFDeM/trkdiag");
  mychain->Add("genReco.hist");

#include "TrkDiag/test/KalFit.C+"

  TFile * outFile = new TFile("save_ce.root","NEW");

  KalFit fit(mychain);
  fit.Acc();
  fit.Res();

  fit.acan->Print("acan_ce.pdf");
  fit.rcan->Print("rcan_ce.pdf");

  outFile->cd();
  outFile->Write();
  outFile->Close();
}
