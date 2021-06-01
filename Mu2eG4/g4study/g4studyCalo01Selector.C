#define g4studyCalo01Selector_cxx

//
// Original author K. Genser (of the content after TTree::MakeSelector() was used)
//

// the macro is NOT meant to be run using PROOF due to the way the
// TEDS map is handled; see below

// the original tree/selector were generated using
// TTree*   tstree = 0x0;
// _file0->GetObject("checkhits/tstep", tstree);
// tstree->MakeSelector("g4studyCalo01Selector");
// and the input file file was generated using Mu2eG4/test/g4studyCalo_01.fcl

// The class definition in g4studyCalo01Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("g4studyCalo01Selector.C")
// Root > T->Process("g4studyCalo01Selector.C","some options")
// Root > T->Process("g4studyCalo01Selector.C+")
//

#include "g4studyCalo01Selector.h"
#include <TH2.h>
#include <TStyle.h>
#include "TCanvas.h"

#include <iostream>
#include <sstream>

#include <vector>
#include <map>

using namespace std;

// data used in several functions

Int_t debug = 1;

bool includeOverflows = true;

Long_t maxPrintEvent = 5;
ULong_t maxPrintHit = 1;

Long_t nEvent;
Long_t fEvent;

vector<TCanvas*> canvases;

Long_t canvasox;
Long_t canvasoy;

bool     histogramsBooked(false);

Int_t    scMaxVol = 200;
Int_t    scMinVol =   2;

TH2D*    h_TED_vs_Vol;

TH2D*    h_TEDS_vs_Vol;

ULong_t  nBinTED = 500;
Double_t lBinTED = 0.0;
Double_t hBinTED = 5.0;
Double_t hBinTEDS = 50.0;

ULong_t  nBinVol = 100;
Double_t lBinVol = 100.0;
Double_t hBinVol = 200.0;

map<Int_t,Double_t> TEDS;
map<Int_t,Double_t>::const_iterator TEDSCIterator;

template <class myTH> void book(myTH* &hp, char const* const htitle, 
	      Int_t nBinX, Double_t lBinX, Double_t hBinX,
	      Int_t nBinY, Double_t lBinY, Double_t hBinY,
	      Color_t const color=kRed) {

  TString hname  = "h_";
  hname += htitle;
  hp =  new myTH(hname, htitle,
	     nBinX, lBinX, hBinX,
	     nBinY, lBinY, hBinY);

  hp->SetMarkerColor(color);
  hp->SetLineColor(color);  

  hp->Sumw2();
  hp->StatOverflows(includeOverflows);
  return;

}

template <class myTH> void book(myTH* &hp, char const* const htitle, 
				Int_t nBinX, Double_t lBinX, Double_t hBinX,
				Int_t nBinY, const Double_t* ybins,
				Color_t const color=kRed) {

  TString hname  = "h_";
  hname += htitle;
  hp =  new myTH(hname, htitle,
	     nBinX, lBinX, hBinX,
	     nBinY, ybins);

  hp->SetMarkerColor(color);
  hp->SetLineColor(color);  

  hp->Sumw2();
  hp->StatOverflows(includeOverflows);
  return;

}

template <class myTH> void book(myTH* &hp, char const* const htitle, 
	      Int_t nBinX, Double_t lBinX, Double_t hBinX,
	      Color_t const color=kRed) {

  TString hname  = "h_";
  hname += htitle;
  hp =  new myTH(hname, htitle,
		 nBinX, lBinX, hBinX);
  hp->SetMarkerColor(color);
  hp->SetLineColor(color);  
  hp->Sumw2();
  hp->StatOverflows(includeOverflows);
  return;
}


void g4studyCalo01Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void g4studyCalo01Selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   // reset the counter
   nEvent = -1;

}

Bool_t g4studyCalo01Selector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either g4studyCalo01Selector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  if (debug>3) {
    TString const option = GetOption();
    cout << "g4studyCalo01Selector::Process called with option " << option.Data() << endl;
  }

  if(!histogramsBooked) {

    histogramsBooked = true;

    cout << "g4studyCalo01Selector::Process booking histograms " << endl;

    book(h_TED_vs_Vol, "TED_vs_Vol",
         nBinVol, lBinVol, hBinVol,
         nBinTED, lBinTED, hBinTED);

    book(h_TEDS_vs_Vol, "TEDS_vs_Vol",
         nBinVol, lBinVol, hBinVol,
         nBinTED, lBinTED, hBinTEDS);

    cout << "g4studyCalo01Selector::Process GetEntries() " << fChain->GetEntries() 
         << endl;

  }

  b_evt  ->GetEntry(entry);
  b_vol  ->GetEntry(entry); 
  b_tedep->GetEntry(entry);

  if (debug>3) {
    cout << "g4studyCalo01Selector::Process evt, vol, tedep " 
         << evt << ", "
         << vol << ", "
         << tedep
         << endl;
  }


  // histogram per hit in the scintillator layers
  if ( vol < scMaxVol && vol >= scMinVol ) h_TED_vs_Vol->Fill(vol,tedep);

  if (nEvent != evt) {

    if (debug>2) {
      cout << "g4studyCalo01Selector::Process evt, nEvent " 
           << evt    << ", "
           << nEvent 
           << endl;
    }

    nEvent = evt;

    // fill the histogram with the collected data;

    // for the very first event the loop will not be entered as the map is empty

    for ( TEDSCIterator = TEDS.begin(); TEDSCIterator!=TEDS.end(); ++TEDSCIterator ) {
        h_TEDS_vs_Vol->Fill(TEDSCIterator->first, TEDSCIterator->second);
    }
    // empty the map
    TEDS.clear();

  } else {

    // collect the data into the map
    if ( vol < scMaxVol && vol >= scMinVol ) {
      TEDS[vol] += tedep;
    }

  }

  return kTRUE;
}

void g4studyCalo01Selector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  // fill the histogram with the collected data for the last event;
  for ( TEDSCIterator = TEDS.begin(); TEDSCIterator!=TEDS.end(); ++TEDSCIterator ) {
    h_TEDS_vs_Vol->Fill(TEDSCIterator->first, TEDSCIterator->second);
  }

}

TString const canvasNamePrefix = "g4studyCalo01";
TString canvasNameFullPrefix;

TString canvasName;
ostringstream forCanvasSuffix;

TCanvas *theCanvas;

ULong_t canvasNumber;

Long_t const canvasOriginX = 10;
Long_t const canvasOriginY = 10;

Long_t const canvasWX = 600;
Long_t const canvasWY = 600;

Long_t const canvasSpace  =  35;
Long_t const canvasShiftX =  25;
Long_t const canvasShiftY =  25;


void prepareNextCanvas( Int_t const logx = 0, Int_t const logy = 0,
			Int_t const gridx = 1, Int_t const gridy = 1
			) {

  Int_t canvasPosID = canvases.size()+1;
  ++canvasNumber;

  forCanvasSuffix.str("");
  forCanvasSuffix.width(3);
  forCanvasSuffix.fill('0');
  forCanvasSuffix << canvasNumber;

  cout << "constructing canvas " << forCanvasSuffix.str() << ", position id: " <<  canvasPosID << endl;

  canvasName = canvasNameFullPrefix + "_c" + forCanvasSuffix.str();
  
  canvasox = canvasOriginX + (canvasPosID%2)*(canvasWX+canvasSpace) + (canvasPosID/2)*canvasShiftX;
  canvasoy = canvasOriginY + (canvasPosID/2)*canvasShiftY;

  theCanvas = new TCanvas(canvasName, canvasName, canvasox, canvasoy, canvasWX, canvasWY);
  canvases.push_back(theCanvas);

  theCanvas->Divide(1,1);

  ULong_t const canvasSlotNumber = 1;

  theCanvas->cd(canvasSlotNumber);

  gPad->SetLogx(logx);
  gPad->SetLogy(logy);

  gPad->SetGridx(gridx);
  gPad->SetGridy(gridy);

}

void displayPrintCanvas() {
  theCanvas->Update();
  theCanvas->Print(canvasName+".png","png");
  //  theCanvas->Write();
  return;
}

void prepareDisplayPrintCanvas(TH1* const hist, char const* const opt="", 
			       Int_t const logx = 0, Int_t const logy = 0,
			       Int_t const gridx = 1, Int_t const gridy = 1
			       ) {

  prepareNextCanvas(logx,logy,gridx,gridy);

  cout << "drawing  " << hist->GetTitle() <<  endl;

  gStyle->SetOptStat("neMRuoi");
  hist->Draw(opt);

  displayPrintCanvas();

}

void prepareDisplay2PrintCanvas(TH1* const hist1, 
				TH1* const hist2, 
				char const* const opt, 
				Int_t const logx = 0, Int_t const logy = 0,
				Int_t const gridx = 1, Int_t const gridy = 1
				) {

  prepareNextCanvas(logx,logy,gridx,gridy);

  cout << "drawing  " << hist1->GetTitle() <<  endl;


  gStyle->SetOptStat("neMRuoi");
  hist1->SetStats(kTRUE);
  hist1->Draw(opt);

  string h2opt("SAME");

  h2opt += opt;

  cout << "drawing  " << hist2->GetTitle() <<  endl;

  hist2->SetStats(kFALSE);
  hist2->Draw(h2opt.c_str());

  displayPrintCanvas();

}

void g4studyCalo01Selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  TString const option = GetOption();
  cout << "g4studyCalo01Selector::Terminate called with option " << option.Data() << endl;

  canvasNameFullPrefix = canvasNamePrefix + "_" + option;

  canvasNumber = 0;

  prepareDisplayPrintCanvas(h_TED_vs_Vol,"");
  prepareDisplayPrintCanvas(h_TEDS_vs_Vol,"");

}
