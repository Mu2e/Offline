//
// c++ (not cint) Root "script" to make some plots based on a root example
// and ReadBack.cc
//
// $Id: Analyzer.C,v 1.12 2011/10/28 18:47:07 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:07 $
//
// Original author KLG
//


#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include <Cintex/Cintex.h>

#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Wrapper.h"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "LTrackerGeom/inc/LTracker.hh"

#include "Analyzer.h"

using namespace std;

// do not make it inline so that root can see the symbols
Analyzer::Analyzer (char const * file,
                    ULong64_t maxevent,
                    ULong64_t maxFullPrint,
                    Double_t minEnergy,
                    char const * cformat):
  _file(file),
  _mevent(maxevent),
  //  _g4ModuleLabel(g4ModuleLabel),
  _minimumEnergy(minEnergy),
  _maxFullPrint(maxFullPrint),
  _nAnalyzed(0),
  _outputFileNamePrefix("Analyzer"),
  _canvasPrintFormat(cformat),
  _hRadius(0),
  _hEnergyDep(0),
  _hTime(0),
  _hMultiplicity(0),
  _hDriftDist(0),
  _hxHit(0),
  _hyHit(0),
  _hzHit(0),
  _hHitNeighbours(0),
  _hCheckPointRadius(0),
  _hMomentumG4(0),
  _hStepLength(0),
  _ntup(0),
  _canvases(0)
{}

// we want to keep the objects after the script exits
// Analyzer::~Analyzer() {
//   delete _hRadius;
//   delete _hEnergyDep;
//   delete _hTime;
//   delete _hMultiplicity;
//   delete _hDriftDist;
//   delete _hxHit;
//   delete _hyHit;
//   delete _hzHit;
//   delete _hHitNeighbours;
//   delete _hCheckPointRadius;
//   delete _hMomentumG4;
//   delete _hStepLength;
//   delete _ntup;
//   delete _canvases;
//   delete _histograms;
// }

void Analyzer::begin(){

  gROOT->SetStyle("Plain");
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerColor(kRed);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetHistLineColor(kRed);
  //  gStyle->SetFuncWidth(1);
  //  gStyle->SetOptFit(1111);
  //gStyle->SetOptStat(111111);
  //gStyle->SetOptStat("neMRuoiSK");

  // Create some 1D histograms.
  _hRadius       = new TH1F("hRadius",       "Radius of Hits;(mm)",       100,  0., 1000. );
  _hEnergyDep    = new TH1F("hEnergyDep",    "Energy Deposited;(keV)",    100,  0.,   10. );
  _hTime         = new TH1F("hTime",         "Pulse Height;(ns)",         100,  0., 2000. );
  _hMultiplicity = new TH1F("hMultiplicity", "Hits per Event",            100,  0.,  100. );
  _hDriftDist    = new TH1F("hDriftDist",    "Crude Drift Distance;(mm)", 100,  0.,    3. );

  _hxHit         = new TH1F("hxHit",  "X of Hit;(mm)",                    100, -1000., 1000. );
  _hyHit         = new TH1F("hyHit",  "Y of Hit;(mm)",                    100, -1000., 1000. );
  _hzHit         = new TH1F("hzHit",  "Z of Hit;(mm)",                    100, -1400., 1400. );

  _hHitNeighbours    = new TH1F("hHitNeighbours",  "Number of hit neighbours",         10,   0., 10. );
  _hCheckPointRadius = new TH1F("hCheckPointRadius","Radius of Reference point;(mm)", 100, 2.25, 2.75 );

  _hMomentumG4 = new TH1F("hMomentumG4",  "Mommenta of particles created inside G4; (MeV)", 100, 0., 100. );
  _hStepLength = new TH1F("hStepLength",  "G4 Step Length in Sensitive Detector; (mm)",     100, 0.,  10. );

  // Create an ntuple.
  _ntup           = new TNtuple("ntup", "Hit ntuple",
                                "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec:pdgId:genId:edep:p:step");
  _canvases = new std::vector<TCanvas*>();
  _canvases->reserve(20);
  return;
}

void Analyzer::plotHist(TH1* hist, char const * opt) {

  // plot a single histogram
  TString canvasTitle =  hist->GetTitle();
  TCanvas* theCanvas = prepareNextCanvas();
  theCanvas->SetTitle(canvasTitle.Data());
  cout << "drawing: " << canvasTitle << " on canvas" << theCanvas->GetName() << endl;
  hist->Draw(opt);
  theCanvas->Update();

}

void Analyzer::plotNHist(vector<TH1*> vhist, char const * opt) {

  // an example on how to plot N histograms on one canvas

  ULong64_t nhist = vhist.size();
  // is nhist even?
  ULong64_t nhalf = nhist/2;
  ULong64_t nx = 2;
  ULong64_t ny = (nhalf*2 == nhist) ? nhalf : nhalf+1;

  TCanvas* theCanvas = prepareNextCanvas(nx,ny);
  TString canvasTitle;

  for (ULong64_t i=0; i<nhist; ++i) {

    ULong64_t ipad = i+1;
    theCanvas->cd(ipad);
    canvasTitle += vhist.at(i)->GetTitle();
    if (ipad<nhist) canvasTitle += ", ";
    cout << "drawing: " << vhist[i]->GetTitle() << " on canvas" <<
      theCanvas->GetName() << " on pad " << ipad << endl;
    vhist.at(i)->Draw(opt);

  }

  theCanvas->SetTitle(canvasTitle.Data());
  theCanvas->Update();

}

void Analyzer::plotNT(char const * nts, char const * cut, char const * opt) {

  TString canvasTitle(nts);
  TCanvas* theCanvas = prepareNextCanvas();
  theCanvas->SetTitle(canvasTitle.Data());
  cout << "drawing: " << canvasTitle << " on canvas" << theCanvas->GetName() << endl;
  // Long64_t TTree::Draw(const char* varexp, const TCut& selection, Option_t* option = "", ...
  _ntup->Draw(nts,cut,opt);
  theCanvas->Update();

}

void Analyzer::plot() {

//   plotHist( _hRadius);
//   plotHist( _hEnergyDep);
//   plotHist( _hTime);
  plotHist( _hMultiplicity);
  // plotHist( _hDriftDist);
//   plotHist( _hxHit);
//   plotHist( _hyHit);
//   plotHist( _hzHit);

  vector<TH1*> htp;
  htp.push_back(_hRadius);
  htp.push_back(_hEnergyDep);
  htp.push_back(_hTime);
  htp.push_back(_hxHit);
  htp.push_back(_hyHit);
  htp.push_back(_hzHit);

  plotNHist(htp);

  //  plotHist( _hHitNeighbours);
  //  plotHist( _hCheckPointRadius);
  plotHist( _hMomentumG4);
  plotHist( _hStepLength);

  plotNT("hx:hy");
  plotNT("hz:time");
  plotNT("hz:step");
  plotNT("sid:trk");
  plotNT("pdgId");
  plotNT("genId");
  plotNT("step:edep");

}

void Analyzer::write() {

  // store canvases/histograms/functions

  TFile*   outputRootFile;
  TString  outputRootFileName(_outputFileNamePrefix);

  outputRootFileName += ".root";

  outputRootFile = new TFile(outputRootFileName,"recreate");

  // get the current objects in memory, here: histogrmas and the nutuple and write them out

  gROOT->GetList()->Write();
  gROOT->GetListOfCanvases()->Write();

  gDirectory->ls("-m");
  gDirectory->ls("-d");

  outputRootFile->Close();

  return;

}

void Analyzer::printOutCanvases() {

  // write out individual canvase files (e.g. png)

  TIter citer(gROOT->GetListOfCanvases());
  while ( TCanvas* nc =  dynamic_cast<TCanvas*>(citer.Next()) ) {
    TString cn = nc->GetName();
    nc->Print(cn+"."+_canvasPrintFormat,_canvasPrintFormat);
  }

}

TCanvas* Analyzer::prepareNextCanvas( Int_t nx, Int_t ny,
                                      Int_t const logx, Int_t const logy,
                                      Int_t const gridx, Int_t const gridy
                                      ) {

  Int_t canvasPosID = _canvases->size()+1;
  Int_t canvasNumber = canvasPosID;

  ostringstream forCanvasSuffix("");
  //  forCanvasSuffix.str("");
  forCanvasSuffix.width(3);
  forCanvasSuffix.fill('0');
  forCanvasSuffix << canvasNumber;

  cout << "constructing canvas " << forCanvasSuffix.str() << ", position id: " <<  canvasPosID << endl;

  TString canvasName = _outputFileNamePrefix + "_c" + forCanvasSuffix.str();

  Long_t canvasox = canvasOriginX + (canvasPosID%2)*(canvasWX+canvasSpace) + (canvasPosID/2)*canvasShiftX;
  Long_t canvasoy = canvasOriginY + (canvasPosID/2)*canvasShiftY;

  TCanvas* theCanvas = new TCanvas(canvasName, canvasName, canvasox, canvasoy, canvasWX, canvasWY);
  _canvases->push_back(theCanvas); // we want to be able to access canvases after the script exits


  // void TPad::Divide(Int_t nx = 1, Int_t ny = 1, Float_t xmargin = 0.01, Float_t ymargin = 0.01, Int_t color = 0)
  // e.g see http://root.cern.ch/root/html526/TCanvas.html
  theCanvas->Divide(nx,ny);

  ULong64_t npads = nx*ny;

  for (ULong_t ipad=1; ipad<=npads; ++ipad) {

    theCanvas->cd(ipad);

    gPad->SetLogx(logx);
    gPad->SetLogy(logy);

    gPad->SetGridx(gridx);
    gPad->SetGridy(gridy);

  }

  return theCanvas;

}

void Analyzer::analyze() {

  // Open the input file that contains Mu2e event data.
  TFile* file = new TFile(Analyzer::_file);

  // Fill a pointer to the Events & Run tree.

  TString EventsTree_name("Events");
  TTree* EventsTree; file->GetObject(EventsTree_name,EventsTree);

  TString RunsTree_name("Runs");
  TTree* RunsTree; file->GetObject(RunsTree_name,RunsTree);

  // Br    3 :mu2eRandomEngineStates_randomsaver__G4Test03.obj :
  // Br    9 :mu2eSimParticles_g4run__G4Test03.obj :
  // Br   41 :mu2eStepPointMCs_g4run__G4Test03.obj :
  // Br   55 :mu2eGenParticles_generate__G4Test03.obj :

  // Br    3 :mu2ePhysicalVolumeInfos_g4run__G4Test03.obj :


  // Branch names without the trailing "." if any

  TString EventAuxiliaryBranchName    ("EventAuxiliary");                   //  the event info / Events tree

  TString SimParticleBranchName       ("mu2eSimParticles_g4run__G4Test03"); // simParticles / Events tree
  TString StepPointMCBranchName       ("mu2eStepPointMCs_g4run__G4Test03"); // hits / Events tree
  TString GenParticleBranchName       ("mu2eGenParticles_generate__G4Test03"); // genParticles / Events tree
  TString PhysicalVolumeInfoBranchName("mu2ePhysicalVolumeInfos_g4run__G4Test03"); // volumes / Runs tree

  // we should create a class/struct for all data needed for one wrapped type
  // here is the running list: BranchName, Wrppd,

  // art::Wrapper<mu2e::SimParticleCollection>   w;
  // art::Wrapper<mu2e::SimParticleCollection>* ww; ww = &w;
  // make sure to not to use cint for the code below, rely on a complier


  art::EventAuxiliary* EventAuxiliaryWrppd = new art::EventAuxiliary(); // this is a very different branch

  art::Wrapper<mu2e::SimParticleCollection>* SimParticleWrppd =
    new art::Wrapper<mu2e::SimParticleCollection>();
  art::Wrapper<mu2e::StepPointMCCollection>* StepPointMCWrppd =
    new art::Wrapper<mu2e::StepPointMCCollection>();
  art::Wrapper<mu2e::GenParticleCollection>* GenParticleWrppd =
    new art::Wrapper<mu2e::GenParticleCollection>();
  art::Wrapper<mu2e::PhysicalVolumeInfoCollection>* PhysicalVolumeInfoWrppd =
    new art::Wrapper<mu2e::PhysicalVolumeInfoCollection>();

  //disable branch ("*");

  EventsTree->SetBranchStatus("*",0); //disable all branches
  RunsTree->SetBranchStatus("*",0);   //disable all branches

  //enable subbranches, relies on the "." in the branch name

  EventsTree->SetBranchStatus(EventAuxiliaryBranchName+"*",1);
  EventsTree->SetBranchAddress(EventAuxiliaryBranchName,&EventAuxiliaryWrppd);

  EventsTree->SetBranchStatus(SimParticleBranchName+"*",1);
  EventsTree->SetBranchAddress(SimParticleBranchName+".",&SimParticleWrppd);
  EventsTree->SetBranchStatus(StepPointMCBranchName+"*",1);
  EventsTree->SetBranchAddress(StepPointMCBranchName+".",&StepPointMCWrppd);
  EventsTree->SetBranchStatus(GenParticleBranchName+"*",1);
  EventsTree->SetBranchAddress(GenParticleBranchName+".",&GenParticleWrppd);

  RunsTree->SetBranchStatus(PhysicalVolumeInfoBranchName+"*",1);
  RunsTree->SetBranchAddress(PhysicalVolumeInfoBranchName+".",&PhysicalVolumeInfoWrppd);

  ULong64_t EventsNEntries = EventsTree->GetEntries();
  ULong64_t RunsNEntries = RunsTree->GetEntries();

  cout << "EventsNEntries/RunsNEntries  " << EventsNEntries << "/" << RunsNEntries  << endl;
  // EventsNEntries/RunsNEntries  200/1

  cout << "EventAuxiliaryWrppd->id().event()  " << EventAuxiliaryWrppd->id().event() << endl;

  cout << "Size of SimParticle branch          " <<
    EventsTree->GetBranch(SimParticleBranchName+".")->GetEntries() << endl;

  ULong64_t EventsMaxEntries = min(EventsNEntries,Analyzer::_mevent);

  for ( ULong64_t i = 0; i<EventsMaxEntries; ++i) {
    EventsTree->GetEntry(i);
    // the object event has been filled at this point for the activated branches
    // Examples of elements of each object: (may not work if there are no
    // entries in a collection, fix it for it later)

    if (i<Analyzer::_maxFullPrint) {

      cout << "EventAuxiliaryWrppd->id().event()  " << EventAuxiliaryWrppd->id().event() << endl;

      key_type zero(0);

      cout << "Event i " << i << " SimParticle _endPosition.dx " <<
        SimParticleWrppd->product()->at(zero).endPosition().x() << endl;
      cout << "Event i " << i << " StepPointMC    _position.dx "
           << StepPointMCWrppd->product()->at(0).position().x() << endl;
      cout << "Event i " << i << " GenParticle _position.dx "
           << GenParticleWrppd->product()->at(0)._position.x() << endl;

      cout << "Size of StepPointMCWrppd->product() " << StepPointMCWrppd->product()->size() << endl;

    }

//     doLTracker(EventAuxiliaryWrppd,
//                StepPointMCWrppd,
//                GenParticleWrppd,
//                SimParticleWrppd,
//                PhysicalVolumeInfoWrppd);

// using those types in the root dict is problematic..., so we use them locally in this function

// void Analyzer::doLTracker(art::EventAuxiliary*                              EventAuxiliaryWrppd,
//                           art::Wrapper<mu2e::StepPointMCCollection>*        StepPointMCWrppd,
//                           art::Wrapper<mu2e::GenParticleCollection>*        GenParticleWrppd,
//                           art::Wrapper<mu2e::SimParticleCollection>*        SimParticleWrppd,
//                           art::Wrapper<mu2e::PhysicalVolumeInfoCollection>* PhysicalVolumeInfoWrppd
//                           ){

    // "aliases/handles"
    art::EventAuxiliary                const & event        = *EventAuxiliaryWrppd;
    mu2e::StepPointMCCollection        const * hits         = StepPointMCWrppd->product();
    mu2e::GenParticleCollection     const * genParticles    = GenParticleWrppd->product();
    mu2e::SimParticleCollection        const * simParticles = SimParticleWrppd->product();
    mu2e::PhysicalVolumeInfoCollection const * volumes      = PhysicalVolumeInfoWrppd->product();

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles!=0 && volumes!=0 );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // ntuple buffer.
    float nt[18];

    // Loop over all hits.
    UInt_t maxhit = hits->size();
    for ( UInt_t i=0; i<maxhit; ++i ){

      // Alias, used for readability.
      const mu2e::StepPointMC& hit = hits->at(i);

      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      //     // Get the straw information:
      //     const mu2e::Straw&      straw = tracker.getStraw( hit.strawIndex() );
      //     const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      //     const CLHEP::Hep3Vector& w     = straw.getDirection();

      //     // Count how many nearest neighbours are also hit.
      //     int nNeighbours = countHitNeighbours( straw, hits );

      //     // Compute an estimate of the drift distance.
      //     TwoLinePCA pca( mid, w, pos, mom);

      //     // Check that the radius of the reference point in the local
      //     // coordinates of the straw.  Should be 2.5 mm.
      //     double s = w.dot(pos-mid);
      //     CLHEP::Hep3Vector point = pos - (mid + s*w);

      // The simulated particle that made this hit.
      key_type trackId = hit.trackId();

      // Default values for these, in case information is not available.
      int pdgId(0);
      mu2e::GenId genId;

      if ( haveSimPart ){
        mu2e::SimParticle const& sim = simParticles->at(trackId);

        // PDG Particle Id of the sim particle that made this hit.
        pdgId = sim.pdgId();

        // If this is a generated particle, which generator did it come from?
        // This default constructs to "unknown".
        if ( sim.fromGenerator() ){
          mu2e::GenParticle const& gen = genParticles->at(sim.generatorIndex());
          genId = gen._generatorId;
        }
      }

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hEnergyDep->Fill(hit.eDep()/CLHEP::keV);
      _hTime->Fill(hit.time());
      // _hHitNeighbours->Fill(nNeighbours);
      // _hCheckPointRadius->Fill(point.mag());

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      //    _hDriftDist->Fill(pca.dca());

      _hStepLength->Fill( hit.stepLength() );

      // Fill the ntuple. (we comment out the elemets requiring the geometry service)
      nt[0]  = event.id().event();
      nt[1]  = hit.trackId().asInt();
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      // nt[6]  = mid.x();
      // nt[7]  = mid.y();
      // nt[8]  = mid.z();
      // nt[9]  = pca.dca();
      nt[10] = hit.time();
      // nt[11] = straw.Id().getDevice();
      // nt[12] = straw.Id().getSector();
      nt[13] = pdgId;
      nt[14] = genId.Id();
      nt[15] = hit.eDep()/CLHEP::keV;
      nt[16] = mom.mag();
      nt[17] = hit.stepLength();

      _ntup->Fill(nt);

      //     // Print out limited to the first few events.
      //     if ( _nAnalyzed < _maxFullPrint ){

      //       cerr << "Readback hit: "
      //            << event.id().event() << " "
      //            << i                  <<  " "
      //            << hit.trackId()      << "   "
      //            << hit.volumeId()     << " "
      //            << straw.Id()         << " | "
      //            << pca.dca()          << " "
      //            << pos                << " "
      //            << mom                << " "
      //            << point.mag()        << " "
      //            << hit.eDep()         << " "
      //            << s
      //            << endl;
      //     }

    } // end loop over hits.


    // Additional printout and histograms about the simulated particles.
    //   if ( haveSimPart && (_nAnalyzed < _maxFullPrint) ){

    //     ConditionsHandle<ParticleDataTable> pdt("ignored");

    for ( mu2e::SimParticleCollection::const_iterator i=simParticles->begin();
          i!=simParticles->end(); ++i ){

      mu2e::SimParticle const& sim = i->second;
      //    for ( size_t i=0; i<simParticles->size(); ++ i){
      //mu2e::SimParticle const& sim = simParticles->at(i);

      if ( sim.madeInG4() ) {

        _hMomentumG4->Fill( sim.startMomentum().rho() );

      }

    }

    //       } else {

    //         // Particle Data group Id number.
    //         int pdgId = sim.pdgId();

    //         // Information about generated particle.
    //         GenParticle const& gen = genParticles->at(sim.generatorIndex());
    //         GenId genId(gen._generatorId);

    //         // Physical volume in which this track started.
    //         PhysicalVolumeInfo const& volInfo = volumes->at(sim.startVolumeIndex());

    //         cerr << "Simulated Particle: "
    //              << i                   << " "
    //              << pdgId               << " "
    //              << genId.name()        << " "
    //              << sim.startPosition() << " "
    //              << volInfo.name()      << " "
    //              << volInfo.copyNo()
    //              << endl;
    //       }

    //     }
    //   }

  } // end event loop

} // end analyze



