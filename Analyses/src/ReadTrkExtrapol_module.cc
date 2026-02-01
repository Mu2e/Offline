//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
// data
#include "Offline/RecoDataProducts/inc/TrkToCaloExtrapol.hh"

// Other includes.
#include "cetlib_except/exception.h"

// Mu2e includes.
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

//root includes
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "cetlib/pow.h"


using namespace std;

namespace mu2e {

  static int ncalls(0);

  class ReadTrkExtrapol : public art::EDAnalyzer {
  public:
    explicit ReadTrkExtrapol(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel")),
      _trkToCaloExtrapolModuleLabel(pset.get<std::string>("trkToCaloExtrapolModuleLabel")),
      _Ntup(0){}

    virtual ~ReadTrkExtrapol() {
    }
    void beginJob();
    void endJob() {}

    void analyze(art::Event const& e );

  private:

    void readTrackExtrapolation(art::Event const& evt, bool skip);

    // Diagnostic level
    int _diagLevel;

    // Label of the extrapolated impact points
    std::string _trkToCaloExtrapolModuleLabel;

    bool _skipEvent;

    TTree* _Ntup;//Ntupla which contains informations about the extrapolation starting from MC

    Int_t _evt,//event Id
      _ntrks;

    Float_t _caloSec[100],
      _trkId[100],
      _trkTime[100],
      _trkTimeErr[100],
      _trkPathLenghtIn[100],
      _trkPathLenghtInErr[100],
      _trkPathLenghtOut[100],
      _trkPathLenghtOutErr[100],
      _trkMom[100],
      _trkMomX[100],
      _trkMomY[100],
      _trkMomZ[100],
      _trkPosX[100],
      _trkPosY[100],
      _trkPosZ[100];

  };



  void ReadTrkExtrapol::beginJob( ) {

    cout << "start ReadTrkExtrapol..."<<endl;

  }



  void ReadTrkExtrapol::analyze(art::Event const& evt ) {

    ++ncalls;

    if (ncalls == 1) {

      art::ServiceHandle<art::TFileService> tfs;
      _Ntup        = tfs->make<TTree>("trkExt", "Extrapolated track info");

      _Ntup->Branch("evt"    , &_evt ,   "evt/F");
      _Ntup->Branch("ntrks"  , &_ntrks , "ntrks/I");
      _Ntup->Branch("caloSec", &_caloSec, "caloSec[ntrks]/F");
      _Ntup->Branch("trkId", &_trkId, "trkId[ntrks]/F");
      _Ntup->Branch("trkTime", &_trkTime, "trkTime[ntrks]/F");
      _Ntup->Branch("trkTimeErr", &_trkTimeErr, "trkTimeErr[ntrks]/F");
      _Ntup->Branch("trkPathLenghtIn", &_trkPathLenghtIn, "trkPathLenghtIn[ntrks]/F");
      _Ntup->Branch("trkPathLenghtInErr", &_trkPathLenghtInErr, "trkPathLenghtInErr[ntrks]/F");
      _Ntup->Branch("trkPathLenghtOut", &_trkPathLenghtOut, "trkPathLenghtOut[ntrks]/F");
      _Ntup->Branch("trkPathLenghtOutErr", &_trkPathLenghtOutErr, "trkPathLenghtOutErr[ntrks]/F");
      _Ntup->Branch("trkMom", &_trkMom, "trkMom[ntrks]/F");
      _Ntup->Branch("trkMomX", &_trkMomX, "trkMomX[ntrks]/F");
      _Ntup->Branch("trkMomY", &_trkMomY, "trkMomY[ntrks]/F");
      _Ntup->Branch("trkMomZ", &_trkMomZ, "trkMomZ[ntrks]/F");
      _Ntup->Branch("trkPosX", &_trkPosX, "trkPosX[ntrks]/F");
      _Ntup->Branch("trkPosY", &_trkPosY, "trkPosY[ntrks]/F");
      _Ntup->Branch("trkPosZ", &_trkPosZ, "trkPosZ[ntrks]/F");

    }


    readTrackExtrapolation(evt, _skipEvent);

  } // end of analyze


  void ReadTrkExtrapol::readTrackExtrapolation(art::Event const& evt, bool skip){

    art::Handle<TrkToCaloExtrapolCollection>  handle;
    evt.getByLabel(_trkToCaloExtrapolModuleLabel, handle);
    const TrkToCaloExtrapolCollection* coll;
    const TrkToCaloExtrapol *trkExt;

    if (handle.isValid()) {
      coll = handle.product();
    } else {
      printf(">>> ERROR in readTrackExtrapolation::doExtrapolation: failed to locate collection");
      printf(". BAIL OUT. \n");
      return;
    }

    _evt   = evt.id().event();
    _ntrks = coll->size();

    for(int i=0; i<_ntrks; ++i){
      trkExt = &coll->at(i);
//       KalRepPtr const& trkPtr = trkExt->trk();
//       const KalRep* trk = trkPtr.get();

      _caloSec[i]    = trkExt->diskId();
      _trkId[i]      = trkExt->trackNumber();
      _trkTime[i]    = trkExt->time();
      _trkTimeErr[i] = trkExt->timeErr();

      _trkPathLenghtIn[i] = trkExt->pathLengthEntrance();
      _trkPathLenghtInErr[i] = trkExt->pathLenghtEntranceErr();
      _trkPathLenghtOut[i] = trkExt->pathLengthExit();
      _trkPathLenghtOutErr[i] = trkExt->pathLenghtEntranceErr();
      _trkMom[i] = trkExt->momentum().mag();
      _trkMomX[i] = trkExt->momentum().x();
      _trkMomY[i] = trkExt->momentum().y();
      _trkMomZ[i] = trkExt->momentum().z();
      _trkPosX[i] = trkExt->entrancePosition().x();
      _trkPosY[i] = trkExt->entrancePosition().y();
      _trkPosZ[i] = trkExt->entrancePosition().z();

    }

    _Ntup->Fill();

    if(evt.id().event() % 100 == 0){
      cout << "Event "<<evt.id().event()<<" ReadTrkExtrapol done..."<<endl;
    }
  }



}



using mu2e::ReadTrkExtrapol;
DEFINE_ART_MODULE(ReadTrkExtrapol)
