//
// Example of accessing fitted tracks via a KalRepPtrCollection.
//
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"

// Need this for the BaBar headers.
using namespace CLHEP;

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"

// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
//#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "TH1D.h"

// C++ includes.
#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class ReadPtr : public art::EDAnalyzer {
  public:

    explicit ReadPtr(fhicl::ParameterSet const& pset);

    void analyze( art::Event const& e) override;
    void beginJob() override;

  private:

    // Information about the data product that contains the fitted tracks.
    TrkParticle     _tpart;
    TrkFitDirection _fdir;
    std::string     _instanceName;
    art::InputTag   _inputTag;

    TH1D* _hNTracks;
    TH1D* _hFitStatus;
    TH1D* _hP;

  };
}   // end namespace mu2e

mu2e::ReadPtr::ReadPtr(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  _tpart((TrkParticle::type)(pset.get<int>("fitparticle"))),
  _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
  _instanceName( _fdir.name() + _tpart.name()),
  _inputTag( pset.get<string>("inputModuleLabel"), _instanceName)
{}

void mu2e::ReadPtr::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  _hNTracks   = tfs->make<TH1D>( "hNTracks",   "Number of tracks per Event",  5,  0.,   5.);
  _hFitStatus = tfs->make<TH1D>( "hFitStatus", "Fit Status",                  10, 0.,  10.);
  _hP         = tfs->make<TH1D>( "hP",         "Reconstructed Momentum",     30, 80., 110.);
}

void mu2e::ReadPtr::analyze(art::Event const& event) {

  auto trksHandle    = event.getValidHandle<KalRepPtrCollection>(_inputTag);
  auto const& tracks = *trksHandle;

  _hNTracks->Fill(tracks.size());

  for ( auto const& trk : tracks ){

    KalRep const& krep = *trk;

    int stat = krep.fitStatus().success();
    _hFitStatus->Fill(stat);

    if ( stat != 1 ) continue;

    double s0   = krep.startValidRange();
    _hP->Fill( krep.momentum(s0).mag() ) ;


  }

}


DEFINE_ART_MODULE(mu2e::ReadPtr);
