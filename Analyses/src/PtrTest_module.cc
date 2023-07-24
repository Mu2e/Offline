//
// Example of how to read a KalRepCollection and create a matching KalRepPtrCollection.
// This can be used as a model for writing a module that merges two KalRepCollections.
//
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

// Need this for the BaBar headers.
using namespace CLHEP;

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

// mu2e tracking
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

// C++ includes.
#include <iostream>
#include <string>

// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

using namespace std;

namespace mu2e {

  class PtrTest : public art::EDProducer {
  public:

    explicit PtrTest(fhicl::ParameterSet const& pset);
    void produce( art::Event& e) override;

  private:

    // Information about the data product that contains the fitted tracks.
    // The instance name is also used as the instance name of the output data product.
    TrkParticle     _tpart;
    TrkFitDirection _fdir;
    std::string     _instanceName;
    art::InputTag   _inputTag;

  };

}  // end namespace mu2e

mu2e::PtrTest::PtrTest(fhicl::ParameterSet const& pset):
  EDProducer{pset},
  _tpart((TrkParticle::type)(pset.get<int>("fitparticle"))),
  _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
  _instanceName( _fdir.name() + _tpart.name()),
  _inputTag( pset.get<string>("inputModuleLabel"), _instanceName)
{
  produces<KalRepPtrCollection>(_instanceName);
}

void mu2e::PtrTest::produce(art::Event& event)
{

  // Access the input data product.
  auto trksHandle  = event.getValidHandle<KalRepCollection>(_inputTag);
  auto const& inputTracks = *trksHandle;

  // Create an empty output data product.
  unique_ptr<KalRepPtrCollection> outputTracks(new KalRepPtrCollection );
  outputTracks->reserve(inputTracks.size());

  // Populate the output data product.
  for ( size_t i=0; i<inputTracks.size(); ++i ){

    KalRep const* krep = inputTracks.get(i);
    if ( !krep ) continue;

    outputTracks->emplace_back(trksHandle,i);
  }

  // Add the output data product to the event.
  event.put(std::move(outputTracks),_instanceName);

}

DEFINE_ART_MODULE(mu2e::PtrTest)
