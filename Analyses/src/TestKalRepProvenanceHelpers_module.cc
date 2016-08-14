//
// Module used to develop, test and illustrate the use of some
// helper functions/classes to access provenance information about
// KalRepCollections.
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/decodeTrackPatRecType.hh"
#include "Mu2eUtilities/inc/KalRepCollectionInfo.hh"
#include "Mu2eUtilities/inc/TrackPatRecType.hh"
#include "Mu2eUtilities/inc/TrkSpecies.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class TestKalRepProvenanceHelpers : public art::EDAnalyzer {
  public:

    explicit TestKalRepProvenanceHelpers(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e) override;

  private:

    // Module label of the module that created a KalRepPtrCollection.
    art::InputTag _tracksTag;

  };
}  // end namespace mu2e

mu2e::TestKalRepProvenanceHelpers::TestKalRepProvenanceHelpers(fhicl::ParameterSet const& pset) :
  art::EDAnalyzer(pset),
  _tracksTag(pset.get<std::string>("tracksTag")){
  TrackPatRecType::printAll(cout);
}

void mu2e::TestKalRepProvenanceHelpers::analyze(const art::Event& event) {

  // Test the collection level helper.
  art::Handle<KalRepCollection> kalReps;
  if ( event.getByLabel(_tracksTag, kalReps) ) {
    KalRepCollectionInfo info( kalReps );
    cout << "KalReps: "
         << kalReps->size() << "  | "
         << info.patRecType() <<  " "
         << info.instanceName() <<  " "
         << info.direction().name() <<  " "
         << info.particleType() << " "
         << info.charge()
         << endl;
  }

  // Test the ptr level helpers.
  auto ptrs = event.getValidHandle<KalRepPtrCollection>(_tracksTag);
  cout << "KalRep ptrs: " << ptrs->size() << endl;
  for ( auto const& ptr : *ptrs ){
    KalRepCollectionInfo info( ptr, event);
    TrackPatRecType type = decodeTrackPatRecType( ptr, event);
    cout << "    : "
         << _tracksTag.label()       << " "
         << ptr                      << " "
         << info.patRecType()        << " "
         << type                     << " | "
         << info.instanceName()      << " "
         << info.direction().name()  << " "
         << info.particleType()      << " "
         << info.charge()
         << endl;
  }

} // end analyze

DEFINE_ART_MODULE(mu2e::TestKalRepProvenanceHelpers);
