//
//  Report the type of pattern recognition used
//  reconstruct a given track.
//

#include <exception>
#include <memory>
#include <string>
#include <typeinfo>

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"

#include "Mu2eUtilities/inc/decodeTrackPatRecType.hh"

namespace mu2e {

  TrackPatRecType decodeTrackPatRecType( KalRepPtr const& ptr, art::Event const& event ){
    art::Handle<KalRepCollection> handle;
    event.get(ptr.id(), handle);
    fhicl::ParameterSet const& pset = handle.provenance()->parameterSet();
    return TrackPatRecType(pset.get<std::string>("module_type"));
  }
}

