//
//  Report the type of pattern recognition used
//  reconstruct a given track.
//

#include <exception>                          // for exception
#include <memory>                                    // for allocator, uniqu...
#include <string>                                    // for string
#include <typeinfo>                                  // for type_info

#include "Mu2eUtilities/inc/decodeTrackPatRecType.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"  // for KalRepCollection
#include "art/Framework/Principal/Handle.h"          // for Handle
#include "art/Framework/Principal/Event.h"           // for Event
#include "art/Framework/Principal/Provenance.h"      // for Provenance
#include "fhiclcpp/ParameterSet.h"                   // for ParameterSet
#include "fhiclcpp/exception.h"                      // for exception

namespace mu2e {

  TrackPatRecType decodeTrackPatRecType( KalRepPtr const& ptr, art::Event const& event ){
    art::Handle<KalRepCollection> handle;
    event.get(ptr.id(), handle);
    fhicl::ParameterSet const& pset = handle.provenance()->parameterSet();
    return TrackPatRecType(pset.get<std::string>("module_type"));
  }
}

