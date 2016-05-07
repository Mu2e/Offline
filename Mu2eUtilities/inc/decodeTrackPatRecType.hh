#ifndef Mu2eUtilities_decodeTrackPatRecType_hh
#define Mu2eUtilities_decodeTrackPatRecType_hh
//
//  Report the type of pattern recognition used
//  reconstruct a given track.
//
//  If the only information you need is this,
//  prefer this function to KalRepCollectionInfo
//  which looks for more information and is slower.
//

#include "Mu2eUtilities/inc/TrackPatRecType.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  TrackPatRecType decodeTrackPatRecType( KalRepPtr const& ptr, art::Event const& event );

}

#endif /* Mu2eUtilities_decodeTrackPatRecType_hh */
