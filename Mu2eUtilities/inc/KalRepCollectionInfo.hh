#ifndef Mu2eUtilities_KalRepCollectionInfo_hh
#define Mu2eUtilities_KalRepCollectionInfo_hh
//
// Convenience class to access provenance information about a KalRepCollection.
//
// If the only information you need is the TrackPatRecType from a KalRepPtr,
// then prefer to use Mu2eUtilities/inc/decodeTrackPatRecType.hh which runs
// faster.
//
// Contact person, Rob Kutschke
//

#include <string>

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "art/Framework/Principal/Event.h"

#include "Mu2eUtilities/inc/KalRepInstanceNameDecoder.hh"
#include "Mu2eUtilities/inc/TrackPatRecType.hh"
#include "Mu2eUtilities/inc/TrkSpecies.hh"

namespace art { class Event; }
namespace art { class ProductID; }
namespace art { template <typename T> class Handle; }
namespace art { template <typename T> class ValidHandle; }

namespace mu2e {

  class KalRepCollectionInfo {
  public:
    KalRepCollectionInfo( art::Handle<KalRepCollection> const&      handle );
    KalRepCollectionInfo( art::ValidHandle<KalRepCollection> const& handle );
    KalRepCollectionInfo( KalRepPtr const& ptr, art::Event const& event );
    KalRepCollectionInfo( art::ProductID const& id,    art::Event const& event );

    TrackPatRecType               patRecType()   const { return patRecType_; }
    TrkFitDirection               direction()    const { return instance_.direction();    }
    TrkSpecies                    particleType() const { return instance_.particleType(); }
    int                           charge()       const { return instance_.charge();       }
    std::string const&            instanceName() const { return instance_.instanceName(); }

  private:
    TrackPatRecType           patRecType_;
    KalRepInstanceNameDecoder instance_;

  };


}

#endif
